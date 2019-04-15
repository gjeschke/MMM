function diagnostics = rigi_flex_engine(restraints,options,handles)
% Compute rigid-body arrangements for a RigiFlex model
%
% function diagnostics = rigi_flex_engine(restraints,options,handles)
%
% restraints    restraint definition
% options       running options
% handles       handles to the RigiFlex GUI, can be missing for modelling
%               on a remote server
%
% G. Jeschke, 20.4.2017

global model

stemloop_mode = false;
trial_echo = false;


% restraints.stemlibs = {};
% restraints.stemloop_rlinks.rba = [];
% restraints.stemloop_rlinks.coor1 = [];
% restraints.stemloop_rlinks.coor2 = [];
% restraints.stemloop_rlinks.lib = [];

solutions_given = false;
solutions = [];

if isfield(restraints,'solutions') && ~isempty(restraints.solutions)
    solutions_given = true;
    poi = strfind(restraints.solutions,'.dat');
    if isempty(poi)
        soln_name = strcat(restraints.solutions,'.dat');
    else
        soln_name = restraints.solutions;
    end
    solutions = load(soln_name);
    solutions = round(solutions);
end

skip_mode = true;

if isfield(restraints,'solution_mode') && ~isempty(restraints.solution_mode)
    switch lower(restraints.solution_mode)
        case 'single'
            add_msg_board('Warning: RigiFlex performs only the trials specified in the solution file');
        case 'all'
            add_msg_board('Warning: RigiFlex performs only the trials specified in the solution file (several combinations possible)');
            skip_mode = false;
        case 'randomized'
            add_msg_board('Warning: RigiFlex only tests random stemloop combinations for specified RBA trials');
            solutions = solutions(:,1:2);
        case 'combisearch'
            add_msg_board('Warning: RigiFlex only searches all stemloop combinations for specified RBA trials');
            solutions = solutions(:,1:2);
            skip_mode = false;
        otherwise
            add_msg_board('Warning: Unspecified solution mode. RigiFlex performs only the trials specified in the solution file');
    end
end

tph = 1000000; % trials per hour, to be replaced with machine-specific value

target_resolution = 3; % resolution with spin labels is not realistically better than 3 Å
res = target_resolution;

diagnostics.success = 0;

% if handles of the GUI figure are supplied, the engine will update the MMM
% GUI during the computation, else this is considered as a run on a remote
% server that does not have access to the GUI
if exist('handles','var')
    interactive = true;
else
    interactive = false;
end

if ~isfield(options,'exhaustive')
    options.exhaustive = false;
end

% setting for the algorithm
maxat = 50000; % maximum number of heavy atoms in rigid body
forgive = 0.8; % scaling factor for van-der-Waals radii in clash tests
min_approach = 5; % minimal approach of two reference points [Å]
max_extension = 180; % maximum distance between any two reference points [Å]
clash_threshold = 1.5*forgive; % a uniform van-der-Waals radius of 1.5 Å is assumed for heavy atoms
clash_fail = 10000; %500; % maximum value of the clash cost function in testing for the unrefined model

% delete or keep SANS output fit files
delete_SANS = 'all'; % can be 'all', 'poor', or 'none', any other choice is interpreted as 'none' 
delete_SAXS = 'all'; % can be 'all', 'poor', or 'none', any other choice is interpreted as 'none' 

if isfield(restraints,'newID') && ~isempty(restraints.newID)
    PDBid = restraints.newID;
else
    PDBid = 'RIFL';
end

maxmodels = restraints.models;
pmodel = restraints.p_model;
pthr = exp(-erfinv(pmodel)^2);

trials = options.max_trials;
maxtime = options.max_time;

SANS_threshold = options.SANS_threshold; 
SAXS_threshold = options.SAXS_threshold;
xlink_threshold = options.xlink_threshold;
xlink_percentage = options.xlink_percentage;
fname = options.fname;

if options.deterministic
    rng(13);
else
    rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs
end

snum = resolve_address(sprintf('[%s]',PDBid));

snum0 = model.current_structure;
if isempty(snum)
    snum = copy_structure(snum0,PDBid);
end
add_msg_board(sprintf('Result is stored in structure %i with identifier %s\n',snum,PDBid));

% extract full atom coordinates and hevay-atom coordinates for the rigid
% bodies



% first find out whether there are stemloop libraries involved
% if so, the corresponding binding motifs must be excluded from the clash
% test
stemlib_chains = zeros(1,length(restraints.stemlibs));
for k = 1:length(restraints.stemlibs)
    adr = sprintf('[%s]%s',PDBid,restraints.stemlibs{k}.chaintag);
    indi = resolve_address(adr);
    if ~isempty(indi)
        stemlib_chains(k) = indi(2);
    end
end
%
kc = 0;
all_chains = zeros(1,100);
chain_coor = cell(length(restraints.rb),20);
heavy_coor = cell(1,length(restraints.rb));
for kr = 1:length(restraints.rb)
    coor_r = zeros(maxat,3);
    crpoi = 0;
    for kcc = 1:length(restraints.rb(kr).chains)
        kc = kc + 1;
        all_chains(kc) = restraints.rb(kr).chains(kcc);
        chain_coor{kr,kcc} = model.structures{snum}(restraints.rb(kr).chains(kcc)).xyz{1};
        if isempty(stemlib_chains) || min(abs(restraints.rb(kr).chains(kcc)-stemlib_chains)) > 0 
            [msg,coorh]=get_object([snum0 restraints.rb(kr).chains(kcc) 1],'xyz_heavy');
            if msg.error
                add_msg_board(sprintf('ERROR in reading %s coordinates',mk_address([snum0 restraints.rb(kr).chains(kcc) 1])));
                add_msg_board(msg.text);
                return
            end
            [mh,~] = size(coorh);
            coor_r(crpoi+1:crpoi+mh,:) = coorh;
            crpoi = crpoi + mh;
        end
    end
    heavy_coor{kr} = coor_r(1:crpoi,:);
end

% convert coordinates to convex hulls
for rba = 1:length(heavy_coor)
    faces = convhulln(heavy_coor{rba});
    [na,~] = size(faces);
    [kpa,ar] = reducepatch(faces,heavy_coor{rba},na);
    hulls(rba).vertices = double(ar);
    hulls(rba).faces = kpa;
end

nc = kc;
all_chains = all_chains(1:nc);

% selection is used for saving partial structures, remove old selection
if isfield(model,'selected')
    model = rmfield(model,'selected');
end
for kc = 1:nc
    model.selected{kc} = [snum all_chains(kc) 1];
end

rbnum = length(restraints.rb);

% augment lower and upper bounds
for k1 = 1:3*rbnum-1
    for k2 = k1+1:3*rbnum
        if restraints.lb(k1,k2) < min_approach
            restraints.lb(k1,k2) = min_approach;
            restraints.lb(k2,k1) = min_approach;
        end
        if restraints.ub(k1,k2) < 0.1 % the unset bounds are zero
            restraints.ub(k1,k2) = max_extension;
            restraints.ub(k2,k1) = max_extension;
        end
    end
end

if options.exhaustive
    [trials,res,trial_pattern] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
    target_trials = tph*maxtime;
    while trials > target_trials
        target_resolution = target_resolution + 0.1;
        [trials,res,trial_pattern] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
    end
    intervals = trial_pattern(:,3);
    digitbase = trial_pattern(:,4);
    irbr = trial_pattern(:,1:2);
end

[lb,ub,err]=triangle(restraints.lb,restraints.ub);

switch err
    case 0
        add_msg_board('Successful bound smoothing with experimental restraints.');
    case 1
        add_ms_board('ERROR: Some distance restraints are inconsistent.');
        success = -1;
        return
    otherwise
        ad_msg_board('Unspecified error in bound smoothing.');
        success = -1;
        return
end

met_err = 0;
embed_err = 0;
bound_err = 0;
clash_err = 0;
aux_fail = 0;
link_fail = 0;
xlink_fail = 0;
label_fail = 0;
stem_fail = 0;
stem2_fail = 0;
stem3_fail = 0;
sans_fail = 0;
saxs_fail = 0;
success = 0;
sans_poi = 0;

if isfield(restraints,'SANS')
    chi_SANS = zeros(length(restraints.SANS),50000);
end
if isfield(restraints,'SAXS')
    chi_SAXS = zeros(length(restraints.SAXS),50000);
end
if isfield(restraints,'xlinks')
    xlink_fulfill = zeros(length(restraints.xlinks),50000);
end

% the following reduces communication overhead for the parfor loop
rb = restraints.rb;
points = restraints.points;
core = restraints.core;
auxiliary = restraints.auxiliary;
links = restraints.links;
xlinks = restraints.xlinks;
sl_links = restraints.stemloop_links;
sl_DEER = restraints.SL_DEER;
stemlibs = restraints.stemlibs;
sl_atoms = 0;
for klib = 1:length(stemlibs)
    sl_coor = stemlibs{klib}.chains{1}.xyz{1};
    [msla,~] = size(sl_coor);
    sl_atoms = sl_atoms + msla;
end
[naux,~] = size(auxiliary);
refines = 0;
[ncore,~] = size(core);

sl_lib_len = zeros(1,length(restraints.stemlibs));
for lib = 1:length(restraints.stemlibs)
    sl_lib_len(lib) = length(restraints.stemlibs{lib}.chains);
end

% Determine the number of additional restraint points in stem loops for the
% rigid bodies

add_rb_points = zeros(1,length(rb));
for kl = 1:length(sl_links)
    add_rb_points(sl_links(kl).rba(1)) = add_rb_points(sl_links(kl).rba(1)) + 1;
    add_rb_points(sl_links(kl).rba(2)) = add_rb_points(sl_links(kl).rba(2)) + 1;
end
for kd = 1:length(sl_DEER)
    add_rb_points(sl_DEER(kd).rba(1)) = add_rb_points(sl_DEER(kd).rba(1)) + 1;
    add_rb_points(sl_DEER(kd).rba(2)) = add_rb_points(sl_DEER(kd).rba(2)) + 1;
end

% if cross-link threshold is larger than a specified cross-link distance,
% make sure that models are kept that conform to the threshold
for kl = 1:length(xlinks)
    if xlink_threshold > xlinks(kl).maxr
        xlinks(kl).maxr = xlink_threshold;
    end
end

final_chi2_SANS = zeros(1,maxmodels);
final_chi2_SAXS = zeros(1,maxmodels);
SANS_curves = cell(length(restraints.SANS),maxmodels);
SAXS_curves = cell(length(restraints.SAXS),maxmodels);
probabilities = zeros(1,maxmodels);
transmats = cell(1,3);

runtime = 0;
worst_res = 0;
parblocks = 0; % ### 0
bask = parblocks*options.granularity;
if options.exhaustive
    maxmodels = 200; % 200;
    maxtime = 600; % ### uncomment for production runs
end
tic,

% maxmodels = 1; % ### comment out for production runs

[pathstr,basname] = fileparts(fname);
solutionname = fullfile(pathstr,strcat(basname,'_solutions.dat'));
repname = fullfile(pathstr,strcat(basname,'_stemloops.dat'));
fid = fopen(solutionname,'wt');
fprintf(fid,'%pblock trial  combi % at granularity %i\n',options.granularity);
fclose(fid);

if ~skip_mode
    combinationname = fullfile(pathstr,strcat(basname,'_combinations.dat'));
    fid = fopen(combinationname,'wt');
    fprintf(fid,'%pblock trial  combi % at granularity %i\n',options.granularity);
    fclose(fid);
end

while runtime <= 3600*maxtime && bask < trials && success < maxmodels
    merr_vec = zeros(1,options.granularity);
    eerr_vec = zeros(1,options.granularity);
    berr_vec = zeros(1,options.granularity);
    aerr_vec = zeros(1,options.granularity);
    rerr_vec = zeros(1,options.granularity);
    lerr_vec = zeros(1,options.granularity);
    cerr_vec = zeros(1,options.granularity);
    xerr_vec = zeros(1,options.granularity);
    stem_vec = zeros(1,options.granularity);
    stem2_vec = zeros(1,options.granularity);
    stem3_vec = zeros(1,options.granularity);
    model_prob = zeros(1,options.granularity);
    scombi = zeros(options.granularity,length(stemlibs));
    all_s_combi = cell(1,options.granularity);
   
    xlink_distances = cell(1,options.granularity); 
    tmats = cell(options.granularity,1);
    % the parallel code section is separate for Monte Carlo and exhaustive
    % Rigi, since model rejection is handled differently
    parblocks = parblocks + 1;
    res_vec = res*ones(1,options.granularity);
    % profile on
    parfor kt = 1:options.granularity % ### parfor
        if solutions_given
            [msoln,nsoln] = size(solutions);
            ksoln = zeros(1,msoln);
            psoln = 0;
            skip = true;
            for kbp = 1:msoln
                if parblocks == solutions(kbp,1) && kt == solutions(kbp,2)
                    psoln = psoln + 1;
                    ksoln(psoln) = kbp;
                    skip = false;
                end
            end
            if skip
                merr_vec(kt) = 1; % skipped trials are reported as metrization errors
                continue
            end
            ksoln = ksoln(1:psoln);
            fprintf(1,'Trial %i.%i: Model will be tested.\n',parblocks,kt);
            if nsoln > 2
                for kksoln = 1:psoln
                    fprintf(1,'combination: ');
                    for kx = 3:nsoln
                        fprintf(1,'%i',solutions(kksoln,kx));
                        if kx < nsoln
                            fprintf(1,', ');
                        else
                            fprintf(1,'.\n');
                        end
                    end
                end
                fprintf(1,'\n');
            end
        end
        m0 = 0;
        fulfill = true;
        if bask + kt <= prod(intervals)
            fractions = get_restraint_fractions(bask+kt,intervals,digitbase);
            [dmatr,err,cres] = metrize_exhaustive(lb,ub,irbr,fractions,intervals);
        else
            err = 1; % report a metrization error for surplus trials
        end
        if err == 1% metrization failed (restraints inconsistent), next trial, increment error counter
            merr_vec(kt) = 1;
        else
            if cres > res_vec(kt)
                res_vec(kt) = cres;
            end
            coor0=dmat2coor(dmatr); % embed distance matrix to Cartesian space
            if isempty(coor0)
                eerr_vec(kt) = 1;
            else
                [coor1,err] = bound_refiner(coor0,lb,ub);
                if err > 0
                    berr_vec(kt) = 1;
                else
                    t_labels = cell(1,length(rb));
                    t_points = cell(1,length(rb));
                    atransmat = zeros(4*length(rb),4);
                    for kr = 1:length(rb)
                        baspoi = 3*(kr-1);
                        baspoi4 = 4*(kr-1);
                        % [~,lab_t,transmat] = rmsd_superimpose(coor1(baspoi+1:baspoi+3,:),rb(kr).ref(1:3,2:4));
                        [~,lab_t,transmat] = superimpose_3points(coor1(baspoi+1:baspoi+3,:),rb(kr).ref(1:3,2:4));
                        t_labels{kr} = lab_t;
                        t_points{kr} = affine_coor_set(points{kr},transmat);
                        atransmat(baspoi4+1:baspoi4+4,:) = transmat;
                    end
                    tmats{kt} = atransmat;
                    paux = 1;
                    for kaux = 1:naux
                        acoor1 = t_points{auxiliary(kaux,1)};
                        acoor2 = t_points{auxiliary(kaux,3)};
                        coor1b = acoor1(auxiliary(kaux,2),:);
                        coor2b = acoor2(auxiliary(kaux,4),:);
                        rsim = norm(coor1b-coor2b);
                        prob = prob_Gaussian(rsim,auxiliary(kaux,5),sqrt(auxiliary(kaux,6)^2+cres^2)); % include resolution into uncertainty
                        paux = paux*prob;
                    end
                    if paux < pthr^naux
                        aerr_vec(kt) = 1;
                        fulfill = false;
                        paux = 0;
                    end
                    % check how well label-to-label distances are reproduced
                    plabels = 1;
                    if fulfill
                        for kcore = 1:ncore
                            acoor1 = t_points{core(kcore,1)};
                            acoor2 = t_points{core(kcore,3)};
                            coor1b = acoor1(core(kcore,2),:);
                            coor2b = acoor2(core(kcore,4),:);
                            rsim = norm(coor1b-coor2b);
                            prob = prob_Gaussian(rsim,core(kcore,5),sqrt(core(kcore,6)^2+cres^2)); % include resolution into uncertainty
                            plabels = plabels*prob;
                        end
                        if plabels < pthr^ncore
                            fulfill = false;
                            rerr_vec(kt) = 1;
                        end
                    end
                    % check for linker lengths
                    if fulfill
                        if ~isempty(links(1).maxr)
                            for kl = 1:length(links)
                                if ~isempty(links(kl).ref_indices)
                                    acoor1 = t_points{links(kl).ref_indices(1)};
                                    acoor2 = t_points{links(kl).ref_indices(3)};
                                    coor1b = acoor1(links(kl).ref_indices(2),:);
                                    coor2b = acoor2(links(kl).ref_indices(4),:);
                                    rlink = norm(coor2b - coor1b);
                                    if rlink > links(kl).maxr + cres % include resolution into uncertainty
                                        fulfill = false;
                                    end
                                end
                            end
                            if ~fulfill
                                lerr_vec(kt) = 1;
                            end
                        end
                    end
                    % check for rigid-body clashes
                    if fulfill
                        clashed = false;
                        for kr1 = 1:length(rb)-1
                            baspoi4 = 4*(kr1-1);
                            atransmat = tmats{kt};
                            transmat = atransmat(baspoi4+1:baspoi4+4,:);
                            hc1 = single(affine_coor_set(heavy_coor{kr1},transmat));
                            hc1h = affine_coor_set(hulls(kr1).vertices,transmat);
                            for kr2 = kr1+1:length(rb)
                                baspoi4b = 4*(kr2-1);
                                transmatb = atransmat(baspoi4b+1:baspoi4b+4,:);
                                hc2 = single(affine_coor_set(heavy_coor{kr2},transmatb));
                                hc2h = affine_coor_set(hulls(kr2).vertices,transmatb);
                                % the fine threshold of -1 prevents the
                                % slow fine-grained clash_cost dunction to
                                % be ever called
                                cost = clash_cost_super_fast(hc1h,hulls(kr1).faces,hc2h,hulls(kr2).faces,hc1,hc2,clash_threshold,-1);
                                if cost > clash_fail
                                    clashed = true;
                                    % fprintf(2,'Clash test failed between RB%i and RB%i at cost function %12.0f\n',kr1,kr2,cost);
                                end
                            end
                        end
                        if clashed
                            cerr_vec(kt) = 1;
                            fulfill = false;
                        end
                    end
                    if fulfill
                        xld = zeros(length(xlinks),1);
                        xlf = 0;
                        for kl = 1:length(xlinks)
                            acoor1 = t_points{xlinks(kl).ref_indices(1)};
                            acoor2 = t_points{xlinks(kl).ref_indices(3)};
                            coor1b = acoor1(xlinks(kl).ref_indices(2),:);
                            coor2b = acoor2(xlinks(kl).ref_indices(4),:);
                            rlink = norm(coor2b - coor1b);
                            xld(kl) = rlink;
                            if rlink > xlinks(kl).maxr + cres % include resolution into uncertainty
                                xlf = xlf + 1;
                            end
                        end
                        xlink_distances{kt} = xld;
                        if 100*(length(xlinks)-xlf)/length(xlinks) < xlink_percentage
                            xerr_vec(kt) = 1;
                            fulfill = false;
                        end
                    end
                    if fulfill && stemloop_mode % check for stemloop linker restraints
                        libind = cell(1,length(sl_lib_len));
                        for lib = 1:length(sl_lib_len)
                            libind{lib} = 1:sl_lib_len(lib);
                        end
                        slr = length(sl_links);
                        atransmat = tmats{kt};
                        libpairs = zeros(slr,2);
                        for sl = 1:slr
                            libpairs(sl,:) = sl_links(sl).lib;
                            baspoi4 = 4*(sl_links(sl).rba(1)-1);
                            transmat1 = atransmat(baspoi4+1:baspoi4+4,:);
                            if sl_links(sl).lib(1)~= 0
                                coor1 = sl_links(sl).coor1(libind{sl_links(sl).lib(1)},:);
                            else
                                coor1 = sl_links(sl).coor1;
                            end
                            coor1b = affine_coor_set(coor1,transmat1);
                            baspoi4 = 4*(sl_links(sl).rba(2)-1);
                            transmat2 = atransmat(baspoi4+1:baspoi4+4,:);
                            if sl_links(sl).lib(2) ~= 0
                                coor2 = sl_links(sl).coor2(libind{sl_links(sl).lib(2)},:);
                            else
                                coor2 = sl_links(sl).coor2;
                            end
                            coor2b = affine_coor_set(coor2,transmat2);
                            [m1,~] = size(coor1b);
                            [m2,~] = size(coor2b);
                            a2 = repmat(sum(coor1b.^2,2),1,m2);
                            b2 = repmat(sum(coor2b.^2,2),1,m1).';
                            pair_dist = sqrt(abs(a2 + b2 - 2*coor1b*coor2b.'));
                            [ind1,ind2] = find(pair_dist <= sl_links(sl).maxr + cres);
                            %                                 if ~isempty(ind1) && ~isempty(ind2)
                            %                                     for k = 1:length(ind1), fprintf(1,'Link(%i), Combi(%i,%i): %4.2f Å\n',sl,ind1(k),ind2(k),pair_dist(ind1(k),ind2(k))); end;
                            %                                 end
                            if sl_links(sl).lib(1) ~= 0
                                oldind1 = libind{sl_links(sl).lib(1)};
                                newind1 = oldind1(ind1);
                                libind{sl_links(sl).lib(1)} = newind1;
                                for osl = 1:sl-1
                                    if libpairs(osl,1) == sl_links(sl).lib(1)
                                        if libpairs(osl,2) ~= 0
                                            oldind3 = libind{libpairs(osl,2)};
                                            newind3 = oldind3(ind1);
                                            libind{libpairs(osl,2)} = newind3;
                                        end
                                    end
                                    if libpairs(osl,2) == sl_links(sl).lib(1)
                                        if libpairs(osl,1) ~= 0
                                            oldind3 = libind{libpairs(osl,1)};
                                            newind3 = oldind3(ind1);
                                            libind{libpairs(osl,1)} = newind3;
                                        end
                                    end
                                end
                            else
                                newind1 = ind1;
                            end
                            if sl_links(sl).lib(2) ~= 0
                                oldind2 = libind{sl_links(sl).lib(2)};
                                newind2 = oldind2(ind2);
                                libind{sl_links(sl).lib(2)} = newind2;
                                for osl = 1:sl-1
                                    if libpairs(osl,1) == sl_links(sl).lib(2)
                                        if libpairs(osl,2) ~= 0
                                            oldind3 = libind{libpairs(osl,2)};
                                            newind3 = oldind3(ind2);
                                            libind{libpairs(osl,2)} = newind3;
                                        end
                                    end
                                    if libpairs(osl,2) == sl_links(sl).lib(2)
                                        if libpairs(osl,1) ~= 0
                                            oldind3 = libind{libpairs(osl,1)};
                                            newind3 = oldind3(ind2);
                                            libind{libpairs(osl,1)} = newind3;
                                        end
                                    end
                                end
                            else
                                newind2 = ind2;
                            end
                            if isempty(newind1) || isempty(newind2)
                                stem_vec(kt) = 1;
                                fulfill = false;
                                break
                            end
                        end
                        if fulfill && ~isempty(libind)
                            ind1 = unique(libind{1});
                            ind1 = 1; % ### this is a hack ###
                            libind{1} = ind1;
                            m0 = length(ind1);
                            combinations = zeros(m0,length(stemlibs));
                            combinations(:,1) = ind1';
                            for sl = 2:length(stemlibs)
                                ind = unique(libind{sl});
                                libind{sl} = ind;
                                m = length(ind);
                                combinations = repmat(combinations,m,1);
                                for sc = 1:m
                                    bas = (sc-1)*m0;
                                    combinations(bas+1:bas+m0,sl) = ind(sc)*ones(m0,1);
                                end
                                m0 = m0*m;
                            end
                        end
                    end
                    if fulfill && stemloop_mode && ~isempty(sl_DEER) % check for stemloop label restraints
                        libpairs = zeros(length(sl_DEER),2);
                        for sld = 1:length(sl_DEER)
                            libpairs(sld,:) = sl_DEER(sld).stemlib;
                            baspoi4 = 4*(sl_DEER(sld).rba(1)-1);
                            transmat1 = atransmat(baspoi4+1:baspoi4+4,:);
                            baspoi4 = 4*(sl_DEER(sld).rba(2)-1);
                            transmat2 = atransmat(baspoi4+1:baspoi4+4,:);
                            if sl_DEER(sld).stemlib(1) == 0
                                ind1 = 1;
                                coor1 = sl_DEER(sld).coor1;
                            else
                                ind1 = libind{sl_DEER(sld).stemlib(1)};
                                coor1 = sl_DEER(sld).coor1(ind1,:);
                            end
                            if sl_DEER(sld).stemlib(2) == 0
                                ind2 = 1;
                                coor2 = sl_DEER(sld).coor2;
                            else
                                ind2 = libind{sl_DEER(sld).stemlib(2)};
                                coor2 = sl_DEER(sld).coor2(ind2,:);
                            end
                            coor1b = affine_coor_set(coor1,transmat1);
                            coor2b = affine_coor_set(coor2,transmat2);
                            [m1,~] = size(coor1b);
                            [m2,~] = size(coor2b);
                            a2 = repmat(sum(coor1b.^2,2),1,m2);
                            b2 = repmat(sum(coor2b.^2,2),1,m1).';
                            pair_dist = sqrt(abs(a2 + b2 - 2*coor1b*coor2b.'));
                            minr = sl_DEER(sld).r - sqrt(sl_DEER(sld).sigr^2+cres^2);
                            maxr = sl_DEER(sld).r + sqrt(sl_DEER(sld).sigr^2+cres^2);
                            [nind1,nind2] = find(pair_dist >= minr & pair_dist <= maxr);
                            %                                 nind1 = unique(nind1);
                            %                                 nind2 = unique(nind2);
                            if sl_DEER(sld).stemlib(1) ~= 0
                                oldind1 = libind{sl_DEER(sld).stemlib(1)};
                                newind1 = oldind1(nind1);
                                libind{sl_DEER(sld).stemlib(1)} = newind1;
                                for osl = 1:sld-1
                                    if libpairs(osl,1) == sl_DEER(sld).stemlib(1)
                                        if libpairs(osl,2) ~= 0
                                            oldind3 = libind{libpairs(osl,2)};
                                            newind3 = oldind3(nind1);
                                            libind{libpairs(osl,2)} = newind3;
                                        end
                                    end
                                    if libpairs(osl,2) == sl_DEER(sld).stemlib(1)
                                        if libpairs(osl,1) ~= 0
                                            oldind3 = libind{libpairs(osl,1)};
                                            newind3 = oldind3(nind1);
                                            libind{libpairs(osl,1)} = newind3;
                                        end
                                    end
                                end
                            else
                                newind1 = nind1;
                            end
                            if sl_DEER(sld).stemlib(2) ~= 0
                                oldind2 = libind{sl_DEER(sld).stemlib(2)};
                                newind2 = oldind2(nind2);
                                libind{sl_DEER(sld).stemlib(2)} = newind2;
                                for osl = 1:sld-1
                                    if libpairs(osl,1) == sl_DEER(sld).stemlib(2)
                                        if libpairs(osl,2) ~= 0
                                            oldind3 = libind{libpairs(osl,2)};
                                            newind3 = oldind3(nind2);
                                            libind{libpairs(osl,2)} = newind3;
                                        end
                                    end
                                    if libpairs(osl,2) == sl_DEER(sld).stemlib(2)
                                        if libpairs(osl,1) ~= 0
                                            oldind3 = libind{libpairs(osl,1)};
                                            newind3 = oldind3(nind2);
                                            libind{libpairs(osl,1)} = newind3;
                                        end
                                    end
                                end
                            else
                                newind2 = nind2;
                            end
                            if isempty(newind1) || isempty(newind2)
                                fulfill = false;
                                stem2_vec(kt) = 1;
                                break
                            end
                        end
                        if fulfill
                            ind1 = unique(libind{1});
                            libind{1} = ind1;
                            m0 = length(ind1);
                            combinations = zeros(m0,length(stemlibs));
                            combinations(:,1) = ind1';
                            for sl = 2:length(stemlibs)
                                ind = unique(libind{sl});
                                libind{sl} = ind;
                                m = length(ind);
                                combinations = repmat(combinations,m,1);
                                for sc = 1:m
                                    bas = (sc-1)*m0;
                                    combinations(bas+1:bas+m0,sl) = ind(sc)*ones(m0,1);
                                end
                                m0 = m0*m;
                            end
                        end
                    end
                    if fulfill
                        atransmat = tmats{kt};
                        if ~isempty(stemlibs) && stemloop_mode
                            % add stemloop anchors and label coordinates to
                            % points and restraints to auxiliary and links,
                            % if necessary
                            points1 = cell(1,length(rb));
                            kr_poi0 = zeros(1,length(rb));
                            for kr = 1:length(rb)
                                [mrp,~] = size(points{kr});
                                points1{kr} = zeros(mrp + add_rb_points(kr),3);
                                kr_poi0(kr) = mrp;
                                points1{kr}(1:mrp,:) = points{kr};
                            end
                            naux1 = naux + length(sl_DEER);
                            auxiliary1 = zeros(naux1,6);
                            auxiliary1(1:naux,:) = auxiliary;
                            links1 = links;
                            min_cost = 1e6;
                            csuccess = 0;
                            atransmat0 = atransmat;
                            batransmat = atransmat;
                            berrvec = zeros(1,5);
                            berrvec(3) = 1;
                            bmprob = 0;
                            bxld = [];
                            ref_indices = zeros(1,4);
                            if solutions_given
                                fprintf(1,'Trial %i.%i: Up to %i combination(s) will be tested\n',parblocks,kt,m0);
                            end
                            kcvec = 1:m0;
                            scramble = rand(1,m0);
                            [~,scrambler] = sort(scramble);
                            kcvec = kcvec(scrambler);
                            refineable = zeros(1,m0);
                            for kcombi = kcvec
                                if solutions_given && nsoln > 2
                                    valid = false;
                                    for kksoln = 1: length(ksoln)
                                        if sum(abs(combinations(kcombi,1:nsoln-2) - solutions(ksoln(kksoln),3:nsoln))) == 0
                                            valid = true;
                                        end
                                    end
                                    if ~valid
                                        continue
                                    end
                                end
                                kr_poi = kr_poi0;
                                auxpoi = naux;
                                linkpoi = length(links);
                                for kl = 1:length(sl_DEER)
                                    rba1 = sl_DEER(kl).rba(1);
                                    rba2 = sl_DEER(kl).rba(2);
                                    sl1 = sl_DEER(kl).stemlib(1);
                                    if sl1 ~= 0
                                        sl1 = combinations(kcombi,sl_DEER(kl).stemlib(1));
                                        coor1 = sl_DEER(kl).coor1(sl1,:);
                                    else
                                        coor1 = sl_DEER(kl).coor1;
                                    end
                                    sl2 = sl_DEER(kl).stemlib(2);
                                    if sl2 ~= 0
                                        sl2 = combinations(kcombi,sl_DEER(kl).stemlib(2));
                                        coor2 = sl_DEER(kl).coor2(sl2,:);
                                    else
                                        coor2 = sl_DEER(kl).coor2;
                                    end
                                    auxpoi = auxpoi + 1;
                                    auxiliary1(auxpoi,1) = rba1;
                                    kr_poi(rba1) = kr_poi(rba1) + 1;
                                    points1{rba1}(kr_poi(rba1),:) = coor1;
                                    auxiliary1(auxpoi,2) = kr_poi(rba1);
                                    auxiliary1(auxpoi,3) = rba2;
                                    kr_poi(rba2) = kr_poi(rba2) + 1;
                                    points1{rba2}(kr_poi(rba2),:) = coor2;
                                    auxiliary1(auxpoi,4) = kr_poi(rba2);
                                    auxiliary1(auxpoi,5) = sl_DEER(kl).r;
                                    auxiliary1(auxpoi,6) = sl_DEER(kl).sigr;
                                end
                                for kl = 1:length(sl_links)
                                    rba1 = sl_links(kl).rba(1);
                                    rba2 = sl_links(kl).rba(2);
                                    sl1 = sl_links(kl).lib(1);
                                    if sl1 ~= 0
                                        sl1 = combinations(kcombi,sl_links(kl).lib(1));
                                        coor1 = sl_links(kl).coor1(sl1,:);
                                    else
                                        coor1 = sl_links(kl).coor1;
                                    end
                                    sl2 = sl_links(kl).lib(2);
                                    if sl2 ~= 0
                                        sl2 = combinations(kcombi,sl_links(kl).lib(2));
                                        coor2 = sl_links(kl).coor2(sl2,:);
                                    else
                                        coor2 = sl_links(kl).coor2;
                                    end
                                    linkpoi = linkpoi + 1;
                                    ref_indices(1) = rba1;
                                    ref_indices(3) = rba2;
                                    kr_poi(rba1) = kr_poi(rba1) + 1;
                                    points1{rba1}(kr_poi(rba1),:) = coor1;
                                    ref_indices(2) = kr_poi(rba1);
                                    kr_poi(rba2) = kr_poi(rba2) + 1;
                                    points1{rba2}(kr_poi(rba2),:) = coor2;
                                    ref_indices(4) = kr_poi(rba2);
                                    links1(linkpoi).ref_indices = ref_indices;
                                    links1(linkpoi).maxr = sl_links(kl).maxr;
                                end
                                % the first refinement is performed without
                                % testing for clashes, as this is much faster
                                % fprintf(1,'Refining #%i (%i,%i,%i)\n',kcombi,combinations(kcombi,:));
                                [atransmat,errvec] = ...
                                    refine_rba(rb,atransmat0,points1,pthr,naux1,auxiliary1,ncore,core,links1,1e6,heavy_coor,xlink_percentage,xlinks,stemlibs,combinations(kcombi,:));
                                if ~sum(errvec)
                                    fprintf(1,'R%i.%i: Secondary refining (%i,%i,%i)\n',parblocks,kt,combinations(kcombi,:));
                                    [atransmat,errvec,mprob,xld,cost] = ...
                                        refine_rba_fast(rb,atransmat,points1,pthr,naux1,auxiliary1,ncore,core,links1,clash_threshold,heavy_coor,xlink_percentage,xlinks,stemlibs,combinations(kcombi,:));
                                    if ~sum(errvec)
                                        fprintf(1,'R%i.%i: Successfully fitted (%i,%i,%i)\n',parblocks,kt,combinations(kcombi,:));
                                        refineable(kcombi) = 1;
                                        if cost < min_cost
                                            bcombi = kcombi;
                                        end
                                        batransmat = atransmat;
                                        berrvec = errvec;
                                        bmprob = mprob;
                                        bxld = xld;
                                        csuccess = csuccess + 1;
                                        if skip_mode
                                            break
                                        end
                                    else
                                        % fprintf(1,'(%i,%i,%i) Initial SL clash score: %8.1f Cost: %8.1f (R%i.%i)\n',combinations(kcombi,:),sl_clash_scores(kcombi),cost,parblocks,kt);
                                    end
                                end
                            end
                            atransmat = batransmat;
                            errvec = berrvec;
                            mprob = bmprob;
                            xld = bxld;
                            if csuccess > 0
                                scombi(kt,:) = combinations(bcombi,:);
                                if ~skip_mode
                                    all_s_combi{kt} = combinations(refineable~=0,:);
                                end
                                % fprintf(1,'R%i.%i: combination (%i, %i, %i) was successfully refined\n',parblocks,kt,scombi(kt,:));
                            end
                            if csuccess > 1 && ~isempty(combinations) && skip_mode
                                fprintf(2,'R%i.%i: %i combinations were successfully refined, but only (%i, %i, %i) is kept\n',parblocks,kt,csuccess,scombi(kt,:));
                            end
                            tmats{kt} = atransmat;
                            xlink_distances{kt} = xld;
                            model_prob(kt) = mprob;
                            aerr_vec(kt) = errvec(1);
                            rerr_vec(kt) = errvec(2);
                            lerr_vec(kt) = errvec(3);
                            cerr_vec(kt) = errvec(4);
                            xerr_vec(kt) = errvec(5);
                        else % treatment in the absence of stemlibs or if stemloop_mode is false
                            % the first refinement is performed without
                            % testing for clashes, as this is much faster
                            if solutions_given && trial_echo
                                fprintf(1,'Trial %i.%i: Model will be refined\n',parblocks,kt);
                            end
                            atransmat0 = atransmat;
                            [atransmat,errvec,mprob,xld] = ...
                                refine_rba(rb,atransmat0,points,pthr,naux,auxiliary,ncore,core,links,1e6,heavy_coor,xlink_percentage,xlinks);
                            if ~sum(errvec)
                                if solutions_given && trial_echo
                                    fprintf(1,'Trial %i.%i: Model will be refined with clash score\n',parblocks,kt);
                                end
                                [atransmat,errvec,mprob,xld,cost,costs] = ...
                                    refine_rba_fast(rb,atransmat,points,pthr,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,xlink_percentage,xlinks);
                                if trial_echo
                                    fprintf(1,'Total cost: %5.2f. aux. %5.2f; core %5.2f; links %5.2f; clash %5.2f; stem %5.2f\n',cost,costs.aux,costs.core,costs.link,costs.clash,costs.stem);
                                end
                            end
                            tmats{kt} = atransmat;
                            xlink_distances{kt} = xld;
                            model_prob(kt) = mprob;
                            aerr_vec(kt) = errvec(1);
                            rerr_vec(kt) = errvec(2);
                            lerr_vec(kt) = errvec(3);
                            cerr_vec(kt) = errvec(4);
                            xerr_vec(kt) = errvec(5);
                        end
                    end
                end
            end
        end
    end
    % profile viewer
    if ~solutions_given
        fprintf(1,'Parallel block %i completed.\n',parblocks);
    end
    % analyze all models generated in the parallel code section
    met_err = met_err + sum(merr_vec);
    embed_err = embed_err + sum(eerr_vec);
    bound_err = bound_err + sum(berr_vec);
    aux_fail = aux_fail + sum(aerr_vec);
    label_fail = label_fail + sum(rerr_vec);
    link_fail = link_fail + sum(lerr_vec);
    clash_err = clash_err + sum(cerr_vec);
    xlink_fail = xlink_fail + sum(xerr_vec);
    stem_fail = stem_fail + sum(stem_vec);
    stem2_fail = stem2_fail + sum(stem2_vec);
    stem3_fail = stem3_fail + sum(stem3_vec);
    fail_vec = merr_vec+eerr_vec+berr_vec+aerr_vec+rerr_vec+lerr_vec+cerr_vec+xerr_vec;
    fail_vec = fail_vec + stem_vec+stem2_vec+stem3_vec;
    if worst_res < max(res_vec)
        worst_res = max(res_vec);
    end
    for k = bask+1:bask+options.granularity
        if fail_vec(k-bask) == 0
            if isempty(all_s_combi{k-bask})
                all_combi = scombi(k-bask,:);
            else
                all_combi = all_s_combi{k-bask};
            end
            atransmat = tmats{k-bask};
            for kr = 1:length(rb)
                baspoi4 = 4*(kr-1);
                transmats{kr} = atransmat(baspoi4+1:baspoi4+4,:);
            end
            [ncombi,~] = size(all_combi);
            for kcombi = 1:ncombi
                success = success + 1;
                ccombi = all_combi(kcombi,:);
                probabilities(success) = model_prob(k-bask)^(1/(naux+ncore));
                tmstd = [];
                if stemloop_mode
                    % replace binding motifs by stemloops from library if stemmloop
                    % libraries are tested
                    for klib = 1:length(stemlibs)
                        for kr = 1:length(restraints.rb)
                            for kcc = 1:length(restraints.rb(kr).chains)
                                if stemlibs{klib}.cind(2) == restraints.rb(kr).chains(kcc)
                                    chain_coor{kr,kcc} = stemlibs{klib}.chains{ccombi(klib)}.xyz{1};
                                    if success == 1
                                        fields = fieldnames(model.structures{snum}(restraints.rb(kr).chains(kcc)));
                                        for kfield = 1:length(fields)
                                            if isfield(stemlibs{klib}.chains{ccombi(klib)},fields{kfield})
                                                model.structures{snum}(restraints.rb(kr).chains(kcc)).(fields{kfield}) = stemlibs{klib}.chains{ccombi(klib)}.(fields{kfield});
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                t_chain_coor = chain_coor;
                for kr = 1:length(restraints.rb)
                    for kc = 1:length(restraints.rb(kr).chains)
                        t_chain_coor{kr,kc} = affine_coor_set(chain_coor{kr,kc},transmats{kr});
                        if isfield(restraints,'superimpose') && restraints.rb(kr).chains(kc) == restraints.superimpose
                            [~,~,tmstd] = rmsd_superimpose(chain_coor{kr,kc},t_chain_coor{kr,kc});
                        end
                    end
                end
                for kr = 1:length(restraints.rb)
                    for kc = 1:length(restraints.rb(kr).chains)
                        if ~isempty(tmstd)
                            t_chain_coor{kr,kc} = affine_coor_set(t_chain_coor{kr,kc},tmstd);
                        end
                        model.structures{snum}(restraints.rb(kr).chains(kc)).xyz{success} = t_chain_coor{kr,kc};
                        if success > 1
                            model.structures{snum}(restraints.rb(kr).chains(kc)).atoms{success} = ...
                                model.structures{snum}(restraints.rb(kr).chains(kc)).atoms{1};
                            model.structures{snum}(restraints.rb(kr).chains(kc)).residues{success} = ...
                                model.structures{snum}(restraints.rb(kr).chains(kc)).residues{1};
                            model.structures{snum}(restraints.rb(kr).chains(kc)).Bfactor{success} = ...
                                model.structures{snum}(restraints.rb(kr).chains(kc)).Bfactor{1};
                            model.structures{snum}(restraints.rb(kr).chains(kc)).Btensor{success} = ...
                                model.structures{snum}(restraints.rb(kr).chains(kc)).Btensor{1};
                        end
                    end
                end
                % now it should be compared to SANS or SAXS restraints,
                % if any
                if solutions_given
                    if isempty(ccombi) || sum(ccombi) == 0
                        fprintf(1,'Trial %i.%i is compared with SAS fits\n',parblocks,k-bask);
                    else
                        fprintf(1,'Trial %i.%i: Combination (%i,%i,%i) is compared with SAS fits\n',parblocks,k-bask,ccombi(:));
                    end
                end
                SANS_chi = 0;
                fulfill = true;
                % now compare stemloop restraints
                if ~stemloop_mode && ~isempty(restraints.stemlibs) && fulfill
                    modnum.block = parblocks;
                    modnum.num = k-bask;
                    modnum.model = success;
                    sl_solutions = fit_stemloop_combinations(restraints,snum0,snum,repname,modnum);
                    if isempty(sl_solutions)
                        success = success - 1;
                        stem_fail = stem_fail + 1;
                        fulfill = false;
                    end
                end
                if fulfill
%                     if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
%                         to_be_deleted = '';
%                         sans_vec = -ones(2,1);
%                         for ks = 1:length(restraints.SANS)
%                             model = rmfield(model,'selected');
%                             ksel = 0;
%                             for kc = 1:length(restraints.SANS(ks).chains)
%                                 if restraints.SANS(ks).chains(kc) > 0
%                                     ksel = ksel + 1;
%                                     model.selected{ksel} = [snum restraints.SANS(ks).chains(kc) success];
%                                 end
%                             end
%                             pdbfile = sprintf('t%i_%i',k,ks);
%                             to_be_deleted = sprintf('t%i_*.*',k);
%                             wr_pdb_selected(pdbfile,'SANS');
%                             [chi2,~,~,result,fit] = fit_SANS_by_cryson(restraints.SANS(ks).data,pdbfile,restraints.SANS(ks).illres);
%                             if isempty(chi2) || isnan(chi2)
%                                 SANS_chi = 1e6;
%                                 if interactive
%                                     fprintf(2,'Warning: SANS fitting failed in trial %i:\n',k);
%                                     fprintf(2,'%s',result);
%                                     if success > 0
%                                         success = success - 1;
%                                     end
%                                     fulfill = false;
%                                 end
%                             else
%                                 SANS_curves{ks,success} = fit;
%                                 sans_vec(ks) = chi2;
%                                 SANS_chi = SANS_chi + chi2;
%                             end
%                         end
%                         if min(sans_vec) > 0
%                             sans_poi = sans_poi + 1;
%                             chi_SANS(:,sans_poi) = sans_vec;
%                             xlink_fulfill(:,sans_poi) = xlink_distances{k-bask};
%                         end
%                         chi2 = SANS_chi/length(restraints.SANS);
%                         if chi2 > SANS_threshold && fulfill
%                             if interactive
%                                 fprintf(1,'SANS chi^2 of %4.2f exceeded threshold of %4.2f in trial %i.\n',chi2,SANS_threshold,k);
%                             end
%                             success = success - 1;
%                             sans_fail = sans_fail + 1;
%                             fulfill = false;
%                             if strcmp(delete_SANS,'all') || strcmp(delete_SANS,'poor')
%                                 delete(strcat(to_be_deleted,'*.*'));
%                             end
%                         else
%                             if strcmp(delete_SANS,'all') && ~isempty(to_be_deleted)
%                                 delete(strcat(to_be_deleted,'*.*'));
%                             end
%                             final_chi2_SANS(success) = chi2;
%                             if interactive && options.display_SANS_fit
%                                 % update multi plot axes
%                                 axes(handles.axes_multi_plot);
%                                 cla;
%                                 hold on;
%                                 for ks = 1:length(restraints.SANS)
%                                     fit = SANS_curves{ks,success};
%                                     plot(fit(:,1),fit(:,2));
%                                     plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
%                                 end
%                                 title(sprintf('SANS fit for model %i (chi^2 = %4.2f, p = %4.2f)',success,final_chi2_SANS(success),probabilities(success)));
%                                 drawnow
%                             end
%                         end
%                     else
%                         sans_poi = sans_poi + 1;
%                         xlink_fulfill(:,sans_poi) = xlink_distances{k-bask};
%                     end
%                     SAXS_chi = 0;
%                     if isfield(restraints,'SAXS') && fulfill
%                         to_be_deleted = '';
%                         for ks = 1:length(restraints.SAXS)
%                             model = rmfield(model,'selected');
%                             ksel = 0;
%                             for kc = 1:length(restraints.SAXS(ks).chains)
%                                 if restraints.SAXS(ks).chains(kc) > 0
%                                     ksel = ksel + 1;
%                                     model.selected{ksel} = [snum restraints.SAXS(ks).chains(kc) success];
%                                 end
%                             end
%                             pdbfile = sprintf('tx%i_%i',k,ks);
%                             to_be_deleted = sprintf('tx%i_*.*',k);
%                             wr_pdb_selected(pdbfile,'SAXS');
%                             SAXS_options.sm = 10*restraints.SAXS(ks).sm;
%                             [chi2,~,~,result,fit] = fit_SAXS_by_crysol(restraints.SAXS(ks).data,pdbfile,SAXS_options);
%                             if isnan(chi2)
%                                 chi2 = 1e6;
%                             end
%                             if isempty(chi2)
%                                 if interactive
%                                     fprintf(2,'Warning: SAXS fitting failed in trial %i:\n',k);
%                                     fprintf(2,'%s',result);
%                                 end
%                                 delete(strcat(pdbfile,'*.*'));
%                                 success = success - 1;
%                                 fulfill = false;
%                                 saxs_fail = saxs_fail + 1;
%                             else
%                                 chi_SAXS(ks,success) = chi2;
%                                 SAXS_curves{ks,success} = fit;
%                                 SAXS_chi = SAXS_chi + chi2;
%                             end
%                         end
%                         chi2 = SAXS_chi/length(restraints.SAXS);
%                         if chi2 > SAXS_threshold
%                             success = success - 1;
%                             saxs_fail = saxs_fail + 1;
%                             fulfill = false;
%                             if strcmp(delete_SAXS,'all') || strcmp(delete_SAXS,'poor')
%                                 delete(strcat(to_be_deleted,'*.*'));
%                             end
%                         elseif ~isempty(chi2)
%                             final_chi2_SAXS(success) = chi2;
%                             if interactive && options.display_SAXS_fit
%                                 % update multi plot axes
%                                 axes(handles.axes_multi_plot);
%                                 cla;
%                                 hold on;
%                                 for ks = 1:length(restraints.SAXS)
%                                     fit = SAXS_curves{ks,success};
%                                     plot(fit(:,1),fit(:,2));
%                                     plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
%                                 end
%                                 title(sprintf('SAXS fit for model %i (chi^2 = %4.2f, p = %4.2f)',success,final_chi2_SAXS(success),probabilities(success)));
%                                 drawnow
%                             end
%                             if strcmp(delete_SAXS,'all') && ~isempty(to_be_deleted)
%                                 delete(strcat(to_be_deleted,'*.*'));
%                             end
%                         end
%                     end
                    if fulfill
                        fid = fopen(solutionname,'at');
                        if ~isempty(ccombi) && sum(ccombi) > 0
                            fprintf(fid,'%8i%6i%6i%6i%6i\n',parblocks,k-bask,ccombi(:));
                        else
                            fprintf(fid,'%8i%6i\n',parblocks,k-bask);
                        end
                        fclose(fid);
                        if restraints.search
                            success = success -1;
                        end
                    end
                end
            end
            if ~skip_mode && fulfill
                fid = fopen(combinationname,'at');
                all_combi = all_s_combi{k-bask};
                [ncombi,~] = size(all_combi);
                for kcombi = 1:ncombi
                    fprintf(fid,'%8i%6i%6i%6i%6i\n',parblocks,k-bask,all_combi(kcombi,:)); 
                end
                fclose(fid);
            end
        end
    end
    bask = bask + options.granularity;
    % save SANS_fitting_RNA chi_SANS sans_poi
    runtime=toc;
    if interactive
        % update multi plot axes
%         if options.display_SANS_chi2
%             axes(handles.axes_multi_plot);
%             cla;
%             hold on;
%             for kl = 1:length(restraints.SANS)
%                 plot(chi_SANS(kl,1:sans_poi),'.');
%             end
%             plot(sum(chi_SANS(:,1:sans_poi))/sans_poi,'k');
%             plot([1,sans_poi],[SANS_threshold,SANS_threshold],':','Color',[0.75,0,0]);
%             title('SANS fulfillment');
%         end
        if options.display_xlinks
            axes(handles.axes_multi_plot);
            cla;
            hold on;
            for kl = 1:length(xlinks)
                plot(xlink_fulfill(kl,1:sans_poi),'.');
            end
            plot([1,sans_poi],[xlink_threshold,xlink_threshold],':','Color',[0.75,0,0]);
            title('Crosslink fulfillment');
        end
        %  update runtime information
        left_time1 = 3600*maxtime - runtime; 
        left_time2 = runtime*(trials-bask)/bask;
        left_time3 = runtime*(maxmodels-success)/success;
        if isnan(left_time3)
            left_time3 = maxtime;
        end
        left_time = min([left_time1 left_time2 left_time3]);
        hours = floor(left_time/3600);
        minutes = round((left_time-3600*floor(left_time/3600))/60);
        handles.text_time_left.String = sprintf('%i h %i min estimated run time to completion.',hours,minutes);
        handles.text_time_left.ForegroundColor = [162,20,47]/256;
        handles.text_max_trials.String = sprintf('%6.2f',100*bask/trials);
        handles.text_percent_time.String = sprintf('%6.2f',100*runtime/(3600*maxtime));
        handles.text_success.String = sprintf('%i',success);
        left_trials = bask;
        dmg_err = met_err + embed_err + bound_err;
        handles.text_dmg_fail.String = sprintf('%6.2f',100*dmg_err/left_trials);
        left_trials = left_trials - dmg_err;
        handles.text_auxiliary_fail.String = sprintf('%6.2f',100*aux_fail/left_trials);
        left_trials = left_trials - aux_fail;
        handles.text_core_fail.String = sprintf('%6.2f',100*label_fail/left_trials);
        left_trials = left_trials - label_fail;
        handles.text_linker_fail.String = sprintf('%6.2f',100*link_fail/left_trials);
        left_trials = left_trials - link_fail;
        handles.text_clash_fail.String = sprintf('%6.2f',100*clash_err/left_trials);
        left_trials = left_trials - clash_err;
        handles.text_xlink_fail.String = sprintf('%6.2f',100*xlink_fail/left_trials);
        left_trials = left_trials - xlink_fail;
        handles.text_SANS_fail.String = sprintf('%6.2f',100*sans_fail/left_trials);
        left_trials = left_trials - sans_fail;
        handles.text_SAXS_fail.String = sprintf('%6.2f',100*saxs_fail/left_trials);
        left_trials = left_trials - saxs_fail;
        % fprintf(1,'Stem link failures: %6.2f%%\n',100*stem_fail/left_trials);
        left_trials = left_trials - stem_fail;
        if left_trials ~= success
            fprintf(2,'Trial dissipation. Expected success: %i. Found success %i.\n',left_trials,success);
        end
        if success > 0
            handles.text_time_per_model.String = sprintf('%12.1f',runtime/success);
        end
        drawnow
    end
end
toc,

if success > maxmodels
    success = maxmodels;
end

diagnostics.runtime = runtime;
diagnostics.trials = bask;
diagnostics.met_err = met_err;
diagnostics.embed_err = embed_err;
diagnostics.bound_err = bound_err;
diagnostics.aux_fail = aux_fail;
diagnostics.label_fail = label_fail;
diagnostics.link_fail = link_fail;
diagnostics.clash_err = clash_err;
diagnostics.xlink_fail = xlink_fail;
% diagnostics.sans_fail = sans_fail;
% diagnostics.saxs_fail = saxs_fail;
diagnostics.success = success;
diagnostics.probabilities = probabilities(1:success);
diagnostics.snum = snum;
diagnostics.resolution = worst_res;
diagnostics.exhaustive = options.exhaustive;
diagnostics.exhaustive_completed = bask >= trials;


% if isfield(restraints, 'SANS')
%     [~,n] = size(SANS_curves);
%     if n >= success
%         diagnostics.final_chi2_SANS = final_chi2_SANS(1:success);
%         diagnostics.SANS_curves = SANS_curves(:,1:success);
%         diagnostics.chi_SANS = chi_SANS(:,1:sans_poi);
%         diagnostics.xlink_fulfill = xlink_fulfill(:,1:sans_poi);
%     else
%         diagnostics.SANS_curves = [];
%         diagnostics.final_chi2_SANS = [];
%         diagnostics.chi_SANS = [];
%         diagnostics.xlink_fulfill = [];
%     end
% end
% 
% if isfield(restraints, 'SAXS')
%     [~,n] = size(SAXS_curves);
%     if n >= success
%         diagnostics.SAXS_curves = SAXS_curves(:,1:success);
%         diagnostics.final_chi2_SAXS = final_chi2_SAXS(1:success);
%     end
% end

return % ###

if success > 0
    model = rmfield(model,'selected');
    spoi = 0;
    for kc = 1:7
        for km = 1:success
            spoi = spoi + 1;
            model.selected{spoi} = [snum kc km];
        end
    end
    message = wr_pdb_selected(fname,PDBid);
    if message.error 
        diagnostics.unsaved = true;
        if interactive
            add_msg_board(sprintf(2,'Warning: Model could not be automatically saved. %s\n',message.text));
        end
    else
        diagnostics.unsaved = false;
    end
    scriptname = fullfile(pathstr,strcat(basname,'.mmm'));
    fid = fopen(scriptname,'wt');
    for kr = 1:length(restraints.rb)
        for kc = 1:length(restraints.rb(kr).chains)
            [~,ctag] = mk_address_parts([snum,restraints.rb(kr).chains(kc)]);
            fprintf(fid,'show [%s](%s){:} ribbon\n',PDBid,ctag);
            if isfield(restraints,'color')
                fprintf(fid,'color [%s](%s){:} %s\n',PDBid,ctag,restraints.color{restraints.rb(kr).chains(kc)});
            end
        end
    end
%     prob = 1./(final_chi2_SANS.^2+final_chi2_SAXS.^2);
%     prob_norm = 2*max(prob); 
    for km = 1:success        
        % alpha = prob(km)/prob_norm;
%         if isnan(alpha)
%             alpha = 0.5;
%         end;
        alpha = probabilities(km);
        fprintf(fid,'transparency [%s](:){%i} %5.3f\n',PDBid,km,alpha);
    end
    fclose(fid);
end


