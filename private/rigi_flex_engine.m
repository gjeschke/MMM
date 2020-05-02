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
% G. Jeschke, 6.9.2019

global model

if ~isfield(restraints,'core')
    restraints.core = [];
end
if ~isfield(restraints,'auxiliary')
    restraints.auxiliary = [];
end

solutions_given = false;
solutions = [];
processed = [];
soln_count = 0;
reference_geometry = cell(1,50000);

if isfield(restraints,'solutions') && ~isempty(restraints.solutions)
    solutions_given = true;
    poi = strfind(restraints.solutions,'.dat');
    if isempty(poi)
        soln_name = strcat(restraints.solutions,'_solutions.dat');
        proc_name = strcat(restraints.solutions,'_processed.dat');
    else
        soln_name = restraints.solutions;
    end
    solutions = load(soln_name);
    solutions = round(solutions);
    if exist(proc_name,'file')
        processed = load(proc_name);
        processed = round(processed);
    else
        processed = [];
    end
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

target_resolution = 3; % resolution with spin labels is not realistically better than 3 ?
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
min_approach = 5; % minimal approach of two reference points [?]
max_extension = 180; % maximum distance between any two reference points [?]
clash_threshold = 1.5*forgive; % a uniform van-der-Waals radius of 1.5 ? is assumed for heavy atoms
clash_fail = 10000; %500; % maximum value of the clash cost function in testing for the unrefined model


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

xlink_threshold = options.xlink_threshold;
xlink_percentage = options.xlink_percentage;
fname = options.fname;

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

hulls(length(heavy_coor)).vertices = 0;
hulls(length(heavy_coor)).faces = 0;
% convert coordinates to convex hulls
for rba = 1:length(heavy_coor)
    faces = convhulln(heavy_coor{rba});
    [na,~] = size(faces);
    [kpa,ar] = reducepatch(faces,heavy_coor{rba},na);
    hulls(rba).vertices = double(ar);
    hulls(rba).faces = kpa;
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
    if restraints.maxtrials > 0
        while trials > restraints.maxtrials
            target_resolution = target_resolution + 0.1;
            [trials,res,trial_pattern] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
        end
    end
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
        diagnostics.success = -1;
        return
    otherwise
        ad_msg_board('Unspecified error in bound smoothing.');
        diagnostics.success = -1;
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

transmats = cell(1,3);

runtime = 0;
worst_res = 0;
parblocks = 0; % ### 0
bask = parblocks*options.granularity;
if options.exhaustive
    maxtime = 100*maxtime; 
end
tic,

probabilities = zeros(1,maxmodels);

[pathstr,basname] = fileparts(fname);
solutionname = fullfile(pathstr,strcat(basname,'_solutions.dat'));
repname = fullfile(pathstr,strcat(basname,'_stemloops.dat'));
geometryname = fullfile(pathstr,strcat(basname,'_geometry.mat'));
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
    if parblocks >= restraints.search_range(1) && parblocks <= restraints.search_range(2)
        res_vec = res*ones(1,options.granularity);
        % profile on
        parfor kt = 1:options.granularity % ### parfor
            if solutions_given
                [msoln,nsoln] = size(solutions);
                [mproc,~] = size(processed);
                ksoln = zeros(1,msoln);
                psoln = 0;
                skip = true;
                for kbp = 1:msoln
                    if parblocks == solutions(kbp,1) && kt == solutions(kbp,2)
                        to_be_processed = true;
                        for kbpp = 1:mproc
                            if parblocks == processed(kbpp,1) && kt == processed(kbpp,2)
                                to_be_processed = false;
                            end
                        end
                        if to_be_processed
                            psoln = psoln + 1;
                            ksoln(psoln) = kbp;
                            skip = false;
                        end
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
                        if fulfill
                            atransmat = tmats{kt};
                            % the first refinement is performed without
                            % testing for clashes, as this is much faster
                            fprintf(1,'Refining T%i.%i\n',parblocks,kt);
                            atransmat0 = atransmat;
                            [atransmat,errvec,mprob,xld] = ...
                                refine_rba(rb,atransmat0,points,pthr,naux,auxiliary,ncore,core,links,1e6,heavy_coor,xlink_percentage,xlinks);
                            if ~sum(errvec)
                                fprintf(1,'Secondary refining T%i.%i\n',parblocks,kt);
                                [atransmat,errvec,mprob,xld] = ...
                                    refine_rba_fast(rb,atransmat,points,pthr,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,xlink_percentage,xlinks);
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
                    sfulfill = true;
                    % now compare stemloop restraints
                    if ~isempty(restraints.stemlibs)
                        modnum.block = parblocks;
                        modnum.num = k-bask;
                        modnum.model = success;
                        [sl_solutions,~,~,~,ccombi] = fit_stemloop_combinations(restraints,snum0,snum,repname,modnum);
                        if isempty(sl_solutions)
                            success = success - 1;
                            stem_fail = stem_fail + 1;
                            sfulfill = false;
                        end
                    end
                    if sfulfill
                        fid = fopen(solutionname,'at');
                        if ~isempty(ccombi) && sum(ccombi) > 0
                            fprintf(fid,'%8i%6i',parblocks,k-bask);
                            for ksl = 1:length(ccombi)
                                fprintf(fid,'%6i',ccombi(ksl));
                            end
                            fprintf(fid,'\n');
                        else
                            fprintf(fid,'%8i%6i\n',parblocks,k-bask);
                        end
                        fclose(fid);
                        soln_count = soln_count + 1;
                        t_points = cell(1,length(rb));
                        for kr = 1:length(rb)
                            t_points{kr} = affine_coor_set(points{kr},transmats{kr});
                        end
                        reference_geometry{soln_count} = t_points;
                        if restraints.search
                            success = success -1;
                        end
                        if solutions_given
                            fid = fopen(proc_name,'at');
                            fprintf(fid,'%8i%6i\n',parblocks,k-bask);
                            fclose(fid);
                        end
                        if restraints.save_rigi
                            if isfield(model,'selected')
                                model = rmfield(model,'selected');
                            end
                            model.selected = cell(1,length(model.structures{snum}));
                            for kc = 1:length(model.structures{snum})
                                model.selected{kc} = [snum kc success];
                            end
                            conformer = fullfile(pathstr,sprintf('%s_T%i_%i.pdb',basname,parblocks,k-bask));
                            wr_pdb_selected(conformer,PDBid);
                        end
                    end
                end
                if ~skip_mode && sfulfill
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
    end
    bask = bask + options.granularity;
    % save SANS_fitting_RNA chi_SANS sans_poi
    runtime=toc;
    if interactive
        % update multi plot axes
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
reference_geometry = reference_geometry(1:soln_count);
save(geometryname,'reference_geometry');
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
diagnostics.success = success;
diagnostics.probabilities = probabilities(1:success);
diagnostics.snum = snum;
diagnostics.resolution = worst_res;
diagnostics.exhaustive = options.exhaustive;
diagnostics.exhaustive_completed = bask >= trials;

