function success = rigi_flex_machine(fname,PDBid,trials,maxtime,SANS_threshold,deterministic)

global model

rfile = 'rigiflex_restraints_cheat_orig.dat';

options.granularity = 10000;

PDBid = strtrim(PDBid);

diagnostics = true;

maxat = 50000; % maximum number of heavy atoms in rigid body

maxmodels = 100;

pmodel = 0.5;
pthr = exp(-erfinv(pmodel)^2);

forgive = 0.8;
min_approach = 5; % minimal approach of two reference points [Å]
max_extension = 180; % maximum distance between any two reference points [Å]
clash_threshold = 1.5*forgive;

delete_SANS = 'poor'; % can be 'all', 'poor', or 'none', any other choice is interpreted as 'none' 
delete_SAXS = 'all'; % can be 'all', 'poor', or 'none', any other choice is interpreted as 'none' 

if ~exist('fname','var') || isempty(fname)
    fname = 'PTBP1_model';
end;

if ~exist('PDBid','var') || isempty(PDBid)
    PDBid = '7PTB';
end;

if ~exist('trials','var') || isempty(trials)
    trials = 100000;
end;

if ~exist('maxtime','var') || isempty(maxtime)
    maxtime = 2*3600;
end;

if ~exist('SANS_threshold','var') || isempty(SANS_threshold)
    SANS_threshold = 14; % 2.5^2;
end;

if ~exist('SAXS_threshold','var') || isempty(SANS_threshold)
    SAXS_threshold = 2.5^2;
end;

if ~exist('deterministic','var') || isempty(deterministic)
    deterministic = true;
end;

if deterministic
    rng(13);
else
    rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs
end;

[restraints,failed] = rd_restraints_rigiflex(rfile);

if failed
    success = -1;
    add_msg_board('ERROR: Restraint file could not be read or interpreted.');
    return
end

snum0 = model.current_structure;
if isfield(restraints,'PDB') && ~isempty(restraints.PDB)
    snum0 = resolve_address(sprintf('[%s]',restraints.PDB));
    if isempty(snum0)
        [message,snum0]=add_pdb(restraints.PDB,restraints.PDB);
        if message.error
            success = -1;
            add_msg_board('ERROR: Requested template structure not present and could not be loaded.');
            return
        end
    end
end

if isfield(restraints,'newID') && ~isempty(restraints.newID)
    PDBid = newID;
end

snum = resolve_address(sprintf('[%s]',PDBid));

if isempty(snum)
    snum = copy_structure(snum0,PDBid);
end
add_msg_board(sprintf('Result is stored in structure %i with identifier %s\n',snum,PDBid));

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
        [msg,coorh]=get_object([snum0 restraints.rb(kr).chains(kcc) 1],'xyz_heavy');
        if msg.error,
            add_msg_board(sprintf('ERROR in reading %s coordinates',mk_address([snum0 restraints.rb(kr).chains(kcc) 1])));
            add_msg_board(msg.text);
            success = -1;
            return
        end
        [mh,~] = size(coorh);
        coor_r(crpoi+1:crpoi+mh,:) = coorh;
        crpoi = crpoi + mh;
    end
    heavy_coor{kr} = coor_r(1:crpoi,:);
end
nc = kc;
all_chains = all_chains(1:nc);


model = rmfield(model,'selected');
for kc = 1:nc
    model.selected{kc} = [snum all_chains(kc) 1];
end;


% augment lower and upper bounds
for k1 = 1:8
    for k2 = k1+1:9
        if restraints.lb(k1,k2) < min_approach
            restraints.lb(k1,k2) = min_approach;
            restraints.lb(k2,k1) = min_approach;
        end;
        if restraints.ub(k1,k2) < 0.1 % the unset bounds are zero
            restraints.ub(k1,k2) = max_extension;
            restraints.ub(k2,k1) = max_extension;
        end
    end
end

[lb,ub,err]=triangle(restraints.lb,restraints.ub,true);

switch err
    case 0
        add_msg_board('Successful bound smoothing with experimental restraints.');
    case 1
        add_ms_board('ERROR: Some distance restraints are inconsistent.');
        success = -1;
        return
    case 2
        add_msg_board('ERROR: At least one bound matrix is not square.');
        success = -1;
        return
    case 3
        add_msg_board('ERROR: Upper and lower bound matrices differ in size.');
        success = -1;
        return
    otherwise
        ad_msg_board('Unspecified error in bound smoothing.');
        success = -1;
        return
end;

met_err = 0;
embed_err = 0;
bound_err = 0;
clash_err = 0;
aux_fail = 0;
link_fail = 0;
label_fail = 0;
sans_fail = 0;
saxs_fail = 0;
success = 0;
sans_poi = 0;

if trials > 50000
    diagnostics = false;
end;
if diagnostics
    prob_labels = zeros(1,trials);
    rms_trace = zeros(length(restraints.rb),trials);
end;
if isfield(restraints,'SANS')
    chi_SANS = zeros(length(restraints.SANS),50000);
end
if isfield(restraints,'SAXS')
    chi_SAXS = zeros(length(restraints.SAXS),50000);
end
final_chi2_SANS = zeros(1,maxmodels);
final_chi2_SAXS = zeros(1,maxmodels);
transmats = cell(1,3);

runtime = 0;
k = 0;
silent = false;
rb = restraints.rb;
tic,
while runtime<=maxtime && k < trials && success < maxmodels
    % k = k + 1;
    merr_vec = zeros(1,options.granularity);
    eerr_vec = zeros(1,options.granularity);
    berr_vec = zeros(1,options.granularity);
    rerr_vec = zeros(1,options.granularity);
    lerr_vec = zeros(1,options.granularity);
    cerr_vec = zeros(1,options.granularity);
    tmats = cell(options.granularity,3);
    parfor kt = 1:options.granularity
        [dmatr,err]=metrize(lb,ub,silent);
        if err==1 % metrization failed (restraints inconsistent), next trial, increment error counter
            merr_vec(kt) = 1;
        else
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
                    for kr = 1:length(rb)
                        baspoi = 3*(kr-1);
                        [~,lab_t,transmat] = rmsd_superimpose(coor1(baspoi+1:baspoi+3,:),rb(kr).ref(1:3,2:4));
                        transmats{kr} = transmat;
                        t_labels{kr} = lab_t;
                        t_points{kr} = affine_coor_set(restraints.points{kr},transmat);
                    end;
                    fulfill = true;
                    paux = 1;
                    [naux,~] = size(restraints.auxiliary);
                    for kaux = 1:naux
                        acoor1 = t_points{restraints.auxiliary(kaux,1)};
                        acoor2 = t_points{restraints.auxiliary(kaux,3)};
                        coor1b = acoor1(restraints.auxiliary(kaux,2),:);
                        coor2b = acoor2(restraints.auxiliary(kaux,4),:);
                        rsim = norm(coor1b-coor2b);
                        prob = prob_Gaussian(rsim,restraints.auxiliary(kaux,5),restraints.auxiliary(kaux,6));
                        paux = paux*prob;
                    end;
                    if paux < pthr^naux
                        fulfill = false;
                        aux_fail = aux_fail + 1;
                    end;
                    % check how well label-to-label distances are reproduced
                    plabels = 1;
                    [ncore,~] = size(restraints.core);
                    if fulfill
                        for kcore = 1:ncore
                            acoor1 = t_points{restraints.core(kcore,1)};
                            acoor2 = t_points{restraints.core(kcore,3)};
                            coor1b = acoor1(restraints.core(kcore,2),:);
                            coor2b = acoor2(restraints.core(kcore,4),:);
                            rsim = norm(coor1b-coor2b);
                            prob = prob_Gaussian(rsim,restraints.core(kcore,5),restraints.core(kcore,6));
                            plabels = plabels*prob;
                        end;
                        if diagnostics
                            prob_labels(k) = plabels;
                        end
                        if plabels < pthr^ncore
                            fulfill = false;
                            label_fail = label_fail + 1;
                        end;
                    end
                    % check for linker lengths
                    if fulfill
                        for kl = 1:length(restraints.links)
                            acoor1 = t_points{restraints.links(kl).ref_indices(1)};
                            acoor2 = t_points{restraints.links(kl).ref_indices(3)};
                            coor1b = acoor1(restraints.links(kl).ref_indices(2),:);
                            coor2b = acoor2(restraints.links(kl).ref_indices(4),:);
                            rlink = norm(coor2b - coor1b);
                            if rlink > restraints.links(kl).maxr
                                fulfill = false;
                            end
                        end
                        if ~fulfill
                            link_fail = link_fail + 1;
                        end;
                    end
                    % check for rigid-body clashes
                    if fulfill
                        clashed = false;
                        for kr1 = 1:length(restraints.rb)-1
                            hc1 = affine_coor_set(heavy_coor{kr1},transmats{kr1});
                            for kr2 = kr1+1:length(restraints.rb)
                                hc2 = affine_coor_set(heavy_coor{kr2},transmats{kr2});
                                mdist = min_dist(hc1,hc2);
                                if mdist < clash_threshold
                                    clashed = true;
                                end;
                            end;
                        end
                        if clashed
                            clash_err = clash_err + 1;
                            fulfill = false;
                        end
                    end
                end
            end
        end
    end
    % the model fulfills all restraints and should be copied
    % into the structure
    met_err = met_err + sum(merr_vec);
    embed_err = embed_err + sum(eerr_vec);
    bound_err = bound_err + sum(berr_vec);

    if fulfill
        success = success + 1;
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
                    model.structures{snum}(restraints.rb(kr).chains(kc)).residues{success} = ...
                        model.structures{2}(restraints.rb(kr).chains(kc)).residues{1};
                    model.structures{snum}(restraints.rb(kr).chains(kc)).Bfactor{success} = ...
                        model.structures{2}(restraints.rb(kr).chains(kc)).Bfactor{1};
                    model.structures{snum}(restraints.rb(kr).chains(kc)).Btensor{success} = ...
                        model.structures{2}(restraints.rb(kr).chains(kc)).Btensor{1};
                end;
            end
        end
        % now it should be compared to SANS or SAXS restraints,
        % if any
        SANS_chi = 0;
        if isfield(restraints,'SANS')
            h_SANS = figure(10+success); clf;
            hold on
            sans_vec = -ones(2,1);
            for ks = 1:length(restraints.SANS)
                model = rmfield(model,'selected');
                for kc = 1:length(restraints.SANS(ks).chains)
                    model.selected{kc} = [snum restraints.SANS(ks).chains(kc) success];
                end
                pdbfile = sprintf('t%i_%i',k,ks);
                to_be_deleted = sprintf('t%i_*.*',k);
                wr_pdb_selected(pdbfile,'SANS');
                [chi2,~,~,result,fit] = fit_SANS_by_cryson(restraints.SANS(ks).data,pdbfile,restraints.SANS(ks).illres);
                if isempty(chi2) || isnan(chi2)
                    SANS_chi = 1e6;
                    fprintf(2,'Warning: SANS fitting failed in trial %i:\n',k);
                    % fprintf(2,'%s',result);
                else
                    sans_vec(ks) = chi2;
                    plot(fit(:,1),fit(:,2));
                    hold on;
                    plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
                    SANS_chi = SANS_chi + chi2;
                end
            end
            if min(sans_vec) > 0
                sans_poi = sans_poi + 1;
                chi_SANS(:,sans_poi) = sans_vec;
            end;
            chi2 = SANS_chi/length(restraints.SANS);
            if chi2 > SANS_threshold
                fprintf(1,'SANS chi^2 of %4.2f exceeded threshold of %4.2f in trial %i.\n',chi2,SANS_threshold,k);
                close(h_SANS);
                success = success - 1;
                sans_fail = sans_fail + 1;
                fulfill = false;
                if strcmp(delete_SANS,'all') || strcmp(delete_SANS,'poor')
                    delete(strcat(to_be_deleted,'*.*'));
                end;
            else
                title(sprintf('SANS fit for model %i (chi^2 = %4.2f)',success,chi2));
                drawnow
                if strcmp(delete_SANS,'all')
                    delete(strcat(to_be_deleted,'*.*'));
                end;
                final_chi2_SANS(success) = chi2;
            end
        end
        SAXS_chi = 0;
        if isfield(restraints,'SAXS') && fulfill
            h_SAXS = figure(1010+success); clf;
            hold on
            for ks = 1:length(restraints.SAXS)
                model = rmfield(model,'selected');
                for kc = 1:length(restraints.SAXS(ks).chains)
                    model.selected{kc} = [snum restraints.SAXS(ks).chains(kc) success];
                    pdbfile = sprintf('t%i_%i',k,ks);
                    wr_pdb_selected(pdbfile,'SAXS');
                    [chi2,~,~,result,fit] = fit_SAXS_by_crysol(restraints.SAXS(ks).data,pdbfile);
                    if isempty(chi2)
                        fprintf(2,'Warning: SAXS fitting failed in trial %i:\n',k);
                        fprintf(2,'%s',result);
                        delete(strcat(pdbfile,'*.*'));
                        success = success - 1;
                        saxs_fail = saxs_fail + 1;
                    else
                        chi_SAXS(ks,success) = chi2;
                        plot(fit(:,1),fit(:,2));
                        hold on;
                        plot(fit(:,1),fit(:,3),':','Color',[0.75,0,0]);
                        SAXS_chi = SAXS_chi + chi2;
                    end
                end
            end
            chi2 = SAXS_chi/length(restraints.SAXS);
            if chi2 > SAXS_threshold
                close(h_SAXS);
                success = success - 1;
                saxs_fail = saxs_fail + 1;
                if strcmp(delete_SAXS,'all') || strcmp(delete_SAXS,'poor')
                    delete(strcat(pdbfile,'*.*'));
                end;
            else
                title(sprintf('SAXS fit for model %i (chi^2 = %4.2f)',success,chi2));
                drawnow
                if strcmp(delete_SAXS,'all')
                    delete(strcat(pdbfile,'*.*'));
                end;
            end
        end
    end;
    runtime=toc;
end;
toc,

trials = k;

final_chi2_SANS = final_chi2_SANS(1:success);
final_chi2_SAXS = final_chi2_SAXS(1:success);

if isfield(restraints, 'SANS')
    chi_SANS = chi_SANS(:,1:sans_poi);
    figure(500); clf;
    plot(chi_SANS(1,:),'bo');
    hold on;
    plot(chi_SANS(2,:),'r*');
end;
save all_chi_SANS chi_SANS final_chi2_SANS
if isfield(restraints, 'SAXS')
    chi_SAXS = chi_SAXS(:,1:success);
    figure(1500); clf;
    plot(chi_SAXS,'o');
end;

if success > 0
    model = rmfield(model,'selected');
    spoi = 0;
    for kc = 1:7
        for km = 1:success
            spoi = spoi + 1;
            model.selected{spoi} = [2 kc km];
        end;
    end;
    message = wr_pdb_selected(fname,PDBid);
    if message.error
        fprintf(2,'Warning: Model could not be automatically saved. %s\n',message.text);
    end
    [pathstr,basname] = fileparts(fname);
    scriptname = fullfile(pathstr,strcat(basname,'.mmm'));
    fid = fopen(scriptname,'wt');
    for kr = 1:length(restraints.rb)
        for kc = 1:length(restraints.rb(kr).chains)
            [~,ctag] = mk_address_parts([snum,restraints.rb(kr).chains(kc)]);
            fprintf(fid,'show [%s](%s){:} ribbon\n',PDBid,ctag);
            fprintf(fid,'color [%s](%s){:} %s\n',PDBid,ctag,restraints.color{restraints.rb(kr).chains(kc)});
        end
    end
    prob = 1./(final_chi2_SANS.^2+final_chi2_SAXS.^2);
    prob_norm = 2*max(prob); 
    for km = 1:success        
        fprintf(fid,'transparency [%s](:){%i} %5.3f\n',PDBid,km,prob(km)/prob_norm);
    end
    fclose(fid);
end

left_trials = trials;
fprintf(1,'%5.2f%% of %i trials were successful. %i models were obtained.\n',100*success/trials,trials,success);
fprintf(1,'%5.2f%% of all trials failed in metrization.\n',100*met_err/trials);
left_trials = left_trials - met_err;
fprintf(1,'%5.2f%% of all trials (%i) failed in embedding (%5.2f%% of still active trials).\n',100*embed_err/trials,embed_err,100*embed_err/left_trials);
left_trials = left_trials - embed_err;
fprintf(1,'%5.2f%% of all trials (%i) failed on bound refinement (%5.2f%% of still active trials).\n',100*bound_err/trials,bound_err,100*bound_err/left_trials);
left_trials = left_trials - bound_err;
fprintf(1,'%5.2f%% of all trials (%i) failed on auxiliary restraints (%5.2f%% of still active trials).\n',100*aux_fail/trials,aux_fail,100*aux_fail/left_trials);
left_trials = left_trials - aux_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed on label restraints (%5.2f%% of still active trials).\n',100*label_fail/trials,label_fail,100*label_fail/left_trials);
left_trials = left_trials - label_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed on linker length (%5.2f%% of still active trials).\n',100*link_fail/trials,link_fail,100*link_fail/left_trials);
left_trials = left_trials - link_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed in clash tests (%5.2f%% of still active trials).\n',100*clash_err/trials,clash_err,100*clash_err/left_trials);
left_trials = left_trials - clash_err;
fprintf(1,'%5.2f%% of all trials (%i) failed due to poor SANS fit (%5.2f%% of still active trials).\n',100*sans_fail/trials,sans_fail,100*sans_fail/left_trials);
left_trials = left_trials - sans_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed due to poor SAXS fit (%5.2f%% of still active trials).\n',100*saxs_fail/trials,saxs_fail,100*saxs_fail/left_trials);
left_trials = left_trials - saxs_fail;
if left_trials ~= success
    fprintf(2,'Trial dissipation. Expected success: %i. Found success %i.\n',left_trials,success);
end

if diagnostics
%     figure(777); clf;
%     plot(rna_dist,'k');
%     hold on;
%     plot([1,trials],[restrain_rna.r,restrain_rna.r],'g');
%     plot([1,trials],[restrain_rna.r-restrain_rna.sigr,restrain_rna.r-restrain_rna.sigr],'r:');
%     plot([1,trials],[restrain_rna.r+restrain_rna.sigr,restrain_rna.r+restrain_rna.sigr],'r:');
    
    figure(333); clf;
    plot(rms_trace(1,:));
    hold on;
    plot(rms_trace(2,:));
    plot(rms_trace(3,:));
    
    figure(777); clf;
    plot(prob_labels,'k');
    hold on;
    plot([1,trials],[pthr^9,pthr^9],'r:');
    
end;

