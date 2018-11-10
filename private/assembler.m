function diagnostics = assembler(all_flex_models,restraints,options,handles)
% Assembles complete models from rigid-body arrangements and flexible
% linker conformations
% the model with all rigid-body arrangements must be the current model
%
% all_flex_models   array of structure numbers for all flexible sections
%                   (second index) of all rigid-body arrangements (first
%                   index)
% restraints        original restraint definition, needed for checking
%                   SANS, SAXS, crosslinks
% options           .maxtime (maximum time)
%                   .fname_bas (basis file name of the modeling run)
% handles           handles of the GUI, optional, needed for interactive
%                   mode
%
% G. Jeschke, 2.5.2017

global model

maxatoms = 50000;
maxtime = 3600*options.max_time;

forgive = 0.8; % scaling factor for van-der-Waals radii in clash tests
clash_threshold = 1.5*forgive; % a uniform van-der-Waals radius of 1.5 Å is assumed for heavy atoms

diagnostics.success = 0;

% if handles of the GUI figure are supplied, the engine will update the MMM
% GUI during the computation, else this is considered as a run on a remote
% server that does not have access to the GUI
if exist('handles','var')
    interactive = true;
else
    interactive = false;
end

% all_flex_models = all_flex_models(:,1:2);

diagnostics.report_name = strcat(options.fname_bas,'_assembler_report.txt');
fidr = fopen(diagnostics.report_name,'wt');
fprintf(fidr,'--- Assembling models from rigid-body arangements and flexible sections for structure %s ---\n',options.fname_bas);

[arrangements,sections] = size(all_flex_models);

adr_rigid = mk_address(model.current_structure);
nummod = zeros(arrangements,sections);
tic;
for ka = 1:arrangements
    combi = 1;
    for ks = 1:sections
        secnum = all_flex_models(ka,ks);
        if secnum == 0
            combi = 0;
        else
            nummod(ka,ks) = length(model.structures{secnum}(1).residues);
            combi = combi*nummod(ka,ks);
        end
    end
    if combi == 0 % at least for one section no models were found
        continue;
    end;
    fprintf(fidr,'\nRigid-body arrangement %i has %i potential combinations of flexible sections.\n',ka,combi);
    if interactive
        handles.uipanel_runtime.Title = sprintf('Assembling flexible sections for rigid-body arrangement %i',ka);
        handles.text_max_trials = sprintf('%i',combi);
    end
    adr_arr = sprintf('%s(:){%i}',adr_rigid,ka);
    [~,rigid_all_xyz] = get_object(adr_arr,'xyz');
    rigid_xyz = zeros(maxatoms,3);
    poi = 0;
    for km = 1:length(rigid_all_xyz)
        [m,~] = size(rigid_all_xyz{km});
        rigid_xyz(poi+1:poi+m,:) = rigid_all_xyz{km};
        poi = poi + m;
    end;
    rigid_xyz = rigid_xyz(1:poi,:);
    all_combi = zeros(1,2*sections);
    for ks = 1:sections
        adr_sec{ks} = mk_address(all_flex_models(ka,ks));
        [mc,~] = size(all_combi);
        ns = length(model.structures{all_flex_models(ka,ks)}(1).residues);
        all_combi = repmat(all_combi,ns,1);
        modnumvec = 1:length(model.structures{all_flex_models(ka,ks)}(1).residues);
        secnumvec = all_flex_models(ka,ks)*ones(mc,1);
        for kc = 1:ns
            basnum = (kc-1)*mc;
            all_combi(basnum+1:basnum+mc,2*ks-1) = secnumvec;
            all_combi(basnum+1:basnum+mc,2*ks) = kc*ones(mc,1);
        end
    end;
    [mc,~] =size(all_combi);
    min_dist = 1e6*ones(1,mc);
    for kc = 1:mc % test all combinations for clashes
        runtime = toc;
        if runtime > maxtime
            break
        end
        for ks1 = 1:sections-1
            s1 = [all_combi(kc,2*ks1-1) 1 all_combi(kc,2*ks1)]; % indices for first section in pair
            [~,xyz1] = get_chain_model(s1,'xyz_heavy');
            [m1,~] = size(xyz1); % get sizes of the coordinates arrays
            for ks2 = ks1+1:sections
                s2 = [all_combi(kc,2*ks2-1) 1 all_combi(kc,2*ks2)]; % indices for second section in pair
                [~,xyz2] = get_chain_model(s2,'xyz_heavy');
                [m2,~] = size(xyz2);
                a2 = repmat(sum(xyz1.^2,2),1,m2);
                b2 = repmat(sum(xyz2.^2,2),1,m1).';
                pair_dist = sqrt(abs(a2 + b2 - 2*xyz1*xyz2.'));
                min_dist_combi = min(min(pair_dist));
                if min_dist_combi < min_dist(kc)
                    min_dist(kc) = min_dist_combi;
                end
           end
        end
    end
    scores = zeros(1,mc);
    for kc = 1:mc
        if min_dist(kc) < clash_threshold
            fprintf(fidr,'Arrangement %i, section combination %i clashes at minimum distance of %4.1f Å.\n',ka,kc,min_dist(kc));
            scores(kc) = 1e6;
        end
    end
    if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
        fprintf(fidr,'\n    --- SANS fitting ---\n');
        SANS_scores= zeros(1,mc);
        SANS_curves{1} = [];
        for kc = 1:mc
            runtime = toc;
            if runtime > maxtime
                break
            end
            if min_dist(kc) > clash_threshold
                if isfield(model, 'selected')
                    model = rmfield(model,'selected');
                end
                ksel = 0;
                SANS_chi = 0;
                SANS_scores(kc) = 0;
                for ks = 1:length(restraints.SANS)
                    % ### needs to be fixed for RNA-connected model ###
                    for kch = 1:length(restraints.SANS(ks).chains)
                        ksel = ksel + 1;
                        if restraints.SANS(ks).chains(kch) > 0 % chain in the rigid-body arrangement
                            model.selected{ksel} = [model.current_structure restraints.SANS(ks).chains(kch) ka];
                        else % flexible segment
                            bas = 2*(abs(restraints.SANS(ks).chains(kch))-1);
                            secstruct = all_combi(kc,bas+1);
                            secmod = all_combi(kc,bas+2);
                            model.selected{ksel} = [secstruct 1 secmod]; % the model that specifies this segment has only one chain
                        end
                    end
                    % model.selected{1} = [all_combi(kc,1) 1 all_combi(kc,2)];
                    pdbfile = sprintf('c%i_%i',ka,kc);
                    to_be_deleted = sprintf('c%i_*.*',ka);
                    wr_pdb_selected(pdbfile,'SANS');
                    [chi2,~,~,result,fit] = fit_SANS_by_cryson(restraints.SANS(ks).data,pdbfile,restraints.SANS(ks).illres);
                    if isempty(chi2) || isnan(chi2)
                        SANS_chi = 1e6;
                        if interactive
                            fprintf(2,'Warning: SANS fitting failed for combination %i:\n',kc);
                            fprintf(2,'%s',result);
                        end;
                        fprintf(fidr,'SANS fitting of curve %i failed in arrangement %i for combination %i.\n',ks,ka,kc);
                    else
                        SANS_curves{ks} = fit;
                        SANS_chi = SANS_chi + chi2;
                    end
                    delete(to_be_deleted);
                end;
                chi2 = SANS_chi/length(restraints.SANS);
                SANS_scores(kc) = SANS_scores(kc)+chi2;
                if chi2 < 1e5
                    fprintf(fidr,'SANS fit for arrangement %i with flexible section combination %i has chi^2 = %5.2f.\n',ka,kc,chi2);
                    if interactive && options.display_SANS_fit
                        % update multi plot axes
                        axes(handles.axes_multi_plot);
                        cla;
                        hold on;
                        for ks = 1:length(restraints.SANS)
                            fit = SANS_curves{ks};
                            plot(fit(:,1),fit(:,2));
                            plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
                        end;
                        title(sprintf('SANS fit for arrangement %i combi %i (chi^2 = %4.2f)',ka,kc,chi2));
                        drawnow
                    end
                end
            end
        end;
        all_SANS{ka} = SANS_curves;
        scores = scores + SANS_scores;
    end
    if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS)
        fprintf(fidr,'\n    --- SAXS fitting ---\n');
        SAXS_scores= zeros(1,mc);
        SAXS_curves{1} = [];
        for kc = 1:mc
            runtime = toc;
            if runtime > maxtime
                break
            end
            if min_dist(kc) > clash_threshold
                if isfield(model, 'selected')
                    model = rmfield(model,'selected');
                end
                ksel = 0;
                SAXS_chi = 0;
                SAXS_scores(kc) = 0;
                for ks = 1:length(restraints.SAXS)
                    for kch = 1:length(restraints.SAXS(ks).chains)
                        ksel = ksel + 1;
                        if restraints.SAXS(ks).chains(kch) > 0 % chain in the rigid-body arrangement
                            model.selected{ksel} = [model.current_structure restraints.SAXS(ks).chains(kch) ka];
                        else % flexible segment
                            bas = 2*(abs(restraints.SAXS(ks).chains(kch))-1);
                            secstruct = all_combi(kc,bas+1);
                            secmod = all_combi(kc,bas+2);
                            model.selected{ksel} = [secstruct 1 secmod]; % the model that specifies this segment has only one chain
                        end
                    end
                    pdbfile = sprintf('c%i_%i',ka,kc);
                    to_be_deleted = sprintf('c%i_*.*',ka);
                    wr_pdb_selected(pdbfile,'SAXS');
                    [chi2,~,~,result,fit] = fit_SAXS_by_crysol(restraints.SAXS(ks).data,pdbfile);
                    if isempty(chi2) || isnan(chi2)
                        SAXS_chi = 1e6;
                        if interactive
                            fprintf(2,'Warning: SAXS fitting failed for combination %i:\n',kc);
                            fprintf(2,'%s',result);
                        end
                        fprintf(fidr,'SAXS fitting of curve %i failed in arrangement %i for combination %i.\n',ks,ka,kc);
                    else
                        SAXS_curves{ks} = fit;
                        SAXS_chi = SAXS_chi + chi2;
                    end
                    delete(to_be_deleted);
                end
                chi2 = SAXS_chi/length(restraints.SAXS);
                SAXS_scores(kc) = SAXS_scores(kc)+chi2;
                if chi2 < 1e5
                    fprintf(fidr,'SAXS fit for arrangement %i with flexible section combination %i has chi^2 = %5.2f.\n',ka,kc,chi2);
                    if interactive && options.display_SANS_fit
                        % update multi plot axes
                        axes(handles.axes_multi_plot);
                        cla;
                        hold on;
                        for ks = 1:length(restraints.SAXS)
                            fit = SAXS_curves{ks};
                            plot(fit(:,1),fit(:,2));
                            plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
                        end
                        title(sprintf('SAXS fit for arrangement %i combi %i (chi^2 = %4.2f)',ka,kc,chi2));
                        drawnow
                    end
                end
            end
        end
        all_SAXS{ka} = SAXS_curves;
        scores = scores + SAXS_scores;
    end
    full_scores{ka} = scores;
    full_combi{ka} = all_combi;
    runtime = toc;
    if runtime > maxtime
        break
    end
end

if ~exist('full_scores','var')
    diagnostics.success = 0;
    fprintf(fidr,'\n### ERROR: No valid model could be assembled. Aborting. ###\n');
    add_msg_board('ERROR: Assembler failed. No valid model.');
    fclose(fidr);
    return
end

diagnostics.scores = full_scores;
diagnostics.combi = full_combi;
if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
    diagnostics.SANS_curves = all_SANS;
end
if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS)
    diagnostics.SAXS_curves = all_SAXS;
end


% combine and sort scores
num_models = 0;
for ka = 1:length(full_scores)
    num_models = num_models + length(full_scores{ka});
end
all_scores = zeros(1,num_models);
all_combi = zeros(num_models,1+2*sections);
poi = 0;
for ka = 1:length(full_scores)
    nma = length(full_scores{ka});
    if nma > 0
        all_scores(poi+1:poi+nma) = full_scores{ka};
        all_combi(poi+1:poi+nma,1) = ka*ones(nma,1);
        all_combi(poi+1:poi+nma,2:end) = full_combi{ka};
        poi = poi + nma;
    end
end    
[sorted_scores,score_indices] = sort(all_scores);
listed_models = 5*restraints.ensemble;
valid_models = length(find(sorted_scores < 1e5));
if listed_models > valid_models
    listed_models = valid_models;
end
fprintf(fidr,'\n--- List of the %i best scoring models ---\n\n',listed_models);
for km = 1:listed_models
    fprintf(fidr,'Arrangement %i with flexible section models ',all_combi(score_indices(km),1));
    for ks = 1:sections
        fprintf(fidr,'%i',all_combi(score_indices(km),2*(ks-1)+3));
        if ks < sections
            fprintf(fidr,', ');
        end
    end
    fprintf(fidr,' has a score of %5.2f\n',sorted_scores(km));
end

if valid_models < restraints.ensemble
    fprintf(fidr,'\n### Warning! Only %i valid models while an ensemble size of %i models was requested. ###\n',listed_models,restraints.ensemble);
end

% Save models

nummod = restraints.ensemble;
if length(sorted_scores) < nummod
    nummod = length(sorted_scores);
end

for km = 1:nummod
    if sorted_scores(km) < 1e5
        pdbfile = strcat(options.fname_bas,sprintf('_model_%i',km));
        if isfield(model, 'selected')
            model = rmfield(model,'selected');
            ksel = 0;
        end
        ka = all_combi(score_indices(km),1);
        for kch = 1:length(model.structures{model.current_structure}) % select chain models of the rigid-body arrangement
            ksel = ksel + 1;
            model.selected{ksel} = [model.current_structure kch ka];
        end
        for ks = 1:sections
            ksel = ksel + 1;
            bas = 2*(ks-1)+1;
            secstruct = all_combi(score_indices(km),bas+1);
            secmod = all_combi(score_indices(km),bas+2);
            model.selected{ksel} = [secstruct 1 secmod]; % the model that specifies this segment has only one chain
        end
        wr_pdb_selected(pdbfile,['M' mk_modnum(km)],sorted_scores(km),true);
    end
end

info.idCode = 'PTBX';
info.class = 'MMM model';
info.depDate = date;
info.title = 'RigiFlex model';
info.chain_tags = ':A:';
info.chain_ids = 1;
info.center = [0,0,0];
info.B_range = [0,20];
info.authors = 'MMM';
info.molecule = 'chimera';
info.organism = 'Dragon musculus';
info.remarks = [4 5 6];
info.SSbonds = [];
info.Modeller_obj = [];
info.Modeller_sid = [];
info.rotamers = false;
info.missing = {};
info.site_tags = ':';
info.sites = [];
info.keywords = '';
info.metal = [];
info.resolution = [];
info.references = [];
info.insertions = [];
info.alternate = false;

pctag = 'A';

snum_rigid = model.current_structure;
adr_rigid = mk_address(snum_rigid);
included = zeros(1,length(model.structures{snum_rigid}));

rna_binding_motifs = ':';
for kRNA = 1:length(restraints.RNA)
    for kbind = 1:length(restraints.RNA(kRNA).bind)
        ctag = restraints.RNA(kRNA).bind(kbind).anchora;
        [cind,msg] = resolve_address(sprintf('%s%s',adr_rigid,ctag));
        if ~isempty(cind) && ~msg.error
            included(cind(2)) = 1;
        end
        pa = strfind(ctag,'(');
        pe = strfind(ctag,')');
        rna_binding_motifs = sprintf('%s%s:',rna_binding_motifs,ctag(pa+1:pe-1));
    end
end
for km = 1:nummod
    ka = all_combi(score_indices(km),1);
    clear indices
    found = true;
    cch = 0;
    while found
        cch = cch+1;
        ctag = id2tag(cch,restraints.peptide_tags{1});
        rnaid = tag2id(ctag,rna_binding_motifs);
        if isempty(ctag) || ~isempty(rnaid)
            found = false;
        else
            [cind,msg] = resolve_address(sprintf('%s(%s)',adr_rigid,ctag));
            if ~isempty(cind) && ~msg.error % this is a rigid-body chain (but RNA sections must be removed)
                indices{cch} = [snum_rigid cind(2) ka];
                included(cind(2)) = 1;
            else
                flexnum = str2double(ctag);
                bas = 2*(flexnum-1)+1;
                secstruct = all_combi(score_indices(km),bas+1);
                secmod = all_combi(score_indices(km),bas+2);
                indices{cch} = [secstruct 1 secmod];
            end            
        end
    end
    
    chain = combine_chains(indices,pctag);
    
    info.atoms = chain.atoms{1};
    info.residues = length(chain.residues{1});
    
    if km == 1
        [~,snum] = add_pdb(chain,info);
        tba = 1;
        
        % add the chains that are not in contiguous segments
        for kc = 1:length(included)
            if ~included(kc)
                tba = tba + 1;
                names = fieldnames(model.structures{snum_rigid}(kc));
                for kf = 1:length(names)
                    model.structures{snum}(tba).(names{kf}) = model.structures{snum_rigid}(kc).(names{kf});
                end
                iso =  model.structures{snum_rigid}(kc).isotopes;
                [cat,~] = size(iso);
                model.structures{snum}(tba).atoms{1} = cat;
                model.structures{snum}(tba).residues{1} = model.structures{snum_rigid}(kc).residues{ka};
                model.structures{snum}(tba).xyz{1} = model.structures{snum_rigid}(kc).xyz{ka};
                model.structures{snum}(tba).Bfactor{1} = model.structures{snum_rigid}(kc).Bfactor{ka};
                model.structures{snum}(tba).Btensor{1} = model.structures{snum_rigid}(kc).Btensor{ka};
            end
        end
    else
        copy_structure(snum,'+mod',[],km,snum); % make a copy of the first model
        model.structures{snum}(1).atoms{km} = chain.atoms{1};
        model.structures{snum}(1).residues{km} = chain.residues{1};
        model.structures{snum}(1).xyz{km} = chain.xyz{1};
        model.structures{snum}(1).Bfactor{km} = chain.Bfactor{1};
        model.structures{snum}(1).Btensor{km} = chain.Btensor{1};
        % replace cocordinates of chains that are not in contiguous
        % segments
        tba = 1;
        for kc = 1:length(included)
            if ~included(kc)
                tba = tba + 1;
                iso =  model.structures{snum_rigid}(kc).isotopes;
                [cat,~] = size(iso);
                model.structures{snum}(tba).atoms{km} = cat;
                model.structures{snum}(tba).residues{km} = model.structures{snum_rigid}(kc).residues{ka};
                model.structures{snum}(tba).xyz{km} = model.structures{snum_rigid}(kc).xyz{ka};
                model.structures{snum}(tba).Bfactor{km} = model.structures{snum_rigid}(kc).Bfactor{ka};
                model.structures{snum}(tba).Btensor{km} = model.structures{snum_rigid}(kc).Btensor{ka};
            end
        end
    end
end
model.current_structure = snum;
pdbfile = strcat(options.fname_bas,'_ensemble');
wr_pdb(pdbfile,restraints.newID);

fclose(fidr);

function modnum = mk_modnum(km)
modnum = sprintf('%i',km);
while length(modnum) < 3
    modnum = ['0' modnum];
end

