function PTBP1_quality_check(fname0,restraint_file,list_file)

global model
global general

interactive = true;
unprocessed = true;

if exist('fname0','var')
    [~,snum0] = add_pdb(fname0);
    fname = sprintf('%s_q.pdb',fname0);
    options.fname = fname;
    [~,name,~] = fileparts(fname);
else
    snum0 = model.current_structure;
end
adr0 = mk_address(snum0);


score_SAS = 0;
score_DEER = 1;
n_DEER = 0;
n_SAS = 0;

[restraints,failed] = rd_restraints_rigiflex(restraint_file,unprocessed); % 181105

if failed
    add_msg_board('Restraint file could not be read.');
    return
end

if ~exist('fname0','var')
    [filename, pathname] = uiputfile(['PTB1' '.pdb'], 'Save final model as PDB');
    if isequal(filename,0) || isequal(pathname,0)
        add_msg_board('Quality check cancelled by user');
        return
    else
        reset_user_paths(pathname);
        general.pdb_files=pathname;
        fname = fullfile(pathname, filename);
        [~,name,~] = fileparts(fname);
        options.fname = fname;
    end
end

fid = fopen(sprintf('%s_diagnosis.dat',name),'wt');
if fid == -1
    add_msg_board('Quality report file could not be written.');
    return
end

oname = sprintf('%s_diagnosis_data.mat',name); 

model.current_structure = snum0;

% sadr = mk_address(snum0);

hfig=gcf;
set(hfig,'Pointer','watch');

fprintf(fid,'DEER restraints for model %s\n\n',name);
fprintf(fid,'>>> core and RNA\n\n');
for k = 1:length(restraints.DEER)
    cpoi = strfind(restraints.DEER(k).adr1,'(');
    switch restraints.DEER(k).adr1(cpoi+1)
        case {'A','C','E'}
            restraints.DEER(k).adr1(cpoi+1) = 'B';
        otherwise
            restraints.DEER(k).adr1(cpoi+1) = 'A';
    end
    cpoi = strfind(restraints.DEER(k).adr2,'(');
    switch restraints.DEER(k).adr2(cpoi+1)
        case {'A','C','E'}
            restraints.DEER(k).adr2(cpoi+1) = 'B';
        otherwise
            restraints.DEER(k).adr2(cpoi+1) = 'A';
    end
    [rax,distr] = mk_distance_distribution(restraints.DEER(k).adr1,restraints.DEER(k).adr2,restraints.DEER(k).label);
    if isfield(restraints.DEER(k),'file') && ~isempty(restraints.DEER(k).file) && ~isempty(rax)
        fit_options.depth = [0.005,0.65];
        fit_options.dim = 3;
        fit_options.kdec = [];
        fit_options.cutoff = 0.9;
        [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax,distr,strcat('deer\',restraints.DEER(k).file),fit_options);
        figure(1000+k); clf;
        plot(texp,vexp,'k');
        hold on;
        plot(fit_options.cutoff*[max(texp),max(texp)],[min(vexp),max(vexp)],'b');
        plot(texp,deer,'Color',[0.75,0,0]);
        plot(texp,bckg,'Color',[0,0.6,0]);
        if restraints.DEER(k).r ~=0 &&  restraints.DEER(k).sigr ~=0
            ftype = 'fit';
        else
            ftype = 'control';
        end
        title(sprintf('%s (%s): %s-%s. rmsd: %6.4f',name,ftype,restraints.DEER(k).adr1,restraints.DEER(k).adr2,param.rmsd));
        fprintf(fid,'%s-%s rmsd: %6.4f, mod. depth: %6.3f, kdec: %6.4f\n',...
            restraints.DEER(k).adr1,restraints.DEER(k).adr2,...
            param.rmsd,param.depth,param.kdec);
    end
    restraints.DEER(k).rax = rax;
    restraints.DEER(k).distr = distr;
    if isempty(rax)
        figure(k); clf;
        title(sprintf('%s: %s-%s. Labeling failure',name,restraints.DEER(k).adr1,restraints.DEER(k).adr2));
        fprintf(fid,'%s-%s. Labeling failure\n',restraints.DEER(k).adr1,restraints.DEER(k).adr2);
        continue
    end
    rax = 10*rax;
    restraints.DEER(k).rax = rax;
    restraints.DEER(k).r_model = sum(restraints.DEER(k).distr.*rax);
    restraints.DEER(k).sigr_model = sqrt(sum(restraints.DEER(k).distr.*(rax-restraints.DEER(k).r_model).^2));
    fprintf(fid,'%s-%s requested: %4.1f (%4.1f) Å found: %4.1f (%4.1f) Å; ',...
        restraints.DEER(k).adr1,restraints.DEER(k).adr2,...
        restraints.DEER(k).r,restraints.DEER(k).sigr,...
        restraints.DEER(k).r_model,restraints.DEER(k).sigr_model);
    figure(k); clf;
    title(sprintf('%s: %s-%s',name,restraints.DEER(k).adr1,restraints.DEER(k).adr2));
    hold on
    plot(rax,distr,'Color',[0.7,0,0]);
    xlabel('r [Å]');
    ylabel('P(r)');
    hold on
    argr = (restraints.DEER(k).r-rax)/(sqrt(2)*restraints.DEER(k).sigr);
    distr_sim = exp(-argr.^2);
    distr_sim = distr_sim/sum(distr_sim);
    restraints.DEER(k).overlap = sum(min([distr_sim;distr]));
    if ~isnan(restraints.DEER(k).overlap)
        score_DEER = score_DEER*restraints.DEER(k).overlap;
        n_DEER = n_DEER + 1;
    end
    fprintf(fid,'Overlap: %5.3f\n',restraints.DEER(k).overlap);
    plot(rax,distr_sim,'Color',[0,0.6,0]);
    if isfield(restraints.DEER(k),'file') && ~isempty(restraints.DEER(k).file)
        dfname = strcat(restraints.DEER(k).file,'_distr.dat');
        dfname = strcat('deer_analysis\',dfname);
        Pdata = load(dfname);
        rexp = Pdata(:,1).';
        distr_exp = Pdata(:,2).';
        [~,rpoi] = min(abs(rexp-10));
        rexp = rexp(1:rpoi)*10;
        distr_exp = distr_exp(1:rpoi);
        distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
        distr_exp = distr_exp/sum(distr_exp);
        plot(rax,distr_exp,'k');
    end
end

core = length(restraints.DEER);

fprintf(fid,'\n>>> flexible peptide linkers\n\n');
for kl = 1:length(restraints.pflex)
    for k = 1:length(restraints.pflex(kl).DEER)
        cpoi = strfind(restraints.pflex(kl).DEER(k).adr1,'(');
        switch restraints.pflex(kl).DEER(k).adr1(cpoi+1)
            case {'A','C','E'}
                restraints.pflex(kl).DEER(k).adr1(cpoi+1) = 'B';
            otherwise
                restraints.pflex(kl).DEER(k).adr1(cpoi+1) = 'A';
        end
        cpoi = strfind(restraints.pflex(kl).DEER(k).adr2,'(');
        switch restraints.pflex(kl).DEER(k).adr2(cpoi+1)
            case {'A','C','E'}
                restraints.pflex(kl).DEER(k).adr2(cpoi+1) = 'B';
            otherwise
                restraints.pflex(kl).DEER(k).adr2(cpoi+1) = 'A';
        end
        label = [restraints.pflex(kl).DEER(k).label1 '|' restraints.pflex(kl).DEER(k).label2];
        [rax,distr] = mk_distance_distribution(restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2,label);
        if isfield(restraints.pflex(kl).DEER(k),'file') && ~isempty(restraints.pflex(kl).DEER(k).file) && ~isempty(rax)
            fit_options.depth = [0.005,0.65];
            fit_options.dim = 3;
            fit_options.kdec = [];
            fit_options.cutoff = 0.9;
            [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax,distr,strcat('deer\',restraints.pflex(kl).DEER(k).file),fit_options);
            figure(1000+core+k); clf;
            plot(texp,vexp,'k');
            hold on;
            plot(fit_options.cutoff*[max(texp),max(texp)],[min(vexp),max(vexp)],'b');
            plot(texp,deer,'Color',[0.75,0,0]);
            plot(texp,bckg,'Color',[0,0.6,0]);
            if restraints.DEER(k).r ~=0 &&  restraints.DEER(k).sigr ~=0
                ftype = 'fit';
            else
                ftype = 'control';
            end
            title(sprintf('%s (%s): %s-%s. rmsd: %6.4f',name,ftype,restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2,param.rmsd));
            fprintf(fid,'%s-%s rmsd: %6.4f, mod. depth: %6.3f, kdec: %6.4f\n',...
                restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2,...
                param.rmsd,param.depth,param.kdec);
        end
        restraints.pflex(kl).DEER(k).rax = rax;
        restraints.pflex(kl).DEER(k).distr = distr;
        if isempty(rax)
            figure(core+k); clf;
            title(sprintf('%s: %s-%s. Labeling failure',name,restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2));
            fprintf(fid,'%s-%s. Labeling failure\n',restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2);
            continue
        end
        rax = 10*rax;
        restraints.pflex(kl).DEER(k).rax = rax;
        restraints.pflex(kl).DEER(k).r_model = sum(restraints.pflex(kl).DEER(k).distr.*rax);
        restraints.pflex(kl).DEER(k).sigr_model = sqrt(sum(restraints.pflex(kl).DEER(k).distr.*(rax-restraints.pflex(kl).DEER(k).r_model).^2));
        fprintf(fid,'%s-%s requested: %4.1f (%4.1f) Å found: %4.1f (%4.1f) Å; ',...
            restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2,...
            restraints.pflex(kl).DEER(k).r,restraints.pflex(kl).DEER(k).sigr,...
            restraints.pflex(kl).DEER(k).r_model,restraints.pflex(kl).DEER(k).sigr_model);
        figure(core+k); clf;
        title(sprintf('%s: %s-%s',name,restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2));
        hold on
        plot(rax,distr,'Color',[0.7,0,0]);
        xlabel('r [Å]');
        ylabel('P(r)');
        hold on
        argr = (restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*restraints.pflex(kl).DEER(k).sigr);
        distr_sim = exp(-argr.^2);
        distr_sim = distr_sim/sum(distr_sim);
        restraints.pflex(kl).DEER(k).overlap = sum(min([distr_sim;distr]));
        if ~isnan(restraints.DEER(k).overlap)
            score_DEER = score_DEER*restraints.DEER(k).overlap;
            n_DEER = n_DEER + 1;
        end
        fprintf(fid,'Overlap: %5.3f\n',restraints.pflex(kl).DEER(k).overlap);
        plot(rax,distr_sim,'Color',[0,0.6,0]);
        if isfield(restraints.pflex(kl).DEER(k),'file') && ~isempty(restraints.pflex(kl).DEER(k).file)
            dfname = strcat(restraints.pflex(kl).DEER(k).file,'_distr.dat');
            dfname = strcat('deer_analysis\',dfname);
            Pdata = load(dfname);
            rexp = Pdata(:,1).';
            distr_exp = Pdata(:,2).';
            [~,rpoi] = min(abs(rexp-10));
            rexp = rexp(1:rpoi)*10;
            distr_exp = distr_exp(1:rpoi);
            distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
            distr_exp = distr_exp/sum(distr_exp);
            plot(rax,distr_exp,'k');
        end
    end
end

score_DEER = 1 - score_DEER^(1/n_DEER);

% Now test for SANS restraints

SANS_chi = 0;
if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
    add_msg_board('Fitting SANS restraints');
    fprintf(fid,'\nSANS restraints\n\n');
    to_be_deleted = '';
    for ks = 1:length(restraints.SANS)
        if isfield(model,'selected')
            model = rmfield(model,'selected');
        end
        model.selected{1} = [snum0 1];
        model.selected{2} = [snum0 2];
        pdbfile = sprintf('t_%i',ks);
        to_be_deleted = 't*.*';
        wr_pdb_selected(pdbfile,'SANS');
        [chi2,~,~,result,fit] = fit_SANS_by_cryson(restraints.SANS(ks).data,pdbfile,restraints.SANS(ks).illres);
        if isempty(chi2) || isnan(chi2)
            SANS_chi = 1e6;
            if interactive
                fprintf(2,'Warning: SANS fitting failed\n');
                fprintf(fid,'Warning: SANS fitting of curve %i failed\n',ks);
                fprintf(2,'%s',result);
            end
        else
            restraints.SANS(ks).fit = fit;
            restraints.SANS(ks).chi2 = chi2;
            fprintf(fid,'SANS curve %i fitted with chi^2 of %6.3f\n',ks,chi2);
            SANS_chi = SANS_chi + chi2;
            score_SAS = score_SAS + chi2;
            n_SAS = n_SAS + 1;
        end
    end
    chi2 = SANS_chi/length(restraints.SANS);
    fprintf(fid,'Mean SANS curve chi^2: %6.3f\n',chi2);
    delete(to_be_deleted);
    figure(101); clf
    hold on;
    for ks = 1:length(restraints.SANS)
        fit = restraints.SANS(ks).fit;
        plot(fit(:,1),fit(:,2));
        plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
    end
    title(sprintf('SANS fit for model (chi^2 = %4.2f)',chi2));
end

SAXS_chi = 0;
if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS)
    add_msg_board('Fitting SAXS restraints');
    fprintf(fid,'\nSAXS restraints\n\n');
    to_be_deleted = '';
    for ks = 1:length(restraints.SAXS)
        if isfield(model,'selected')
            model = rmfield(model,'selected');
        end
        model.selected{1} = [snum0 1];
        model.selected{2} = [snum0 2];
        pdbfile = sprintf('t_%i',ks);
        to_be_deleted = 't*.*';
        wr_pdb_selected(pdbfile,'SAXS');
        SAXS_curve = load_SAXS_curve(restraints.SAXS(ks).data);
        options.sm = 10*max(SAXS_curve(:,1));
        options.lm = 50;
        options.fb = 18;
        [chi2,~,~,result,fit] = fit_SAXS_by_crysol(restraints.SAXS(ks).data,pdbfile,options);
        if isempty(chi2) || isnan(chi2)
            SAXS_chi = 1e6;
            if interactive
                fprintf(2,'Warning: SAXS fitting failed\n');
                fprintf(fid,'Warning: SAXS fitting of curve %i failed\n',ks);
                fprintf(2,'%s',result);
            end
        else
            restraints.SAXS(ks).fit = fit;
            restraints.SAXS(ks).chi2 = chi2;
            fprintf(fid,'SAXS curve %i fitted with chi^2 of %6.3f\n',ks,chi2);
            SAXS_chi = SAXS_chi + chi2;
            score_SAS = score_SAS + chi2;
            n_SAS = n_SAS + 1;
        end
    end
    chi2 = SAXS_chi/length(restraints.SAXS);
    fprintf(fid,'Mean SAXS curve chi^2: %6.3f\n',chi2);
    delete(to_be_deleted);
    figure(102); clf;
    hold on;
    for ks = 1:length(restraints.SAXS)
        fit = restraints.SAXS(ks).fit;
        plot(fit(:,1),fit(:,2));
        plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
    end
    title(sprintf('SAXS fit for model (chi^2 = %4.2f)',chi2));
end

score_SAS = score_SAS/n_SAS;

% test for superposition of RRMs
fprintf(fid,'\nTesting for RRM superposition with NMR models\n\n');

adr1 = '[PTB7](A)58-154';
adr2 = '(B)58-154';

[rmsd,~,msg] = backbone_overlap_peptide(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

adr1 = '[PTB7](C)182-283';
adr2 = '(B)182-283';

[rmsd,~,msg] = backbone_overlap_peptide(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

adr1 = '[PTB7](E)337-430';
adr2 = '(B)337-430';

[rmsd,~,msg] = backbone_overlap_peptide(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

adr1 = '[PTB7](E)454-531';
adr2 = '(B)454-531';

[rmsd,~,msg] = backbone_overlap_peptide(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

fprintf(fid,'\nTesting for RNA binding motif superposition\n\n');

adr1 = '[PTB7](B)327-329';
adr2 = '(A)327-329';

[rmsd,~,msg] = overlap_nucleic_acid(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

adr1 = '[PTB7](D)358-360';
adr2 = '(A)358-360';

[rmsd,~,msg] = overlap_nucleic_acid(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

adr1 = '[PTB7](F)342-344';
adr2 = '(A)342-344';

[rmsd,~,msg] = overlap_nucleic_acid(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

adr1 = '[PTB7](G)302-304';
adr2 = '(A)302-304';

[rmsd,~,msg] = overlap_nucleic_acid(adr1,adr2);
if msg.error
    fprintf(2,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
    fprintf(fid,'Sections %s and %s could not be superimposed (%s)\n',adr1,adr2,msg.text);
else
    fprintf(fid,'%s overlaps with %s rmsd %5.2f Å\n',adr2,adr1,rmsd);
end

if ~isempty(restraints.xlinks)
    [KK_pairs,KE_pairs] = find_crosslinkable_pairs([snum0 1]);
    
    [nKK,~] = size(KK_pairs);
    [nKE,~] = size(KE_pairs);
    
    fprintf(fid,'\nAnalysis of crosslink restraints\n\n');
    
    for k = 1:length(restraints.xlinks)
        adr1 = restraints.xlinks(k).adr1;
        adr2 = restraints.xlinks(k).adr2;
        cpoi = strfind(adr1,'(');
        switch adr1(cpoi+1)
            case {'A','C','E'}
                adr1(cpoi+1) = 'B';
            otherwise
                adr1(cpoi+1) = 'A';
        end
        switch adr2(cpoi+1)
            case {'A','C','E'}
                adr2(cpoi+1) = 'B';
            otherwise
                adr2(cpoi+1) = 'A';
        end
        ind1 = resolve_address(adr1);
        ind2 = resolve_address(adr2);
        for k1 = 1:nKK
            det1 = sum(abs(ind1 - KK_pairs(k1,1:4))) + sum(abs(ind2 - KK_pairs(k1,5:8)));
            det2 = sum(abs(ind2 - KK_pairs(k1,1:4))) + sum(abs(ind1 - KK_pairs(k1,5:8)));
            if det1 == 0 || det2 == 0
                fprintf(fid,'KK crosslink %s-%s has CA-CA distance %4.1f Å\n',adr1,adr2,KK_pairs(k1,9));
            end
        end
        for k1 = 1:nKE
            det1 = sum(abs(ind1 - KE_pairs(k1,1:4))) + sum(abs(ind2 - KE_pairs(k1,5:8)));
            det2 = sum(abs(ind2 - KE_pairs(k1,1:4))) + sum(abs(ind1 - KE_pairs(k1,5:8)));
            if det1 == 0 || det2 == 0
                fprintf(fid,'KE crosslink %s-%s has CA-CA distance %4.1f Å\n',adr1,adr2,KE_pairs(k1,9));
            end
        end
    end
    
    fprintf(fid,'\nList of potentially crosslinkable pairs (<= 20 Å CA-CA distance)\n\n');
    
    for k1 = 1:nKK
        adr1 = mk_address(KK_pairs(k1,1:4));
        adr2 = mk_address(KK_pairs(k1,5:8));
        if KK_pairs(k1,9) <= 20
            fprintf(fid,'Potential KK crosslink %s-%s has CA-CA distance %4.1f Å\n',adr1,adr2,KK_pairs(k1,9));
        end
    end
    
    for k1 = 1:nKE
        adr1 = mk_address(KE_pairs(k1,1:4));
        adr2 = mk_address(KE_pairs(k1,5:8));
        if KE_pairs(k1,9) <= 20
            fprintf(fid,'Potential KE crosslink %s-%s has CA-CA distance %4.1f Å\n',adr1,adr2,KE_pairs(k1,9));
        end
    end
end

fprintf(fid,'\nFinal scores\n\n');
fprintf(fid,'Mean DEER overlap deficiency: %5.3f\n',score_DEER);
fprintf(fid,'Mean SAS chi_squared        : %5.2f\n',score_SAS);

fclose(fid);

save(oname,'restraints');

% Transform to superposition frame and store

adr2 = strcat('%s%s',adr0,restraints.superimpose);
adr1 = '[PTB7](E)454-531';
[rmsd,transmat,msg] = backbone_overlap_peptide(adr1,adr2);
if msg.error
    add_msg_board(sprintf('ERROR: Superposition of model onto template failed (%s)',msg.text));
    set(hfig,'Pointer','arrow');
    return
else
    add_msg_board(sprintf('Superposition rmsd is %5.2f Ã',rmsd));
end

transform_structure(snum0,transmat);

model.selected = cell(1,2);
model.selected{1} = [snum0 1];
model.selected{2} = [snum0 2];

message = wr_pdb_selected(fname,'PTB1');

if message.error
    add_msg_board(strcat('ERROR (wr_pdb_selected): ',message.text));
else
    add_msg_board(sprintf('Structure file %s written',fname));
end


fid = fopen(list_file,'at');
if fid == -1
    add_msg_board('Solution score could not be recorded.');
    return
end

fprintf(fid,'%5.3f    %5.2f   %% %s\n',score_DEER,score_SAS,name);

fclose(fid);

set(hfig,'Pointer','arrow');
