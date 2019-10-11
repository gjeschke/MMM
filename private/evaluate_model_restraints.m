function [restraints,score_DEER,score_SAS] = evaluate_model_restraints(fname0,restraints,options)
% [restraints,score_DEER,score_SAS] = evaluate_model_restraints(fname0,restraints,options)
%
% Evaluates restraint data for one conformation (model) of a protein or
% complex, the evaluated data are added to structure restraints
%
% fname0        name of the PDB file containing the model
% restraints    RigiFlex or domain ensemble restraints, to be read with
%               rd_restraints_rigiflex or rd_restraints
% options       processing options, optional structure with fields
%               .depth  min/max modulation depth for DEER 
%               .dim    DEER background dimension
%               .kdec   min/max DEER decay time constant
%               .cutoff fractional cutoff of DEER primary trace 
%               .lm     max order of harmonics in SAXS curve computation
%               .fb     order of Fibonacci grid in SAXS curve computation
%               .xlmax  max cross-link distance
%
% score_DEER    total DEER overlap score
% score_SAS     sum of small-angle scattering chi^2 values (SANS & SAXS)
%
% G. Jeschke, 28.6.2019

global model

interactive = true;

if ~exist('options','var')
    % DEER time-domain fits
    options.depth = [0.005,0.65];
    options.dim = 3;
    options.kdec = [];
    options.cutoff = 0.9;
    % SAXS curve fits
    options.lm = 50;
    options.fb = 18;
    % cross-links
    options.xlmax = 20;
end

if ~isfield(options,'depth')
    options.depth = [0.005,0.65];
end

if ~isfield(options,'dim')
    options.dim = 3;
end

if ~isfield(options,'kdec')
    options.kdec = [];
end

if ~isfield(options,'cutoff')
    options.cutoff = 0.9;
end

if ~isfield(options,'lm')
    options.lm = 50;
end

if ~isfield(options,'fb')
    options.fb = 18;
end

if ~isfield(options,'xlmax')
    options.xlmax = 20;
end

options.err = true;

[~,snum0] = add_pdb(fname0);

score_SAS = 0;
score_DEER = 1;
n_DEER = 0;
n_SAS = 0;

fid = fopen(sprintf('%s_diagnosis.dat',fname0),'wt');
if fid == -1
    add_msg_board('Quality report file could not be written.');
    return
end

model.current_structure = snum0;

hfig=gcf;
set(hfig,'Pointer','watch');

% Evaluate DEER restraints in rigid bodies

fprintf(fid,'DEER restraints for model %s\n\n',fname0);
fprintf(fid,'>>> core and RNA\n\n');
for k = 1:length(restraints.DEER)
    restraints.DEER(k).adr1 = correct_address(restraints.DEER(k).adr1,restraints);
    restraints.DEER(k).adr2 = correct_address(restraints.DEER(k).adr2,restraints);
    [rax,distr] = mk_distance_distribution(restraints.DEER(k).adr1,restraints.DEER(k).adr2,restraints.DEER(k).label);
    restraints.DEER(k).rax = rax;
    restraints.DEER(k).distr = distr;
    if isfield(restraints.DEER(k),'file') && ~isempty(restraints.DEER(k).file) && ~isempty(rax)
        [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax,distr,strcat('deer\',restraints.DEER(k).file),options);
        restraints.DEER(k).texp = texp;
        restraints.DEER(k).vexp = vexp;
        restraints.DEER(k).deer = deer;
        restraints.DEER(k).bckg = bckg;
        restraints.DEER(k).param = param;
        fprintf(fid,'%s-%s rmsd: %6.4f, mod. depth: %6.3f, kdec: %6.4f\n',...
            restraints.DEER(k).adr1,restraints.DEER(k).adr2,...
            param.rmsd,param.depth,param.kdec);
    else
        restraints.DEER(k).texp = [];
        restraints.DEER(k).vexp = [];
        restraints.DEER(k).deer = [];
        restraints.DEER(k).bckg = [];
        restraints.DEER(k).param = [];
    end
    if isempty(rax)
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
    argr = (restraints.DEER(k).r-rax)/(sqrt(2)*restraints.DEER(k).sigr);
    distr_sim = exp(-argr.^2);
    distr_sim = distr_sim/sum(distr_sim);
    restraints.DEER(k).overlap = sum(min([distr_sim;distr]));
    if ~isnan(restraints.DEER(k).overlap)
        score_DEER = score_DEER*restraints.DEER(k).overlap;
        n_DEER = n_DEER + 1;
    end
    fprintf(fid,'Overlap: %5.3f\n',restraints.DEER(k).overlap);
    if isfield(restraints.DEER(k),'file') && ~isempty(restraints.DEER(k).file)
        dfname = strcat(restraints.DEER(k).file,'_distr.dat');
        dfname = strcat('deer_analysis\',dfname);
        Pdata = load(dfname);
        rexp = Pdata(:,1).';
        distr_exp = Pdata(:,2).';
        distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
        distr_exp = distr_exp/sum(distr_exp);
        restraints.DEER(k).distr_exp = distr_exp;
    else
        restraints.DEER(k).distr_exp = [];
    end
end

% Evaluate DEER restraints in flexible linkers

if isfield(restraints,'pflex') && ~isempty(restraints.pflex) && isfield(restraints.pflex,'DEER')
    fprintf(fid,'\n>>> flexible peptide linkers\n\n');
    for kl = 1:length(restraints.pflex)
        for k = 1:length(restraints.pflex(kl).DEER)
            restraints.pflex(kl).DEER(k).adr1 = correct_address(restraints.pflex(kl).DEER(k).adr1,restraints);
            restraints.pflex(kl).DEER(k).adr2 = correct_address(restraints.pflex(kl).DEER(k).adr2,restraints);
            label = [restraints.pflex(kl).DEER(k).label1 '|' restraints.pflex(kl).DEER(k).label2];
            [rax,distr] = mk_distance_distribution(restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2,label);
            restraints.pflex(kl).DEER(k).rax = rax;
            restraints.pflex(kl).DEER(k).distr = distr;
            if isfield(restraints.pflex(kl).DEER(k),'file') && ~isempty(restraints.pflex(kl).DEER(k).file) && ~isempty(rax)
                [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax,distr,strcat('deer\',restraints.pflex(kl).DEER(k).file),options);
                restraints.pflex(kl).DEER(k).texp = texp;
                restraints.pflex(kl).DEER(k).vexp = vexp;
                restraints.pflex(kl).DEER(k).deer = deer;
                restraints.pflex(kl).DEER(k).bckg = bckg;
                restraints.pflex(kl).DEER(k).param = param;
                fprintf(fid,'%s-%s rmsd: %6.4f, mod. depth: %6.3f, kdec: %6.4f\n',...
                    restraints.pflex(kl).DEER(k).adr1,restraints.pflex(kl).DEER(k).adr2,...
                    param.rmsd,param.depth,param.kdec);
            else
                restraints.pflex(kl).DEER(k).texp = [];
                restraints.pflex(kl).DEER(k).vexp = [];
                restraints.pflex(kl).DEER(k).deer = [];
                restraints.pflex(kl).DEER(k).bckg = [];
                restraints.pflex(kl).DEER(k).param = [];
            end
            if isempty(rax)
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
            argr = (restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*restraints.pflex(kl).DEER(k).sigr);
            distr_sim = exp(-argr.^2);
            distr_sim = distr_sim/sum(distr_sim);
            restraints.pflex(kl).DEER(k).overlap = sum(min([distr_sim;distr]));
            if ~isnan(restraints.DEER(k).overlap)
                score_DEER = score_DEER*restraints.DEER(k).overlap;
                n_DEER = n_DEER + 1;
            end
            fprintf(fid,'Overlap: %5.3f\n',restraints.pflex(kl).DEER(k).overlap);
            if isfield(restraints.pflex(kl).DEER(k),'file') && ~isempty(restraints.pflex(kl).DEER(k).file)
                dfname = strcat(restraints.pflex(kl).DEER(k).file,'_distr.dat');
                dfname = strcat('deer_analysis\',dfname);
                Pdata = load(dfname);
                rexp = Pdata(:,1).';
                distr_exp = Pdata(:,2).';
                distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
                distr_exp = distr_exp/sum(distr_exp);
                restraints.pflex(kl).DEER(k).distr_exp = distr_exp;
            else
                restraints.pflex(kl).DEER(k).distr_exp = [];
            end
        end
    end
end
score_DEER = 1 - score_DEER^(1/n_DEER);

% Evaluate SANS restraints

SANS_chi = 0;
if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
    add_msg_board('Fitting SANS restraints');
    fprintf(fid,'\nSANS restraints\n\n');
    to_be_deleted = '';
    for ks = 1:length(restraints.SANS)
        if isfield(model,'selected')
            model = rmfield(model,'selected');
        end
        model.selected{1} = snum0;
        pdbfile = sprintf('t_%i',ks);
        to_be_deleted = 't*.*';
        wr_pdb_selected(pdbfile,'SANS');
        options.D2O = restraints.SANS(ks).D2O;
        [chi2,~,~,result,fit] = fit_SANS_by_cryson(restraints.SANS(ks).data,pdbfile,restraints.SANS(ks).illres,options);
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
        model.selected{1} = snum0;
        pdbfile = sprintf('t_%i',ks);
        to_be_deleted = 't*.*';
        wr_pdb_selected(pdbfile,'SAXS');
        SAXS_curve = load_SAXS_curve(restraints.SAXS(ks).data);
        if ~isempty(restraints.SAXS(ks).sm)
            options.sm = restraints.SAXS(ks).sm;
        elseif max(SAXS_curve(:,1)) < 1
            options.sm = 10*max(SAXS_curve(:,1));
        else
            options.sm = max(SAXS_curve(:,1));
        end
        smin = min(SAXS_curve(:,1));
        if options.sm > 1
            smin = smin/10;
        end
        options.smin = smin;
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
end

score_SAS = score_SAS/n_SAS;

if isfield(restraints,'xlinks') && ~isempty(restraints.xlinks)
    [KK_pairs,KE_pairs] = find_crosslinkable_pairs([snum0 1]);
    
    [nKK,~] = size(KK_pairs);
    [nKE,~] = size(KE_pairs);
    
    fprintf(fid,'\nAnalysis of crosslink restraints\n\n');
    
    for k = 1:length(restraints.xlinks)
        restraints.xlinks(k).adr1 = correct_address(restraints.xlinks(k).adr1,restraints);
        restraints.xlinks(k).adr2 = correct_address(restraints.xlinks(k).adr2,restraints);
        adr1 = restraints.xlinks(k).adr1;
        adr2 = restraints.xlinks(k).adr2;
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
    
    fprintf(fid,'\nList of potentially crosslinkable pairs (<= %4.1f Å CA-CA distance)\n\n',options.xlmax);
    
    for k1 = 1:nKK
        adr1 = mk_address(KK_pairs(k1,1:4));
        adr2 = mk_address(KK_pairs(k1,5:8));
        if KK_pairs(k1,9) <= options.xlmax
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

fname = sprintf('%s_restraints.mat',fname0);
save(fname,'restraints','score_DEER','score_SAS');

set(hfig,'Pointer','arrow');

function newadr = correct_address(oldadr,restraints)

Aminus1 = 64; % ASCII code of the last character before 'A'

if isfield(restraints,'substitute') && ~isempty(restraints.substitute)
    pa = strfind(oldadr,'(');
    pe = strfind(oldadr,')');
    if pa > 1
        pre = oldadr(1:pa-1);
    else
        pre = '';
    end
    if pe < length(oldadr)
        suff = oldadr(pe+1:end);
    else
        suff = '';
    end
    if ~isempty(pa) && ~isempty(pe)
       chain = oldadr(pa:pe);
       for kc = 1:length(restraints.substitute)
           tbs = false;
           for ks = 1:length(restraints.substitute(kc).chains)
               if strcmpi(chain,restraints.substitute(kc).chains{ks})
                   tbs = true;
                   break
               end
           end
           if tbs
               chain = sprintf('(%c)',char(Aminus1+kc));
               break
           end
       end
       newadr = [pre chain suff];
    else
        newadr = oldadr;
    end
else
    newadr = oldadr;
end

