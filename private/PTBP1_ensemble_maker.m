function PTBP1_ensemble_maker(fname)

global model

SAS_fit_file = sprintf('%s_SAS_fit.mat',fname);
DEER_fit_file = sprintf('%s_DEER_fit.mat',fname);
definition_file = sprintf('%s_ensemble_definition.mat',fname);
pdb_file = sprintf('%s_ensemble.pdb',fname);
script_file = sprintf('%s_ensemble.mmm',fname);
model_list = 'best_models.dat';

score_file = sprintf('%s_scores.dat',fname);
scores0 = load(score_file);
scores = scores0;

[~,n] = size(scores);
% normalize scores
for k = 1:n
    scores(:,k) = scores(:,k)/min(scores(:,k));
end
total_score = sum(scores,2);
[~,score_poi] = sort(total_score);
sorted_scores = scores(score_poi,:);
sorted_scores0 = scores0(score_poi,:);

fid=fopen(score_file);
if fid==-1
    add_msg_board('ERROR: Ensemble file list does not exist');
    return;
end
flist = cell(1,200);
nf = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break,
    end
    if ~isempty(tline)
        nf = nf + 1;
        [~,remain] = strtok(tline,'%');
        fname = strtok(remain(2:end));
        flist{nf} = fname;
    end
end
fclose(fid);
flist = flist(1:nf);

SAS_fits = load(SAS_fit_file);
all_fits_SAS = SAS_fits.all_fits;

% load PTB1_DEER_coeff
DEER_fits = load(DEER_fit_file);
all_fits_DEER = DEER_fits.all_fits;

def = load(definition_file);
v = def.v;

figure(100); clf;
plot(v(1:end-3),'k.');
% remove insignificant contributions
coeff = v(1:end-3);
coeff = coeff/max(coeff);
coeff(coeff<1e-2) = 0;
sc_coeff = sum(coeff);
v(1:end-3) = coeff;
v = v/sc_coeff;

hold on
plot(v(1:end-3),'ro');

fprintf(1,'Fit with %6.4f\n',def.fom);

[fom_SAS,fits] = sim_multi_SAS(v,all_fits_SAS);

fprintf(1,'Total SAS chi2 %6.4f\n',fom_SAS);

figure(102); clf;
plot(fits(1).s,fits(1).curve);
hold on;
plot(fits(1).s,fits(1).sim);
plot(fits(2).s,fits(2).curve);
plot(fits(2).s,fits(2).sim);
title(sprintf('SANS curve fits with chi2 = %5.2f and %5.2f',fits(1).chi2,fits(2).chi2));

figure(103); clf;
plot(fits(3).s,fits(3).curve);
hold on;
plot(fits(3).s,fits(3).sim);
title(sprintf('SAXS curve fit with chi2 = %5.2f',fits(3).chi2));

[fom_DEER,fits_DEER] = sim_multi_DEER(v(1:end-3),all_fits_DEER);

fprintf(1,'DEER overlap deficiency %6.4f\n',fom_DEER);

n_DEER = length(fits_DEER);
for k = 1:n_DEER
    figure(k+3); clf;
    plot(fits_DEER(k).rax,fits_DEER(k).restraint,'g');
    hold on;
    plot(fits_DEER(k).rax,fits_DEER(k).distr,'r');
    title(sprintf('DEER fit with overlap %5.3f',sum(min([fits_DEER(k).restraint';fits_DEER(k).distr']))));
end

flists = flist(score_poi);
coeff = v(1:end-3);
coeff1 = coeff(coeff>0);
flist1 = flists(coeff>0);
s_scores1 = sorted_scores0(coeff>0,:);

[coeff2,poi] = sort(coeff1,'descend');
flist2 = flist1(poi);
s_scores2 = s_scores1(poi,:);

% Make ensemble PDB file
fprintf(1,'Population of model  1: %5.2f\n',coeff2(1));
[msg,snum0] = add_pdb(flist2{1});
if msg.error
    add_msg_board(sprintf('ERROR: First PDB file %s could not be read (%s).',flist{score_poi(1)},msg.text));
    return
end
for kf = 2:length(flist2)
    copy_structure(snum0,'+mod',[],kf,snum0); % make a copy of the first model
    fprintf(1,'Population of model %2i: %5.2f\n',kf,coeff2(kf));
    [msg,snumc] = add_pdb(flist2{kf});
    if msg.error
        add_msg_board(sprintf('ERROR: PDB file #%i (%s) could not be read (%s).',kf,flist{score_poi(kf)},msg.text));
        return
    end
    for kc = 1:2
        model.structures{snum0}(kc).atoms{kf} = model.structures{snumc}(kc).atoms{1};
        model.structures{snum0}(kc).residues{kf} =model.structures{snumc}(kc).residues{1};
        model.structures{snum0}(kc).xyz{kf} = model.structures{snumc}(kc).xyz{1};
        model.structures{snum0}(kc).Bfactor{kf} =  model.structures{snumc}(kc).Bfactor{1};
        model.structures{snum0}(kc).Btensor{kf} = model.structures{snumc}(kc).Btensor{1};
    end
end
model.selected = cell(1,1);
model.selected{1} = snum0;
message = wr_pdb_selected(pdb_file,'PTB1');

fid = fopen(model_list,'wt');
if fid == -1
    add_msg_board('Ensemble score file could not be recorded.');
    return
end

for km = 1:length(flist2)
    fprintf(fid,'%5.3f    %5.2f   %% %s\n',s_scores2(km,:),flist2{km});
end
fclose(fid);

if message.error
    add_msg_board(strcat('ERROR (wr_pdb_selected): ',message.text));
else
    add_msg_board('Structure file PTBP1_EMCV_IRES_ensemble_101.pdb written');
end

fid = fopen(script_file,'wt');
fprintf(fid,'show [PTB1] ribbon\n');
fprintf(fid,'color [PTB1](A){:} grey\n');
fprintf(fid,'color [PTB1](B){:}58-154 crimson\n');
fprintf(fid,'color [PTB1](B){:}155-181 goldenrod\n');
fprintf(fid,'color [PTB1](B){:}182-283 darkgreen\n');
fprintf(fid,'color [PTB1](B){:}284-336 slateblue\n');
fprintf(fid,'color [PTB1](B){:}337-531 darkblue\n');

for km = 1:length(coeff2)
    alpha = coeff2(km)/max(coeff2);
    fprintf(fid,'transparency [PTB1](:){%i} %5.3f\n',km,alpha);
end
fclose(fid);

rname = sprintf('%s_diagnosis_data.mat',flist2{1});

s = load(rname,'restraints');

% make DEER distance distributions

fit_options.depth = [0.005,0.65];
fit_options.dim = 3;
fit_options.kdec = [];
fit_options.cutoff = 0.85;

core = length(s.restraints.DEER);
rax = s.restraints.DEER(1).rax;
flex = 0;

for kl = 1:length(s.restraints.pflex)
    flex = flex + length(s.restraints.pflex(kl).DEER);
end

all_distr = zeros(core+flex,length(rax));
all_distr_sim = zeros(core+flex,length(rax));
all_distr_exp = zeros(core+flex,length(rax));
all_fnames = cell(1,core+flex);
all_adr1 = cell(1,core+flex);
all_adr2 = cell(1,core+flex);
all_flags = zeros(1,core+flex);

for k = 1:length(s.restraints.DEER)
    distr = s.restraints.DEER(k).distr;
    if ~isempty(distr)
        figure(k); clf; hold on;
        argr = (s.restraints.DEER(k).r-rax)/(sqrt(2)*s.restraints.DEER(k).sigr);
        distr_sim = exp(-argr.^2);
        distr_sim = distr_sim/sum(distr_sim);
        all_distr_sim(k,:) = distr_sim;
        plot(rax,distr_sim,'Color',[0,0.6,0]);
        title(sprintf('%s-%s',s.restraints.DEER(k).adr1,s.restraints.DEER(k).adr2));
        plot(rax,coeff2(1)*distr,'Color',[0.2,0.2,1]);
        all_distr(k,:) = all_distr(k,:) + coeff2(1)*distr;
        all_fnames{k} = s.restraints.DEER(k).file;
        if s.restraints.DEER(k).r ~=0 &&  s.restraints.DEER(k).sigr ~=0
            all_flags(k) = 1;
        end
        all_adr1{k} = s.restraints.DEER(k).adr1;
        all_adr2{k} = s.restraints.DEER(k).adr2;
    end
    if isfield(s.restraints.DEER(k),'file') && ~isempty(s.restraints.DEER(k).file)
        dfname = strcat(s.restraints.DEER(k).file,'_distr.dat');
        dfname = strcat('deer_analysis\',dfname);
        Pdata = load(dfname);
        rexp = Pdata(:,1).';
        distr_exp = Pdata(:,2).';
        [~,rpoi] = min(abs(rexp-10));
        rexp = rexp(1:rpoi)*10;
        distr_exp = distr_exp(1:rpoi);
        distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
        distr_exp = distr_exp/sum(distr_exp);
        all_distr_exp(k,:) = distr_exp;
    end
end

pflex = core;
for kl = 1:length(s.restraints.pflex)
    for k = 1:length(s.restraints.pflex(kl).DEER)
        pflex = pflex + 1;
        distr = s.restraints.pflex(kl).DEER(k).distr;
        if ~isempty(distr)
            figure(pflex); clf; hold on;
            title(sprintf('%s-%s',s.restraints.pflex(kl).DEER(k).adr1,s.restraints.pflex(kl).DEER(k).adr2));
            argr = (s.restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*s.restraints.pflex(kl).DEER(k).sigr);
            distr_sim = exp(-argr.^2);
            distr_sim = distr_sim/sum(distr_sim);
            all_distr_sim(pflex,:) = distr_sim;
            plot(rax,distr_sim,'Color',[0,0.6,0]);
            plot(rax,coeff2(1)*distr,'Color',[0.2,0.2,1]);
            all_distr(pflex,:) = all_distr(pflex,:) + coeff2(1)*distr;
            all_fnames{pflex} = s.restraints.pflex(kl).DEER(k).file;
            if s.restraints.pflex(kl).DEER(k).r ~=0 &&  s.restraints.pflex(kl).DEER(k).sigr ~=0
                all_flags(pflex) = 1;
            end
            all_adr1{pflex} = s.restraints.pflex(kl).DEER(k).adr1;
            all_adr2{pflex} = s.restraints.pflex(kl).DEER(k).adr2;
        end
        if isfield(s.restraints.pflex(kl).DEER(k),'file') && ~isempty(s.restraints.pflex(kl).DEER(k).file)
            dfname = strcat(s.restraints.pflex(kl).DEER(k).file,'_distr.dat');
            dfname = strcat('deer_analysis\',dfname);
            Pdata = load(dfname);
            rexp = Pdata(:,1).';
            distr_exp = Pdata(:,2).';
            [~,rpoi] = min(abs(rexp-10));
            rexp = rexp(1:rpoi)*10;
            distr_exp = distr_exp(1:rpoi);
            distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
            distr_exp = distr_exp/sum(distr_exp);
            all_distr_exp(pflex,:) = distr_exp;
        end
    end
end

for kf = 2:length(flist2)
    rname = sprintf('%s_diagnosis_data.mat',flist2{kf});
    
    s = load(rname,'restraints');
    
    for k = 1:length(s.restraints.DEER)
        distr = s.restraints.DEER(k).distr;
        if ~isempty(distr)
            figure(k); hold on;
            title(sprintf('%s-%s',s.restraints.DEER(k).adr1,s.restraints.DEER(k).adr2));
            plot(rax,coeff2(kf)*distr,'Color',[0.2,0.2,1]);
            all_distr(k,:) = all_distr(k,:) + coeff2(kf)*distr;
        end
    end
    pflex = core;
    for kl = 1:length(s.restraints.pflex)
        for k = 1:length(s.restraints.pflex(kl).DEER)
            pflex = pflex + 1;
            distr = s.restraints.pflex(kl).DEER(k).distr;
            if ~isempty(distr)
                figure(pflex); hold on;
                title(sprintf('%s-%s',s.restraints.pflex(kl).DEER(k).adr1,s.restraints.pflex(kl).DEER(k).adr2));
                plot(rax,coeff2(kf)*distr,'Color',[0.2,0.2,1]);
                all_distr(pflex,:) = all_distr(pflex,:) + coeff2(kf)*distr;
            end
        end
    end
end

score_DEER = 1;
n_DEER = 0;

for k = 1:core+flex
    all_distr(k,:) = all_distr(k,:)/sum(all_distr(k,:));
    figure(k); hold on;
    plot(rax,all_distr(k,:),'Color',[0.6,0,0]);
    if sum(all_distr_exp(k,:)) > 10*eps
        plot(rax,all_distr_exp(k,:),'k');
    end
    overlap = sum(min([all_distr_sim(k,:);all_distr(k,:)]));
    if ~isnan(overlap) && overlap > 10*eps
        score_DEER = score_DEER*overlap;
        n_DEER = n_DEER + 1;
    end
    if all_flags(k)
        ftype = 'fit';
    else
        ftype = 'control';
    end
    title(sprintf('%s: %s-%s. overlap: %5.3f',ftype,all_adr1{k},all_adr2{k},overlap));
    [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax/10,all_distr(k,:),strcat('deer\',all_fnames{k}),fit_options);
    figure(1000+k); clf;
    plot(texp,vexp,'k');
    hold on;
    plot(fit_options.cutoff*[max(texp),max(texp)],[min(vexp),max(vexp)],'b');
    plot(texp,deer,'Color',[0.75,0,0]);
    plot(texp,bckg,'Color',[0,0.6,0]);
    title(sprintf('%s: %s-%s. rmsd: %6.4f',ftype,all_adr1{k},all_adr2{k},param.rmsd));
end

score_DEER = 1 - score_DEER^(1/n_DEER);
fprintf(1,'Total DEER score : %5.3f\n',score_DEER);