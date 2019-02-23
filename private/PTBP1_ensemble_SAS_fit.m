function PTBP1_ensemble_SAS_fit(ensemble_size)

if ~exist('ensemble_size','var')
    ensemble_size = 20;
end

scores0 = load('PTBP1_solution_scores.dat');
[m,n] = size(scores0);
scores = scores0;
if ensemble_size > m
    ensemble_size = m;
end

% normalize scores
for k = 1:n
    scores(:,k) = scores(:,k)/min(scores(:,k));
end
total_score = sum(scores,2);
[sorted_score,score_poi] = sort(total_score);


fid=fopen('PTBP1_solution_scores.dat');
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
%         if ~contains(fname,'.pdb')
%             fname=strcat(fname,'.pdb'); 
%         end
        flist{nf} = fname;
    end
end
fclose(fid);
flist = flist(1:nf);

% Make SAS information
fprintf(1,'Minimum normalized score     : %5.2f at model %i\n',sorted_score(1),score_poi(1));
[chi2_vec,fits] = get_SAS(flist{score_poi(1)});
all_fits = cell(1,length(fits));
all_chi2 = zeros(ensemble_size,length(fits));
all_msq = zeros(ensemble_size,length(fits));
all_chi2(1,:) = chi2_vec; 
np = zeros(1,length(fits));
min_fom = 1e6;
fom = 0;
best_model = 0;
for k = 1:length(fits)
    dat0 = fits{k};
    [m,~] = size(dat0);
    np(k) = m;
    data = zeros(m,3+ensemble_size);
    data(:,1:2) = dat0(:,1:2);
    data(:,3) = dat0(:,4);
    data(:,4) = dat0(:,3);
    my_chi2 = sum(((dat0(:,2)-dat0(:,3))./dat0(:,4)).^2);
    all_msq(1,k) = my_chi2/(m-1);
    fom = fom + my_chi2/(m-1);
    all_fits{k} = data;
end
if min_fom > fom
    min_fom = fom;
    best_model = 1;
end

for km = 2:ensemble_size
    fprintf(1,'Normalized score         : %5.2f at model %i\n',sorted_score(km),score_poi(km));
    [chi2_vec,fits] = get_SAS(flist{score_poi(km)});
    all_chi2(km,:) = chi2_vec;
    fom = 0;
    for k = 1:length(fits)
        dat0 = fits{k};
        data = all_fits{k};
        data(:,3+km) = dat0(:,3);
        my_chi2 = sum(((dat0(:,2)-dat0(:,3))./dat0(:,4)).^2);
        all_msq(km,k) = my_chi2/(np(k)-1);
        fom = fom + my_chi2/(np(k)-1);
        all_fits{k} = data;
    end
    if min_fom > fom
        min_fom = fom;
        best_model = km;
    end
end

fprintf(1,'Best total SAS score was %5.2f at model %i (%s)\n',min_fom,best_model,flist{score_poi(best_model)});

fanonym=@(v_opt)fit_multi_SAS(v_opt,all_fits);

l=[zeros(1,ensemble_size) -0.005 -0.005 -0.01];
u=[ones(1,ensemble_size) 0.005 0.005 0.01];
v0 = [1 zeros(1,ensemble_size-1) 0 0 0];

fom0 = fit_multi_SAS(v0,all_fits);

tic,
v = patternsearch(fanonym,v0,[],[],[],[],l,u);
toc,

[fom,fits] = sim_multi_SAS(v,all_fits);

coeff = v(1:ensemble_size);
coeff = coeff/sum(coeff);

figure(1); clf;
plot(coeff,'k.');

figure(2); clf;
plot(fits(1).s,fits(1).curve);
hold on;
plot(fits(1).s,fits(1).sim);
plot(fits(2).s,fits(2).curve);
plot(fits(2).s,fits(2).sim);
title(sprintf('SANS curve fits with chi2 = %5.2f and %5.2f',fits(1).chi2,fits(2).chi2));

figure(3); clf;
plot(fits(3).s,fits(3).curve);
hold on;
plot(fits(3).s,fits(3).sim);
title(sprintf('SAXS curve fit with chi2 = %5.2f',fits(3).chi2));

fprintf(1,'Offset of SANS(1.2 m) data: %8.4f\n',v(ensemble_size+1));
fprintf(1,'Offset of SANS(4 m) data  : %8.4f\n',v(ensemble_size+2));
fprintf(1,'Offset of SAXS data       : %8.4f\n',v(ensemble_size+3));
fprintf(1,'Total SAS chi2 improved from %6.2f to %6.2f\n',fom0,fom);

save PTB1_SAS_fits v coeff fom fits all_fits

function [chi2_vec,fits] = get_SAS(fname)

global model

options.err = true;
options.lm = 50;
options.fb = 18;
options.D2O = 0.66;

fprintf(1,'Model % s:\n',fname);

[~,snum0] = add_pdb(fname);

chi2_vec = zeros(1,3);
fits = cell(1,3);

rname = sprintf('%s_diagnosis_data.mat',fname);

s = load(rname,'restraints');

n_SANS = 0;
if isfield(s.restraints,'SANS') && ~isempty(s.restraints.SANS)
    n_SANS = length(s.restraints.SANS);
    for ks = 1:length(s.restraints.SANS)
        if isfield(model,'selected')
            model = rmfield(model,'selected');
        end
        model.selected{1} = [snum0 1];
        model.selected{2} = [snum0 2];
        pdbfile = sprintf('t_%i',ks);
        to_be_deleted = 't*.*';
        wr_pdb_selected(pdbfile,'SANS');
        [chi2,~,~,result,fit] = fit_SANS_by_cryson(s.restraints.SANS(ks).data,pdbfile,s.restraints.SANS(ks).illres,options);
        if isempty(chi2) || isnan(chi2)
            fprintf(2,'Warning: SANS fitting failed\n');
            fprintf(fid,'Warning: SANS fitting of curve %i failed\n',ks);
            fprintf(2,'%s',result);
        else
            fits{ks} = fit;
            chi2_vec(ks) = chi2;
        end
       fprintf(1,'SANS fit for model (chi^2 = %4.2f)\n',chi2_vec(ks));
       delete(to_be_deleted);
    end
end



if isfield(s.restraints,'SAXS') && ~isempty(s.restraints.SAXS)
    for ks = 1:length(s.restraints.SAXS)
        if isfield(model,'selected')
            model = rmfield(model,'selected');
        end
        model.selected{1} = [snum0 1];
        model.selected{2} = [snum0 2];
        pdbfile = sprintf('t_%i',ks);
        to_be_deleted = 't*.*';
        wr_pdb_selected(pdbfile,'SAXS');
        SAXS_curve = load_SAXS_curve(s.restraints.SAXS(ks).data);
        sm = max(SAXS_curve(:,1));
        smin = min(SAXS_curve(:,1));
        options.smin = smin;
        options.sm = 10*sm;
        [chi2,~,~,result,fit] = fit_SAXS_by_crysol(s.restraints.SAXS(ks).data,pdbfile,options);
        if isempty(chi2) || isnan(chi2)
            fprintf(2,'Warning: SAXS fitting failed\n');
            fprintf(2,'%s',result);
        else
            fits{n_SANS+ks} = fit;
            chi2_vec(n_SANS+ks) = chi2;
        end
        fprintf(1,'SAXS fit chi^2: %6.3f\n',chi2_vec(n_SANS+ks));
        delete(to_be_deleted);
    end
end

