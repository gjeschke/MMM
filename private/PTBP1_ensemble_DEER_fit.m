function PTBP1_ensemble_DEER_fit(ensemble_size)

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

% Make DEER information
fprintf(1,'Minimum normalized score     : %5.2f at model %i\n',sorted_score(1),score_poi(1));
[overlaps,fits] = get_DEER(flist{score_poi(1)});
if isempty(overlaps)
    fprintf(2,'ERROR: Not all restraints could be evaluated for first model %s. Aborting\n',flist{score_poi(1)});
    return
end
n_DEER = length(overlaps);
score_DEER = 1 - prod(overlaps)^(1/n_DEER);
fid = fopen('PTBP1_DEER_scores.dat','at');
if fid == -1
    add_msg_board('Solution score could not be recorded.');
    return
end
fprintf(fid,'%5.3f    %5.2f   %% %s\n',score_DEER,sorted_score(1),flist{score_poi(1)});
fclose(fid);

all_fits = cell(1,length(fits));
for k = 1:length(fits)
    dat0 = fits{k};
    [m,~] = size(dat0);
    data = zeros(m,2+ensemble_size);
    data(:,1:3) = dat0(:,1:3);
    all_fits{k} = data;
end

for km = 2:ensemble_size
    fprintf(1,'Normalized score         : %5.2f at model %i\n',sorted_score(km),score_poi(km));
    [overlaps,fits] = get_DEER(flist{score_poi(km)});
    if isempty(overlaps)
        fprintf(2,'Warning: For model %s not all DEER restraints could be evaluated. Model skipped.\n',flist{score_poi(km)});
        continue
    else
        score_DEER = 1 - prod(overlaps)^(1/n_DEER);
        fid = fopen('PTBP1_DEER_scores.dat','at');
        if fid == -1
            add_msg_board('Solution score could not be recorded.');
            return
        end      
        fprintf(fid,'%5.3f    %5.2f   %% %s\n',score_DEER,sorted_score(km),flist{score_poi(km)});       
        fclose(fid);
    end
    for k = 1:length(fits)
        dat0 = fits{k};
        data = all_fits{k};
        data(:,2+km) = dat0(:,3);
        all_fits{k} = data;
    end
end

fanonym=@(v_opt)fit_multi_DEER(v_opt,all_fits);

l=zeros(1,ensemble_size);
u=ones(1,ensemble_size);
v0 = ones(1,ensemble_size)/ensemble_size;

fom0 = fit_multi_DEER(v0,all_fits);

options = optimoptions('patternsearch','MaxFunctionEvaluations',10000*length(v0),'MaxIterations',500*length(v0));

tic,
[v,fom1,exitflag,output] = patternsearch(fanonym,v0,[],[],[],[],l,u,[],options);
toc,

[fom,fits] = sim_multi_DEER(v,all_fits);

coeff = v(1:ensemble_size);
coeff = coeff/sum(coeff);

figure(1); clf;
plot(coeff,'k.');

for k = 1:n_DEER
    figure(k+1); clf;
    plot(fits(k).rax,fits(k).restraint,'g');
    hold on;
    plot(fits(k).rax,fits(k).distr,'r');
    title(sprintf('DEER fit with overlap %5.3f',sum(min([fits(k).restraint';fits(k).distr']))));
end
fprintf(1,'Total DEER overlap improved from %6.3f to %6.3f\n',fom0,fom);

save PTB1_DEER_coeff coeff fom1 exitflag output v all_fits

function [overlaps,fits] = get_DEER(fname)

fprintf(1,'Model % s:\n',fname);

rname = sprintf('%s_diagnosis_data.mat',fname);

s = load(rname,'restraints');

% make DEER distance distributions

core = length(s.restraints.DEER);
rax = s.restraints.DEER(1).rax;
flex = 0;

for kl = 1:length(s.restraints.pflex)
    flex = flex + length(s.restraints.pflex(kl).DEER);
end

overlaps = zeros(1,core+flex);
fits = cell(1,core+flex);

n_DEER = 0;
score_DEER = 1;
DEER_valid = true;

for k = 1:length(s.restraints.DEER)
    if s.restraints.DEER(k).r ~=0 &&  s.restraints.DEER(k).sigr ~=0 % this is a restraint to be tested
        distr = s.restraints.DEER(k).distr;
        n_DEER = n_DEER + 1;
        if ~isempty(distr)
            distr = distr/sum(distr);
            data = zeros(length(rax),3);
            data(:,1) = rax.';
            argr = (s.restraints.DEER(k).r-rax)/(sqrt(2)*s.restraints.DEER(k).sigr);
            distr_sim = exp(-argr.^2);
            distr_sim = distr_sim/sum(distr_sim);
            data(:,2) = distr_sim.';
            data(:,3) = distr.';
            fits{n_DEER} = data;
            overlaps(n_DEER) = sum(min([distr_sim;distr]));
            score_DEER = score_DEER * overlaps(n_DEER);
        else
            DEER_valid = false;
            break
        end
    end
end

if ~DEER_valid
    overlaps = [];
    fits = {};
    return
end

pflex = core;
for kl = 1:length(s.restraints.pflex)
    for k = 1:length(s.restraints.pflex(kl).DEER)
        pflex = pflex + 1;
        if s.restraints.pflex(kl).DEER(k).r ~=0 &&  s.restraints.pflex(kl).DEER(k).sigr ~=0 % valid restraint
            n_DEER = n_DEER + 1;
            distr = s.restraints.pflex(kl).DEER(k).distr;
            if ~isempty(distr)
                distr = distr/sum(distr);
                data = zeros(length(rax),3);
                data(:,1) = rax.';
                argr = (s.restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*s.restraints.pflex(kl).DEER(k).sigr);
                distr_sim = exp(-argr.^2);
                distr_sim = distr_sim/sum(distr_sim);
                data(:,2) = distr_sim.';
                data(:,3) = distr.';
                fits{n_DEER} = data;
                overlaps(n_DEER) = sum(min([distr_sim;distr]));
                score_DEER = score_DEER * overlaps(n_DEER);
            else
                DEER_valid = false;
                break
            end
        end
    end
end

score_DEER = 1 - score_DEER^(1/n_DEER);
fprintf(1,'Total DEER score : %5.3f\n',score_DEER);

if ~DEER_valid
    overlaps = [];
    fits = {};
    return
end

overlaps = real(overlaps(1:n_DEER));
fits = fits(1:n_DEER);
