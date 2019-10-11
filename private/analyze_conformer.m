function [overlaps,fits_DDR,score_DDR,chi2_vec,fits_SAS] = analyze_conformer(conformer)

overlaps = [];
fits_DDR = {};
score_DDR = [];
chi2_vec = [];
fits_SAS = {};

fname = sprintf('%s_restraints.mat',conformer);
if exist(fname,'file')
    s = load(fname);
    model_restraints = s.restraints;
else
    fprintf(2,'ERROR: No diagnostics for conformer %s\n',conformer);
    return
end
[overlaps,fits_DDR,score_DDR] = get_DEER(model_restraints);
[chi2_vec,fits_SAS] = get_SAS(model_restraints);


function [overlaps,fits,score_DEER] = get_DEER(restraints)

% make DEER distance distributions

core = length(restraints.DEER);
rax = restraints.DEER(1).rax;
flex = 0;

for kl = 1:length(restraints.pflex)
    flex = flex + length(restraints.pflex(kl).DEER);
end

overlaps = zeros(1,core+flex);
fits = cell(1,core+flex);

n_DEER = 0;
score_DEER = 1;
DEER_valid = true;

for k = 1:length(restraints.DEER)
    if restraints.DEER(k).r ~=0 &&  restraints.DEER(k).sigr ~=0 % this is a restraint to be tested
        distr = restraints.DEER(k).distr;
        n_DEER = n_DEER + 1;
        if ~isempty(distr)
            distr = distr/sum(distr);
            data = zeros(length(rax),3);
            data(:,1) = rax.';
            argr = (restraints.DEER(k).r-rax)/(sqrt(2)*restraints.DEER(k).sigr);
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
for kl = 1:length(restraints.pflex)
    for k = 1:length(restraints.pflex(kl).DEER)
        pflex = pflex + 1;
        if restraints.pflex(kl).DEER(k).r ~=0 &&  restraints.pflex(kl).DEER(k).sigr ~=0 % valid restraint
            n_DEER = n_DEER + 1;
            distr = restraints.pflex(kl).DEER(k).distr;
            if ~isempty(distr)
                distr = distr/sum(distr);
                data = zeros(length(rax),3);
                data(:,1) = rax.';
                argr = (restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*restraints.pflex(kl).DEER(k).sigr);
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

function [chi2_vec,fits,valid] = get_SAS(restraints)

chi2_vec = zeros(1,3);
fits = cell(1,3);
valid = true;

n_SANS = 0;
if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
    n_SANS = length(restraints.SANS);
    for ks = 1:n_SANS
        fits{ks} = restraints.SANS(ks).fit;
        chi2_vec(ks) = restraints.SANS(ks).chi2;
    end
end

if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS)
    n_SAXS = length(restraints.SAXS);
    for ks = 1:n_SAXS
        fits{n_SANS+ks} = restraints.SAXS(ks).fit;
        chi2_vec(n_SANS+ks) = restraints.SAXS(ks).chi2;
    end
end

