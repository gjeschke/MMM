function fom = fit_multi_DEER(v,fit,opt)
% Computes the distance distributions for an ensemble and their geometric
% mean overlap deficiency with respect to restraint distributions
%
% v     vector of populations for the ensemble members
% fit   cell vector of distance distribution information, each cell
%       contains data for one DEER restraint, an array [m,n], where
%       m   number of points in distance domain
%       n-2 number of ensemble members, must match length of v
%       1st column    : distance axis (not used here)
%       2nd column    : restraint distance distribution
%       columns 3... n: distance distributions for ensemble members
% opt   optional information for plot/display
%
% G. Jeschke, 3.7.2019

update_number = 1000; % Number of function calls before plot output is updated

persistent call_count
if isempty(call_count)
    call_count = 0;
end

fom = 1;
for k = 1:length(fit)
    data = fit{k};
    restraint = data(:,2);
    [~,n] = size(data);
    basis = data(:,3:n);
    coeff = v;
    sim = basis*coeff';
    sim = sim/sum(sim);
    restraint = restraint/sum(restraint);
    % overlap = real(sum(sqrt(restraint.*sim)));
    overlap = sum(min([restraint';sim']));
    fom = fom * overlap;
end

fom = 1 - fom^(1/length(fit));

if exist('opt','var')
    if mod(call_count,update_number) == 0
        call_count = 0;
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_ensemble_size,v(1:opt.old_ensemble_size),'o','Color',[0.6,0,0]);
        plot(opt.old_ensemble_size+1:length(v),v(opt.old_ensemble_size+1:length(v)),'.','Color',[0.6,0,0]);
        opt.text_DEER_fom0 = sprintf('%5.3f',fom);
        coeff = v(v > opt.threshold*max(v));
        opt.text_DEER_size = length(coeff);
        drawnow
    end
    call_count = call_count + 1;
end