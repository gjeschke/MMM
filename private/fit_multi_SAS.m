function fom = fit_multi_SAS(v,fit,opt)
% Computes the SAXS and SANS curves for an ensemble and their sum of
% chi^2 with respect to experimental curves
%
% v     vector of populations for the ensemble members and of constant
%       offsets for all SAS curves
% fit   cell vector of SAS curve information, each cell
%       contains data for one SAS curve, an array [m,n], where
%       m   number of points in scattering vector domain
%       n-3 number of ensemble members, must match length of v - number of
%           curves
%       1st column      : scattering vector axis (not used here)
%       2nd column      : scattering intensity
%       3rd column      : errors of scattering intensity
%       columns 4... n-3: simulated SAS curves for ensemble members
% opt   optional information for plot/display
%
% G. Jeschke, 3.7.2019

update_number = 1000; % Number of function calls before plot output is updated

persistent call_count
if isempty(call_count)
    call_count = 0;
end

fom = 0;
for k = 1:length(fit)
    data = fit{k};
    curve = data(:,2);
    errors = data(:,3);
    [m,n] = size(data);
    basis = data(:,4:n);
    coeff = v(1:n-3);
    sim = basis*coeff' + v(n-3+k);
    sc = sum(curve.*sim)/sum(sim.*sim);
    chi2 = sum(((curve-sc*sim)./errors).^2)/(m-1);
    fom = fom + chi2;
end

if exist('opt','var')  && isfield(opt,'interactive') && opt.interactive
    if mod(call_count,update_number) == 0
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_ensemble_size,v(1:opt.old_ensemble_size),'o','Color',[0,0,0.75]);
        plot(opt.old_ensemble_size+1:n-3,v(opt.old_ensemble_size+1:n-3),'.','Color',[0,0,0.75]);
        opt.text_SAS_fom0.String = sprintf('%5.3f',fom);
        coeff = v(v > opt.threshold*max(v));
        opt.text_SAS_size.String = sprintf('%i',length(coeff));
        opt.text_iteration_counter.String = sprintf('%i',call_count/1000);
        drawnow
    end
    call_count = call_count + 1;
end