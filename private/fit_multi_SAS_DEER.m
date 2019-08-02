function fom = fit_multi_SAS_DEER(v,all_fits_SAS,all_fits_DEER,normalize,opt)
% Computes integrative fit loss for combined fitting of populations in an
% ensemble that should fulfill DEER and SAXS/SANS restraints
%
% v     vector of populations for the ensemble members and of constant
%       offsets for all SAS curves
% all_fits_SAS
%       cell vector of SAS curve information, each cell
%       contains data for one SAS curve, an array [m,n], where
%       m   number of points in scattering vector domain
%       n-3 number of ensemble members, must match length of v - number of
%           curves
%       1st column      : scattering vector axis (not used here)
%       2nd column      : scattering intensity
%       3rd column      : errors of scattering intensity
%       columns 4... n-3: simulated SAS curves for ensemble members
% all_fits_DEER
%       cell vector of distance distribution information, each cell
%       contains data for one DEER restraint, an array [m,n], where
%       m   number of points in distance domain
%       n-2 number of ensemble members, must match length of v
%       1st column    : distance axis (not used here)
%       2nd column    : restraint distance distribution
%       columns 3... n: distance distributions for ensemble members
% normalize
%       best-fit figures of merit for individual restraint sets, structure
%       .SAS    sum of chi^2 for all SAXS/SANS curves
%       .DEER   geometric mean overlap deficiency for all distance
%               distributions
% opt   optional information for plot/display
%
% G. Jeschke, 3.7.2019

update_number = 1000; % Number of function calls before plot output is updated

persistent call_count
if isempty(call_count)
    call_count = 0;
end

data = all_fits_SAS{1};
[~,n] = size(data);

fom_SAS = fit_multi_SAS(v,all_fits_SAS);
fom_DEER = fit_multi_DEER(v(1:n-3),all_fits_DEER);

fom = (fom_SAS/normalize.SAS + fom_DEER/normalize.DEER)/2 - 1;

if exist('opt','var')
    if mod(call_count,update_number) == 0
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_ensemble_size,v(1:opt.old_ensemble_size),'o','Color',[0,0.6,0]);
        plot(opt.old_ensemble_size+1:n-3,v(opt.old_ensemble_size+1:n-3),'.','Color',[0,0.6,0]);
        opt.text_DEER_fom.String = sprintf('%5.3f',fom_DEER);
        opt.text_SAS_fom.String = sprintf('%5.3f',fom_SAS);
        opt.text_loss.String = sprintf('%5.3f',fom);
        coeff = v(v > opt.threshold*max(v));
        opt.text_ensemble_size.String = sprintf('%i',length(coeff));
        opt.text_iteration_counter.String = sprintf('%i',call_count/1000);
        drawnow
    end
    call_count = call_count + 1;
end