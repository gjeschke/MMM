function [fom,fits] = sim_multi_DEER(v,fit)

fom = 1;
fits(length(fit)).rax = [];
fits(length(fit)).restraint = [];
fits(length(fit)).distr = [];
fits(length(fit)).overlap = [];
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
    fits(k).rax = data(:,1);
    fits(k).restraint = restraint;
    fits(k).distr = sim;
    fits(k).overlap = overlap;
end

fom = 1 - fom^(1/length(fit));