function fom = fit_multi_SAS(v,fit)

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