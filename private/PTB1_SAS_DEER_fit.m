function PTB1_SAS_DEER_fit(fname)

SAS_fit_file = sprintf('%s_SAS_fit.mat',fname);
DEER_fit_file = sprintf('%s_DEER_fit.mat',fname);
definition_file = sprintf('%s_ensemble_definition.mat',fname);

SAS_fits = load(SAS_fit_file);
normalize.SAS = SAS_fits.fom;
all_fits_SAS = SAS_fits.all_fits;

% load PTB1_DEER_coeff
DEER_fits = load(DEER_fit_file);
ensemble_size = length(DEER_fits.coeff);
normalize.DEER = DEER_fits.fom1;
all_fits_DEER = DEER_fits.all_fits;

fanonym=@(v_opt)fit_multi_SAS_DEER(v_opt,all_fits_SAS,all_fits_DEER,normalize);

l=[zeros(1,ensemble_size) -0.005 -0.005 -0.01];
u=[ones(1,ensemble_size) 0.005 0.005 0.01];
v0 = [ones(1,ensemble_size)/ensemble_size 0 0 0];
options = optimoptions('patternsearch','MaxFunctionEvaluations',50000*length(v0),'MaxIterations',2500*length(v0));

tic,
[v,fom,exitflag,output] = patternsearch(fanonym,v0,[],[],[],[],l,u,[],options);
toc,

figure(1); clf;
plot(v(1:end-3),'k.');

coeff1 = v(1:end-3);
coeff2 = coeff1/max(coeff1);
coeff2 = coeff2(coeff2 >= 1e-2);
fprintf(1,'%i populations are at least 1%% of the maximum population\n',length(coeff2)); 

fprintf(1,'Fit with quality loss %6.4f\n',fom);

[fom_SAS,fits] = sim_multi_SAS(v,all_fits_SAS);

fprintf(1,'Total SAS chi2 %6.4f\n',fom_SAS);


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


[fom_DEER,fits_DEER] = sim_multi_DEER(v(1:end-3),all_fits_DEER);

fprintf(1,'DEER overlap deficiency %6.4f\n',fom_DEER);

save(definition_file,'v','fom','fom_SAS','fom_DEER','fits','fits_DEER','exitflag','output');

n_DEER = length(fits_DEER);
for k = 1:n_DEER
    figure(k+3); clf;
    plot(fits_DEER(k).rax,fits_DEER(k).restraint,'g');
    hold on;
    plot(fits_DEER(k).rax,fits_DEER(k).distr,'r');
    title(sprintf('DEER fit with overlap %5.3f',sum(sqrt(fits_DEER(k).restraint.*fits_DEER(k).distr))));
end
