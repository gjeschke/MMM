function fom = fit_multi_SAS_DEER(v,all_fits_SAS,all_fits_DEER,normalize)

fom_SAS = fit_multi_SAS(v,all_fits_SAS);
fom_DEER = fit_multi_DEER(v(1:end-3),all_fits_DEER);

fom = (fom_SAS/normalize.SAS + fom_DEER/normalize.DEER)/2 - 1;