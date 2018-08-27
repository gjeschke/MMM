function ff=get_formfactor_MMM(distr,kernel,t)
% computes DEER formfactor from distance distribution distr,
% a kernel matrix, and the t axis of the kernel matrix
% care about modulation depth scaling is taken
% (c) G. Jeschke, 2006

distr=0.01*distr/sum(distr);
ff=100*(exp(distr*kernel)-ones(size(t)))+ones(size(t));