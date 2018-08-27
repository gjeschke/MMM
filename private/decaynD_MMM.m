function sim=decaynD_MMM(v,texp,t0,logB)
% function sim=decaynD_MMM(v,x)
%
% Computes DEER decay for n-dimensional homogeneous distribution
% type of distribution is provided in argument logB
%
% v(1)  decay constant
% texp  abscissa values for result (vector)
% t0    abscissa values corresponding to logb
%
sim0=exp(16000*v(1)*logB);
sim=interp1(t0,sim0,texp,'pchip');

