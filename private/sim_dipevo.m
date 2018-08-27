function deer=sim_dipevo(v,ff,logB)
%
% Simulation of a DEER time trace from the form factor and
% a two-parameter vector for three-dimensional background
%
% v	  background parameters: v(1)=Delta  modulation depth
%                            v(2)=k      decay constant
% ff      form factor for the given distance distribution
% logB    logarithm of the background function

bckg=exp(v(2)*16000*logB);
deconvoluted=v(1)*ff+(1-v(1))*ones(size(ff));
deer=deconvoluted.*bckg;


