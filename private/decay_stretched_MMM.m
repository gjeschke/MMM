function sim=decay_stretched_MMM(v,x)
% function sim=decay_stretched_MMM(v,x)
%
% Computes stretched exponential decay function
%
% v(1)  decay constant
% v(2)  3 ksi, where ksi is the stretch exponent
% x     abscissa values (vector)
%
sim=exp(-v(1)*x.^(v(2)/3));
