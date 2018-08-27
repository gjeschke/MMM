function rms=rms_dipevo(v,ff,logB,t,texp,vexp),
%
% Computes root mean square deviation of DEER data simulated
% for form factor ff and background parameters v from experimental 
% DEER data vexp given at times texp
% the logarithm of the background function logB and the time axis t
% of the kernel must be given
% v(1)  modulation depth Delta
% v(2)  decay constant k, is multiplied mit logB
%
% (c) G. Jeschke, 2006
%

n=length(texp);
deer0=sim_dipevo(v,ff,logB);
deer=interp1(t,deer0,texp,'pchip');
sc=sum(vexp.*vexp)/sum(vexp.*deer);
diff=sc*deer-vexp;
rms=sqrt(sum(diff.*diff)/(n-1));
