function rms = rms_stretched_MMM(v,texp,vexp,ff,mod_depth)
%RMS_stretched_MMM	Root mean square error of function exp(-k*t^ksi).
%	rms = rms_stretched_MMM(v,texp,vexp,ff).
%	
%  Parameter: v(1) Modulation depth, v(2) time constant k, v(3) ksi
%
% mod_depth is an optional argument, if present, depth is not fitted

%	Copyright (c) 2009 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;
if v(2)<0, rms=1.0e6; return; end;
if length(v)>2,
    if v(3)<0, rms=1.0e6; return; end; 
    if v(1)>1, rms=1.0e6; return; end;
end;

if nargin>4,
    bckg=decay_stretched_MMM(v,texp);
else
    bckg=decay_stretched_MMM(v(2:3),texp);
    mod_depth=v(1);
end;
deconvoluted=mod_depth*ff+(1-mod_depth)*ones(size(ff));
sim=deconvoluted.*bckg;

diff=sim-vexp;
rms=sqrt(sum(diff.^2)/length(diff));

		

