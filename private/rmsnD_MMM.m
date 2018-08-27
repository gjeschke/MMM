function rms = rmsnD_MMM(v,texp,vexp,ff,t0,logB,mod_depth)
%rmsnd_MMM	Root mean square error of function exp(-k*t^(n/3)).
%	rms = rmsnD_MM(v,texp,vexp,ff,logB)
%	
%  Parameter: v(1) Modulation depth, v(2) time constant

%	Copyright (c) 2009 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;
if length(v)>1,
    if v(2)<0, rms=1.0e6; return; end;
    if v(1)>1, rms=1.0e6; return; end;
    mod_depth=v(1);
    v=v(2);
end;
bckg=decaynD_MMM(v,texp,t0,logB);

deconvoluted=mod_depth*ff+(1-mod_depth)*ones(size(ff));
sim=deconvoluted.*bckg;

diff=sim-vexp;
rms=sqrt(sum(diff.^2)/length(diff));
