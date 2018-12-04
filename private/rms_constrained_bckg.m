function rms = rms_constrained_bckg(v,texp,vexp,ff,constraints)
% rms = rms_constrained_bckg(v,texp,vexp,ff,constraints)
%	
%  Parameter: v(1) rate constant k, v(2) modulation depth, v(3) dim
%
% the parameter vector can have one, two, or three arguments
% 
% v             parameter vector
% texp          experimental time axis
% vexp          experimental data
% ff            form factor
% constraints   array (3,2) with lines
%               minimum decay rate, maximum decay rate
%               minimum modulation depth, maximum modulation depth
%               minimum dimension, maximum dimension
%
%	G. Jeschke, 22.11.2018

for k = 1:length(v)
    if v(k) < constraints(k,1) || v(k) > constraints(k,2)
        rms = 1.0e6;
        return
    end
end

if length(v) < 3
    v(3) = constraints(3,1);
end

if length(v) < 2
    v(2) = constraints(2,1);
end

if v(1)<0, rms=1.0e6; return; end
if v(2)<0, rms=1.0e6; return; end
if v(3)<0, rms=1.0e6; return; end
if v(2)>1, rms=1.0e6; return; end

bckg = decay_stretched_MMM([v(1) v(3)],texp);
deconvoluted  = v(2)*ff+(1-v(2))*ones(size(ff));
sim = deconvoluted.*bckg;
sc = sum(sim.*sim)/sum(sim.*vexp);
vexp = sc*vexp;

diff=sim-vexp;
rms=sqrt(sum(diff.^2)/sum(vexp.^2));

		

