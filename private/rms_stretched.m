function rms = rms_stretched(v,x,y)
%RMS_streched	Root mean square error of function exp(-k*t^ksi).
%	rms = rms_stretched(v,x,y).
%	
%  Parameter: v(2) Amplitude, v(1) Zeitkonstante, v(3) ksi

%	Copyright (c) 2004 by Gunnar Jeschke

if v(1)<0, rms=1.0e6; return; end;
sim=decay_stretched(v,x);
diff=sim-y;
rms=sum(diff.*diff);

		

