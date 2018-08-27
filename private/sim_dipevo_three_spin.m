close all
clear all

if ispc
    bas_path='c:/users/jeschke/Matlab_libs/';
    addpath([bas_path 'easyspin-3.0.0/easyspin/']);
else
    bas_path='/home/guje/';
end;

nKnots=500; % 150 number of knots for orientational grid (EasySpin)
[phi,theta] = sphgrid('Ci',nKnots); % orientational grid (EasySpin)

ra=4.16; % N-O to N-O distances
rb=4.16;
rc=4.16;
sigr=0.190; % standard deviation of distance
s=(ra+rb+rc)/2;
thbc0=2*asin(sqrt((s-rb)*(s-rc)/(rb*rc))); % Halbwinkelsätze
thac0=2*asin(sqrt((s-ra)*(s-rc)/(ra*rc)));
thab0=2*asin(sqrt((s-ra)*(s-rb)/(ra*rb)));

disp('Angles:');
disp(sprintf('alpha= %5.1f°',180*thbc0/pi));
disp(sprintf('beta = %5.1f°',180*thac0/pi));
disp(sprintf('gamma= %5.1f°',180*thab0/pi));

dt=0.016;
ndat=374;
pair=zeros(1,ndat);
triple=zeros(1,ndat);

t=linspace(0,(ndat-1)*dt,ndat); % generate time axes


% tc=max(t1)/2;
% dmptr=t1/tc;
% dmptr=exp(-dmptr.^2);
% dmp=dmptr'*dmptr;

% powder average
for k=1:length(phi),
    if mod(k,1000)==0, disp(k/length(phi)); end;
    st=sin(theta(k));
    ct=cos(theta(k));
    cf=cos(phi(k));
    ras=ra+randn*sigr;
    rbs=rb+randn*sigr;
    rcs=rc+randn*sigr;
    s=(ras+rbs+rcs)/2;
    thbc=2*asin(sqrt((s-rbs)*(s-rcs)/(rbs*rcs))); % Halbwinkelsätze
    thac=2*asin(sqrt((s-ras)*(s-rcs)/(ras*rcs)));
    thab=2*asin(sqrt((s-ras)*(s-rbs)/(ras*rbs)));
	nydda=52.04/ras^3; % dipolar frequency
	wa=2*pi*nydda; % angular dipolar frequency
	nyddb=52.04/rbs^3; % dipolar frequency
	wb=2*pi*nyddb; % angular dipolar frequency
	nyddc=52.04/rcs^3; % dipolar frequency
	wc=2*pi*nyddc; % angular dipolar frequency
    wa1=wa*(3*ct^2-1);
    wb1=wa*(3*ct^2-1);
    wc1=wa*(3*ct^2-1);
	sab=sin(thab);
	cab=cos(thab);
	sac=sin(thac);
	cac=cos(thac);
	sbc=sin(thbc);
	cbc=cos(thbc);
    xab=sab*st*cf+cab*ct; % Jeschke/Bittl, CPL 294, 323-331 (1998), before eqn. (17)
    xac=sac*st*cf+cac*ct; % Jeschke/Bittl, CPL 294, 323-331 (1998), before eqn. (17)
    xbc=sbc*st*cf+cbc*ct; % Jeschke/Bittl, CPL 294, 323-331 (1998), before eqn. (17)
	wab2=wa*(3*xab^2-1);
	wac2=wa*(3*xac^2-1);
	wba2=wb*(3*xab^2-1);
	wbc2=wb*(3*xbc^2-1);
	wca2=wc*(3*xac^2-1);
	wcb2=wc*(3*xbc^2-1);
	tr1=cos(wa1*t);
	pair=pair+tr1;
	tr2=cos(wab2*t);
    pair=pair+tr2;
	triple=triple+tr1.*tr2;
	tr2=cos(wac2*t);
    pair=pair+tr2;
	triple=triple+tr1.*tr2;
	tr1=cos(wb1*t);
	pair=pair+tr1;
	tr2=cos(wba2*t);
    pair=pair+tr2;
	triple=triple+tr1.*tr2;
	tr2=cos(wbc2*t);
    pair=pair+tr2;
	triple=triple+tr1.*tr2;
	tr1=cos(wc1*t);
	pair=pair+tr1;
	tr2=cos(wca2*t);
    pair=pair+tr2;
	triple=triple+tr1.*tr2;
	tr2=cos(wcb2*t);
    pair=pair+tr2;
	triple=triple+tr1.*tr2;
end;
pair=pair/max(pair);
triple=triple/max(triple);

% z=dmp.*z;

awin=apowin('ham+',ndat);
pair1=awin'.*pair;
pair1(1)=pair1(1)/2;
pair1=[pair1 zeros(1,1024-ndat)];
triple1=awin'.*triple;
triple1(1)=triple1(1)/2;
triple1=[triple1 zeros(1,1024-ndat)];

pair_spc=fftshift(real(fft(pair1)));
triple_spc=fftshift(real(fft(triple1)));

fmin=-length(t)/(2*max(t));          % s. Schweiger/Jeschke book p. 106
fmax=(length(t)-1)/(2*max(t));
ny=linspace(fmin,fmax,length(pair_spc));

pair_spc=pair_spc/max(max(pair_spc));
triple_spc=triple_spc/max(max(triple_spc));

figure
plot(ny,pair_spc,'k');
hold on;
plot(ny,triple_spc,'r');

axis([-8,8,-0.05,1.05]);


save sim_corr_TR111_rigid ny pair_spc triple_spc t pair triple