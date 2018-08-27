function [texp,vexp,zt,dt]=pre_process_MMM(t0,v0,zt,dt)
%
% Preprocessing of DEER data
% 
% pre-processed data are shifted on the time axis to minimize
% the first moment, no phase correction
% simplified version for MMM
%
% when data reduction is required, information from all primary data points
% is included by smoothing, if requested data points fall in between
% primary data points, spline interpolation is used, data are normalized to
% the maximum of the real part
%
% Input parameters:
% t0    time axis of original experimental data
% v0    original experimental data (echo integral)
% zt    zero time shift (optional), if not given, an automatic
%       determination is performed that minimizes the first moment around
%       the maximum
% dt    time increment of the pre-processed data, must be a positive
%       multiple of eight (optional), defaults to the lowest positive
%       multiple of eight that leads to a data set with at most 1024 points
%
% Additional output parameters:
% texp  zero-shift corrected experimental time axis
% vexp  interpolated/smoothed experimental data


[m,n]=size(v0); % catches cases where 2D data without explicit summation over nuclear modulation are offered
% if m > 1
%     v0 = sum(v0);
% end

last=v0(length(v0));
phi0=atan2(imag(last),real(last));
im_off_0=0;
v=[phi0,im_off_0];
pstart=round(length(v0)/8); % use only last 7/8 of data for phase/offset correction
fit_offset=false;
if ~fit_offset,
    v=v(1);
end;
v=fminsearch(@rmsd_phi_offset,v,[],v0(pstart:end));
phi=v(1);
if fit_offset,
    imo=v(2);
else
    imo=0;
end;
if sum(real(v0*exp(1i*phi)))<0, phi=phi+pi; end;
vexp=(v0-1i*imo)*exp(1i*phi);
vexp=real(vexp)/max(real(vexp));


t1=min(t0):1:max(t0);
v1=interp1(t0,vexp,t1,'spline',vexp(1));
% get zero time, if not provided
if nargin<3,
    expdat=real(v1);
    % Determine maximum
    [maxi,mp]=max(expdat);
    nmp=1;
    % Determine time zero by moment analysis
    if mp>1 && mp<length(expdat),
        dmi=mp-1;
        dma=length(expdat)-mp;
        dm=dmi; if dma<dm, dm=dma; end;
        maxd=floor(dm/2);
        dev=1e20;
        nmp=mp;
        for k=-maxd+mp:maxd+mp,
            summ=0;
            for l=-maxd:maxd,
                summ=summ+expdat(k+l)*l;
            end;
            if abs(summ)<dev, dev=abs(summ); nmp=k; end;
        end;
    end;
    zt=t1(nmp);
end;

% Correct time axis
t0=t0-zt*ones(size(t0));
ta=t0(1); % real start time
te=t0(length(t0)); % real end time

dt0=t0(2)-t0(1);

% get dt if not provided
if nargin<4,
    dt=dt0;
    if dt<8, dt=8; end;
    if mod(dt,8)~=0,
        dt=8*ceil(dt/8);
    end;
end;
redfac0=dt/dt0; % data reduction factor for requested dt

max_nexp=1024; % no more than 1024 data points
if n>max_nexp,
    redfac=ceil(n/max_nexp);
    if redfac>redfac0, % requested dt is too small
        dt=redfac*dt0;
    end;
    if mod(dt,8)~=0, % next multiple of eight
        dt=8*ceil(dt/8);
    end;
    redfac=round(dt/dt0);
    % Data reduction with smoothing
    pas=pascal(redfac);
    b=zeros(1,redfac);
    for k=1:redfac,
        b(k)=pas(redfac-k+1,k);
    end;
    b=b/sum(b);
    z0=vexp;
    a=1;
    [m2,n2]=size(z0);
    z=zeros(m2,n2);
    for kk=1:m2,
        z0a=z0(kk,:);
        z0a=filter(b,a,z0a);
        z(kk,:)=z0a;
        z(kk,1)=z0(kk,1);
    end;
else
    z=vexp;
end;

% Make new time axis
mint=dt*ceil(ta/dt);
if mint<ta, mint=ta; end;
maxt=dt*floor(te/dt);
% while maxt>te, maxt=maxt-dt; end;
if maxt>te, maxt=te; end;
newn=1+floor((maxt-mint)/dt);
maxt=mint+dt*(newn-1);
texp=mint:dt:maxt;

% Interpolation onto new time axis
vexp=zeros(m,length(texp));
[m2,n2]=size(z);
for kk=1:m2,
    vexp(kk,:)=interp1(t0,z(kk,:),texp,'spline',z(kk,1));
end;

[mi,poi]=min(abs(texp));
texp=texp(poi:end);
vexp=vexp(poi:end);


function rmsd=rmsd_phi_offset(v,tr)
% Computes root mean square deviation of the imaginary part of
% phase-corrected data from zero
% before phase correction, an offset is subtracted from the imaginary part
%
% v(1)  phase correction phi (rad)
% v(2)  offset
% tr    complex data trace
%
% rmsd  root mean square deviation of imaginary part from zero
%
% G. Jeschke, 2009

if length(v)>1,
    tr=tr-i*v(2);
    itr=imag(tr*exp(i*v(1)));
    rmsd=sqrt(sum(itr.*itr));
else
    itr=imag(tr*exp(i*v(1)));
    rmsd=sqrt(sum(itr.*itr));
end;
