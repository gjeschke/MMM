function zt=get_zerotime_MMM(handles,t,expdat)
%
% Estimates zero time by moment analysis
%

t1=min(t):1:max(t);
v1=interp1(t,expdat,t1,'spline',expdat(1));
v1=v1/max(v1);

% Determine maximum
mp=0;
maxi=-1e20;
for k=1:length(v1),
    if v1(k)>maxi, maxi=v1(k); mp=k; end;
end;

nmp=1;
% Determine time zero by moment analysis
if mp>1 && mp<length(v1),
    dmi=mp-1;
    dma=length(v1)-mp;
    dm=dmi; if dma<dm, dm=dma; end;
    maxd=floor(dm/2);
    dev=1e20;
    nmp=mp;
    for k=-maxd+mp:maxd+mp,
        summ=sum(v1(k-maxd:k+maxd).*(-maxd:maxd));
        if abs(summ)<dev, dev=abs(summ); nmp=k; end;
    end;
end;
zt=t1(nmp)-t1(1);
zt_string=sprintf('%d',zt);
set(handles.zt_edit,'String',zt_string);
handles.zerotime=zt;

% Update handles structure
guidata(handles.main_figure, handles);