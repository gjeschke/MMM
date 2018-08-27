function [backbone,rung,normal]=mk_ribbon(spr,Ca_coor,O_coor,sheetflag)
% function [backbone,rung,normal]=mk_ribbon(Ca_coor,O_coor,sheetflag)
%
% Create coordinates of a ribbon guide curve for a protein backbone segment,
% using a slightly modified version of the Carson-Bugg algorithm (J. Mol. 
% Graph. 1986, 4, 121-122), unit vectors in and normal to the ribbon face
% are also supplied
%
% uses: function bezierInterp
%
% spr                   segments per residue
% Ca_coor               coordinates of Calpha atoms
% O_coor                optional, coordinates of peptide oxygen atoms
% sheetflag             0 for loops and alpha-helices, 1 for sheets,
%                       defaults to 0
%
% backbone              coordinates of backbone guide points
% rung                  unit vectors within peptide plane, perpendicular to
%                       backbone, is empty, if O_coor is missing
% normal                unit vectors perpendicular to both peptide plane
%                       and backbone, is empty, if O_coor is missing
% 
% (c) G. Jeschke, 2009

backbone=[];
rung=[];
normal=[];

Tension=0; % assume cardinal splines with zero tension

if nargin<4,
    sheetflag=0;
end;

[m,n]=size(Ca_coor);

if m<2,
    return;
end;

guide_length=(m-1)*spr+1;
backbone=zeros(guide_length,3);
backbone0=[Ca_coor(1,:); Ca_coor; Ca_coor(m,:)];

if nargin>2,
    rung0=zeros(m,n);
    dirsign=1;

    for k=1:m-1,
        a=Ca_coor(k+1,:)-Ca_coor(k,:);
        a=a/norm(a);
        b=O_coor(k,:)-Ca_coor(k,:);
        c=cross_rowvec(a,b);
        c=c/norm(c); % normal to peptide plane
        d=cross_rowvec(c,a);
        d=d/norm(d); % y axis in peptide plane
        % if sheetflag, dirsign=dirsign*-1; end; % flip side, if necessary
        if k>1,
            dirsign=dirsign*sum(d.*d0);
        end;
        rung0(k,:)=dirsign*d;
        d0=d;
    end;
    rung0(m,:)=dirsign*d;
    normal=backbone;
    rung=backbone;
    rung0=[rung0(1,:); rung0; rung0(m,:)];    
end;


for k=1:m-1,
    basnum=(k-1)*spr;
       backbone(basnum+1:basnum+spr+1,:)=cardinal_spline(backbone0(k:k+3,:),Tension,spr);
       backbone(basnum+1:basnum+spr+1,:)=cardinal_spline(backbone0(k:k+3,:),Tension,spr);
    if nargin>2,
       rung(basnum+1:basnum+spr+1,:)=cardinal_spline(rung0(k:k+3,:),Tension,spr);    
    end;
end;
if nargin>2,
    normal=zeros(size(rung));
    for k=1:guide_length-1,
        e=backbone(k+1,:)-backbone(k,:);
        f=rung(k,:);
        g=cross_rowvec(e,f);
        normal(k,:)=g/norm(g);
        rung(k,:)=f/norm(f);
    end;
    normal(guide_length,:)=normal(guide_length-1,:);
    rung(guide_length,:)=rung(guide_length,:)/norm(rung(guide_length,:));
end;

if sheetflag,
    backbone=smooth_backbone(backbone,spr);
end;

function backbone=smooth_backbone(backbone,spr)
% Polynomial fitting (and smoothing) of the backbone curve

[m,n]=size(backbone);
mb=floor(m/spr);
if m<20; return; end;
x=1:m;
xb=linspace(1,m,mb);
xb=xb';
x=x';
for k=1:3,
    y=backbone(:,k);
    xd=[1,m,1,m];
    d=[y(1),y(m),y(2)-y(1),y(m)-y(m-1)];
    js=[0,0,1,1];
    pp=spfit(x,y,xb,3,xd,d,js);
    smoothed=ppval(pp,xb);
    if ~sum(isnan(smoothed)), 
        backbone(:,k)=interp1(xb,smoothed,x,'spline');
    end;
end;
