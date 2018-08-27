function [helix,x,y,z]=plot_helix(backbone,rung,normal,col)
% helix=plot_helix(backbone,rung,normal,spr,col)
%
% Makes ribbon model for an alpha-helix
% the graphic output is sent to the current axes object
%
% backbone  backbone curve as obtained with get_ribbon
% rung      backbone rungs as obtained with get_ribbon
% normal    normal vectors on ribbon as obtained with get_ribbon
% col       color specifier, defaults to 'g' (green)
%
% helix     handles to the graphic objects (tube sections) that define the
%           helix
%
% (c) G. Jeschke, 2009

falpha=1; % opaque

if nargin<4,
    col=[1,0,0]; % use red color if nothing is specified
end;

[m,n]=size(backbone);

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;

[x,y,z]=generalized_tubeplot(backbone',rung',normal',1);
helix=surf(x,y,z);
set(helix,'EdgeColor','none','FaceColor',col','FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off'); % 
maxx=max([maxx,max(max(x))]);
minx=min([minx,min(min(x))]);
maxy=max([maxy,max(max(y))]);
miny=min([miny,min(min(y))]);
maxz=max([maxz,max(max(z))]);
minz=min([minz,min(min(z))]);


x=[minx,maxx];
y=[miny,maxy];
z=[minz,maxz];