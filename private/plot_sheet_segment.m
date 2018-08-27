function [sheet,x,y,z]=plot_sheet_segment(backbone,rung,normal,col,Cterm)
% sheet_handles=plot_sheet_segment(backbone,rung,normal,col,Cterm)
%
% Plots ribbon model for a beta-sheet segment
% the graphic output is sent to the current axes object
%
% backbone  backbone curve as obtained with get_ribbon
% rung      backbone rungs as obtained with get_ribbon
% normal    normal vectors on ribbon as obtained with get_ribbon
% col       color specifier, defaults to 'g' (green)
% Cterm     flag that requires C-terminal arrow if 1, no arrow if 0
%
% sheet     handles to the graphic objects (tube sections) that define the
%           sheet
% x,y,z     extension of graphics objects, x=min(x),max(x) ...
%
% (c) G. Jeschke, 2009

global graph_settings

% Geometric parameters (in Angstroem)
sheet_width=2;
sheet_arrow=2.25;
sheet_height=0.3;
falpha=1; % opaque

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;

if nargin<4,
    col=[0,0,1]; % use blue color if nothing is specified
end;

if Cterm
    [x,y,z]=generalized_tubeplot(backbone',rung',normal',2,graph_settings.spr);
else
    [x,y,z]=generalized_tubeplot(backbone',rung',normal',2);
end;

sheet=surf(x,y,z);
set(sheet,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
maxx=max([maxx,max(max(x))]);
minx=min([minx,min(min(x))]);
maxy=max([maxy,max(max(y))]);
miny=min([miny,min(min(y))]);
maxz=max([maxz,max(max(z))]);
minz=min([minz,min(min(z))]);


x=[minx,maxx];
y=[miny,maxy];
z=[minz,maxz];