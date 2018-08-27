function [loop,x,y,z]=plot_loop(backbone,col,radius)
% [loop,x,y,z]=plot_loop(backbone,col)
%
% Plots ribbon model for a loop (or turn), can also be used for tube model
% the graphic output is sent to the current axes object
%
% backbone  backbone curve as obtained with get_ribbon
% col       color specifier, defaults to 'g' (green)
% radius    tube radius, defaults to 0.4 Angstroem
%
% loop      handles to the graphic objects (tube sections) that define the
%           loop or turn
% x,y,z     extension of graphics object (min and max of coordinates)
%
% (c) G. Jeschke, 2009

res=10; % ten faces

% Geometric parameters (in Angstroem)
if nargin<3,
    radius=0.2;
end;
falpha=1; % opaque

if nargin<2,
    col=[0,1,0]; % use green color if nothing is specified
end;

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;
[x,y,z]=tubeplot(backbone',radius/2);
loop=surf(x,y,z);
maxx=max([maxx,max(max(x))]);
minx=min([minx,min(min(x))]);
maxy=max([maxy,max(max(y))]);
miny=min([miny,min(min(y))]);
maxz=max([maxz,max(max(z))]);
minz=min([minz,min(min(z))]);
set(loop,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
