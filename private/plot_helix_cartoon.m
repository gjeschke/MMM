function [helix,x,y,z,spine]=plot_helix_cartoon(CA,col)
% helix=plot_helix_cartoon(CA,col)
%
% Makes cartoon model for an alpha-helix
% the graphic output is sent to the current axes object
%
% CA        C-alpha coordinates
% col       color specifier, defaults to 'm' (magenta)
%
% helix     handles to the graphic objects (tube sections) that define the
%           helix
%
% (c) G. Jeschke, 2014

global graph_settings

falpha=1; % opaque
contract = 1;

if nargin<2,
    col=[1,0,1]; % use magenta color if nothing is specified
end;

[p0,v]=rmsd_line_3D(CA);
% v = contract*v;
% p0 = p0 + (1-contract)*v;

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;

CA1 = p0-v/2;
if norm(CA1-CA(1,:)) < norm(CA1-CA(end,:)),
    spine = [p0-v/2;p0+v/2];
else
    spine = [p0+v/2;p0-v/2];
end;    
[x,y,z]=tubeplot(spine',graph_settings.helix_cartoon_radius,graph_settings.helix_cartoon_points);
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