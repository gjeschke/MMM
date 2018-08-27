function obj=arrow(p1,p2,color,radius,alpha)
% function obj=arrow(p1,p2,color,radius,alpha)
%
% Plots a capped tube withh arrow pointer from point p1 to point p2
%
% p1        Cartesian coordinates of point 1
% p2        Cartesian coordinates of point 2
% color     RGB color triple, defaults to [0,0,1] (blue)
% radius    radius, defaults to graph_settings.stick_radius
% alpha     transparency value from 0 (invisible) to 1 (opaque), defaults
%           to 1
%
% obj       handles of three graphics objects (cap spheres and tube
%           cylinder)
%
% G. Jeschke, 2010

global graph_settings

if nargin<3,
    color=[0,0,1];
end;

if nargin<4,
    radius=graph_settings.stick_radius;
end;

if nargin<5,
    alpha=1;
end;

[x,y,z,t]=point2trisphere(p1,radius);
obj1=trisurf(t,x,y,z);
set(obj1, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha',alpha,'FaceLighting','gouraud','Clipping','off');
set(obj1, 'CDataMapping','direct','AlphaDataMapping','none');
[x,y,z] = vec2cone(p1,p2,radius,10);
obj2=surface(x,y,z);
set(obj2, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha',alpha,'FaceLighting','gouraud','Clipping','off');
set(obj2, 'CDataMapping','direct','AlphaDataMapping','none');
obj=[obj1,obj2];
end