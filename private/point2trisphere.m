function [x,y,z,t]=point2trisphere(coor,radius,n)
% Creates a triangulated sphere centered at coordinates coor with radius r
% function is just a wrapper for BuildSphere by Luigi Giaccari
%
% coor      [x0,y0,z0] coordinate vector of the center of the sphere
% radius    radius of the sphere
% n         optional quality number, between 0 and inf, be careful, number
%           of points grows dramatically with increasing n, defaults to 0 
%
% x,y,z     coordinates of triangle vertices
% t         face matrix for calling trisurf, use as
%           trisurf(t,p(:,1),p(:,2),p(:,3))
%
% G. Jeschke, 2009

if nargin<3,
    n=1;
end;

[p,t]=BuildSphere(n);

p=radius*p;
offset=ones(size(p(:,1)));
x=p(:,1)+coor(1)*offset;
y=p(:,2)+coor(2)*offset;
z=p(:,3)+coor(3)*offset;

