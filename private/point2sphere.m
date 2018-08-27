function [x,y,z] = point2sphere(coor,radius)
%VEC2TUBE	Create parametric spherical surface at point with coordinates coor
%	[x,y,z] = point2sphere(coor,radius).
%	

%	(c) 2007 by Gunnar Jeschke
%

[x,y,z]=sphere;
[m,n]=size(x);
x=radius*x;
y=radius*y;
z=radius*z;

x=x+coor(1)*ones(m,n);
y=y+coor(2)*ones(m,n);
z=z+coor(3)*ones(m,n);


