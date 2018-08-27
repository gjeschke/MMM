function [x,y,z] = vec2tube(coor1,coor2,radius,res)
%VEC2TUBE	Create parametric tube surface around the vector from coor1 to coor2
%	[x,y,z] = vec2tube(coor1,coor2,radius).
%	

%	(c) 2000 by Gunnar Jeschke
%

if nargin<4,
    res=20;
end;

vec=coor2-coor1;
[fi,th,lv]=cart2sph(vec(1),vec(2),vec(3));
th=pi/2-th;
[x,y,z]=cylinder([1,1],res);
[m,n]=size(x);
z=lv*z;
x=radius*x;
y=radius*y;

eth=zeros(3,3);
eth(1,1)=cos(th);
eth(1,3)=-sin(th);
eth(3,1)=sin(th);
eth(3,3)=cos(th);
eth(2,2)=1;

efi=zeros(3,3);
efi(1,1)=cos(fi);
efi(1,2)=sin(fi);
efi(2,1)=-sin(fi);
efi(2,2)=cos(fi);
efi(3,3)=1;

for k=1:m,
   for l=1:n,
      act=[x(k,l) y(k,l) z(k,l)];
      act=act*eth;
      act=act*efi;
      x(k,l)=act(1);
      y(k,l)=act(2);
      z(k,l)=act(3);
   end;
end;

x=x+coor1(1)*ones(m,n);
y=y+coor1(2)*ones(m,n);
z=z+coor1(3)*ones(m,n);


