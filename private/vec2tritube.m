function [t,x,y,z] = vec2tritube(coor1,coor2,radius)
%VEC2TUBE	Create triangulated cylindrical tube surface around the vector 
%           from coor1 to coor2
%	[x,y,z] = vec2tube(coor1,coor2,radius).
%
% uses MyRobustCrust by Luigi Giaccari 
%	
%
%	(c) 2009 by Gunnar Jeschke
%

n=10;

vec=coor2-coor1;
[fi,th,lv]=cart2sph(vec(1),vec(2),vec(3));
th=pi/2-th;

x=radius*cos(linspace(0,2*pi*(n-1)/n,n));
y=radius*sin(linspace(0,2*pi*(n-1)/n,n));
p=zeros(2*n+2,3);
p(1,:)=[0,0,0];
p(2:1+n,1)=x';
p(2:1+n,2)=y';
p(n+2:2*n+1,1)=x';
p(n+2:2*n+1,2)=y';
p(n+2:2*n+1,3)=lv*ones(n,1);
p(2*n+2,:)=[0,0,lv];

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

rmat=eth*efi;

for k=1:2*n+2,
    p(k,:)=p(k,:)*rmat+coor1;
end;

t = delaunayn(p);

x=p(:,1);
y=p(:,2);
z=p(:,3);