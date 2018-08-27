function [intersect,hit]=intersect_triangle(ray,triangle,n,hit_test)
% function intersect=intersect_triangle(ray,triangle,n)
%
% determines whether a ray intersects a triangle
%
% if optional input argument n is missing, cross_colvec.m is required
% only if n is missing, a test for degeneracy of the triangle is done,
% otherwise calling routine is responsible for triangle not being
% degenerate
%
% for maximum speed call with three input parameters, including n 
% 
% Input:
% ray       ray specified by a 2*3 array of xyz ccordinates of the origin
%           ray(1,:) and a second point ray(2,:) on the ray, ray does not
%           end at ray(2,:), but extends to infinity
% triangle  triangle specified by a 3*3 array of xyz coordinates
%           triangle(i,:) where i=1...3 are point indices
% n         n optional face normal of the triangle, provided for speed
%
% Output:
% intersect flag that is 0 if the ray does not intersect the triangle, 1 if
%           it does intersect, 2 if it is in triangle plane, 
%           -1 if the triangle is degenerate
% hit       indices of points that were hit by the ray, empty if ray
%           intersects inside triangle or does not intersect, length 1 for
%           vertex hit, length 2 for edge hit, length 3 for ray coinciding
%           with triangle plane
% hit_test  optional falg to suppress hit test (for speed), defaults to 1
%           (hit test is performed), if 0, hit is always empty
%
% based on C++ code by SoftSurfer
% Copyright 2001, softSurfer (www.softsurfer.com)
% This code may be freely used and modified for any purpose
% providing that this copyright notice is included with it.
% SoftSurfer makes no warranty for this code, and cannot be held
% liable for any real or imagined damage resulting from its use.
% Users of this code must verify correctness for their application.
%
% Matlab adaptation and handling for full, vertex and edge hits by
% G. Jeschke, 2010
%#eml

u=zeros(1,3);
v=zeros(1,3);
dir=zeros(1,3);

hit=[];

u=triangle(2,:)-triangle(1,:);
v=triangle(3,:)-triangle(1,:);

if nargin<3, % if face normal was not provided, compute it and check for degeneracy
    n=cross_colvec(u,v); % do not use Matlab cross, it's very slow!
    if norm(n)<eps,
        intersect=-1;
        return;
    end;
end;

if nargin<4,
    hit_test=1;
end;

% axes(hMain.axes_model);
% plot3(ray(:,1),ray(:,2),ray(:,3),'r','LineWidth',1);
dir=ray(2,:)-ray(1,:);
w0=ray(1,:)-triangle(1,:);
a=sum(n.*w0);
b=sum(n.*dir);
if abs(b)<eps, % ray parallel to triangle plane
    if a==0,
        intersect=2; % ray coincides with triangle plane
        if hit_test,
            hit=[1,2,3];
        end;
        return;
    else
        intersect=0; % ray disjoint from plane
        return;
    end;
end;

r=a/b; % parameter of ray equation where ray intersects triangle plane
if r>0, % ray goes away from triangle
    intersect=0;
    return;
else
    ipoint=ray(1,:)-r*dir; % intersection point with triangle plane
end;

% check whether intersection point is inside triangle

uu=sum(u.*u);
uv=sum(u.*v);
vv=sum(v.*v);
w=ipoint-triangle(1,:);
wu=sum(w.*u);
wv=sum(w.*v);
D=uv^2-uu*vv;

% compute parametric coordinates
s=(uv*wv-vv*wu)/D;
if s<0 || s>1,
    intersect=0;
    return;
end;
t=(uv*wu-uu*wv)/D;
if t<0 || (s+t)>1,
    intersect=0;
    return;
end;

intersect=1;

if hit_test,
    % now we have to test for vertex or edge hits
    if s==0,
        if t==0,
            hit=1;
        elseif t==1
            hit=3;
        else
            hit=[1,3];
        end;
    else
        if t==0,
            if s==1,
                hit=2;
            else
                hit=[1,2];
            end;
        elseif s+t==1,
            hit=[2,3];
        end;
    end;
end;
