% Spline interpolation for given four control points in 3D space

function coor=spline_interp(fourpoints,N)
% Input:
% fourpoints  : four control points of bezier curve
%               as 4*D array, D is number of dimensions
% N           : number of interpolation points 


[m,n]=size(fourpoints);
x0=fourpoints(:,1);
y0=fourpoints(:,1);
z0=fourpoints(:,1);

t=linspace(0,1,N)';
t0=linspace(0,1,m)';

coor=[interp1(t0,x0,t,'spline'),interp1(t0,y0,t,'spline'),interp1(t0,z0,t,'spline')];
