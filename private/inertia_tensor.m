function inertia=inertia_tensor(xyz,cflag,masses)
% function inertia=inertia_tensor(xyz,cflag,masses)
%
% inertia tensor of an object
%
% xyz       Cartesian coordinates of the mass points
% cflag     optional flag that dermines if center coordinates are
%           subtracted before computation, 0 no subtraction, 1 subtraction,
%           defaults to 1
% masses    optional vector of masses, defaults to a unity vector
%
% G. Jeschke, 2009

[n,m]=size(xyz);

if nargin<2,
    cflag=1;
end;

if cflag,
    xyzc=mean(xyz);
    for k=1:n, xyz(k,:)=xyz(k,:)-xyzc; end;
end;

if nargin<3,
    masses=ones(n,1);
end;

% compute inertia tensor

inertia=zeros(3);

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);

inertia(1,1)=sum(masses.*y.^2)+sum(masses.*z.^2);
inertia(1,2)=-sum(masses.*x.*y);
inertia(1,3)=-sum(masses.*x.*z);
inertia(2,1)=-sum(masses.*y.*x);
inertia(2,2)=sum(masses.*x.^2)+sum(masses.*z.^2);
inertia(2,3)=-sum(masses.*y.*z);
inertia(3,1)=-sum(masses.*z.*x);
inertia(3,2)=-sum(masses.*z.*y);
inertia(3,3)=sum(masses.*x.^2)+sum(masses.*y.^2);
