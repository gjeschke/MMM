function [pst,center] = point_spread_tensor(xyz,weights)
% function [pst,center] = point_spread_tensor(xyz,weights)
%
% points spread tensor of an object
%
% xyz       Cartesian coordinates of the mass points
% weigts    optional weights, defaults to a equal weights 
%
% G. Jeschke, 2017

[n,~]=size(xyz);

if ~exist('weights','var'),
    weights = ones(n,1);
end;
weights = weights/sum(weights);

center = weights'*xyz;
xyz = xyz - repmat(center,n,1);


% compute pst tensor

pst=zeros(3);

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);

pst(1,1)=sum(weights.*x.^2);
pst(1,2)=sum(weights.*x.*y);
pst(1,3)=sum(weights.*x.*z);
pst(2,1)=sum(weights.*y.*x);
pst(2,2)=sum(weights.*y.^2);
pst(2,3)=sum(weights.*y.*z);
pst(3,1)=sum(weights.*z.*x);
pst(3,2)=sum(weights.*z.*y);
pst(3,3)=sum(weights.*z.^2);
