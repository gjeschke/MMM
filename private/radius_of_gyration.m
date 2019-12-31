function Rg = radius_of_gyration(coor)
% Rg = radius_of_gyration(coor)
% 
% Radius of gyration of a set of atom coordinates
%
% coor  (l,3) matrix of Cartesian atom coordinates
% Rg    radius of gyration
%
% G. Jeschke, 23.12.2019

[l,~] = size(coor);
rc = mean(coor); % center coordinate
rel_coor = coor - repmat(rc,l,1);
Rg = sqrt(sum(sum(rel_coor.^2))/l);

