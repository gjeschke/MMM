function Rg = ensemble_radius_of_gyration(allcoor,N,p)
% Rg = ensemble_radius_of_gyration(coor,N)
% 
% Radius of gyration of a set of atom coordinates
%
% allcoor  (N*l,3) matrix of Cartesian atom coordinates
% N         number of ensemble members
% p         population vector (1,N), if empty or missing, all populations
%           are 1/N
% Rg        radius of gyration
%
% G. Jeschke, 23.12.2019

if ~exist('p','var') || isempty(p)
    p = ones(1,N)/N;
end

[m,~] = size(allcoor);
l = m/N;
if abs(round(l)-l) > 5*eps
    warning('ensemble_radius_of_gyration: matrix size does not fit ensemble size');
    return
end
Rg2 = 0;
% Rgvec = zeros(1,N);
for k = 1:N
    coor = allcoor(1+(k-1)*l:k*l,:);
    rc = mean(coor); % center coordinate
    rel_coor = coor - repmat(rc,l,1);
    % Rgvec(k) = sqrt(sum(sum(rel_coor.^2))/l);
    Rg2 = Rg2 + p(k)*sum(sum(rel_coor.^2))/l;
end

% Rg = sum(p.*Rgvec);
Rg = sqrt(Rg2);
