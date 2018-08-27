function out = get_energyLJ_UFF_attract(coor1,coor2,run_opt)

% out = get_energyLJ_UFF(coor1,coor2,run_opt)

% accepts coor1 and coor2 - Nx4 arrays for the molecule 1 and molecule 2
% (:,1) - atom types (6 for carbon etc); (:,2:4) - xyz coordinates for each atom

% run_opt - various run options
% orig_ind - optinal. Here original index (better - numbers) of the original coordinate array

global uff

enhance_attraction = 1; %3.5;


conv_factor=(4.185*1000); % conversion to SI, - if working per mol

%--------------------------------------------------------------------------
used_potential='UFF Towhee';    % potential type used
%--------------------------------------------------------------------------
% 0 - no statistics made (should be selected first)
% 1 - statistics for every atom pairs (additional information, but a huge file, time consuming)
forgive=run_opt.forgive; % get softening factor for the LJ potential
if length(forgive) > 1
    enhance_attraction = forgive(2);
    forgive = forgive(1);
end
%--------------------------------------------------------------------------

[m1,~] = size(coor1); % get sizes of the coordinates arrays
[m2,~] = size(coor2);

a = coor1(:,2:4);
b = coor2(:,2:4);
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));


Rmin2_i0 = uff.LJ_r(coor1(:,1));
eps_i = uff.LJ_D(coor1(:,1))*conv_factor;
Rmin2_j0 = uff.LJ_r(coor2(:,1));
eps_j = uff.LJ_D(coor2(:,1))*conv_factor;

eps = sqrt(eps_i(:)*eps_j);
Rmin2_i = Rmin2_i0*forgive; % soften LJ with the forgive factor
Rmin2_j = Rmin2_j0*forgive; % soften LJ with the forgive factor
Rmin = sqrt(Rmin2_i(:)*Rmin2_j); % UFF has geometric average combination rule
% Rmin = sqrt(repmat(Rmin2_i(:),1,m2).*repmat(Rmin2_j,m1,1)); % UFF has geometric average combination rule
            % ATTENTION! Rmin2_i - is actually Rmin_i/2
q = (Rmin./pair_dist).^2;
q = q.*q.*q;
Energy = eps.*(q.^2-2*q*enhance_attraction);
Energy(pair_dist>10) = 0;
TotalEnergy = sum(sum(Energy));

out.energy = TotalEnergy;
out.used = used_potential;
