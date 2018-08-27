function out = get_energyLJ_UFF_Bavg(coor1,coor2,run_opt,Bfac1,Bfac2)

% out = get_energyLJ_UFF_Bavg(coor1,coor2,run_opt,Bfac1,Bfac2)

% accepts coor1 and coor2 - Nx4 arrays for the molecule 1 and molecule 2
% (:,1) - atom types (6 for carbon etc); (:,2:4) - xyz coordinates for each atom

% run_opt - various run options
% orig_ind - optinal. Here original index (better - numbers) of the original coordinate array
% Bfac1     B factors for molecule 1
% Bfac2     B factors for molecule 2

global uff
global avg_LJ_pot

reduce_B = 0.000001;
conv_factor=(4.185*1000); % conversion to SI, - if working per mol
gas_un=8.314472;

%--------------------------------------------------------------------------
used_potential='UFF Towhee';    % potential type used
%--------------------------------------------------------------------------
% 0 - no statistics made (should be selected first)
% 1 - statistics for every atom pairs (additional information, but a huge file, time consuming)
forgive=run_opt.forgive; % get softening factor for the LJ potential
%--------------------------------------------------------------------------
[m1,~] = size(coor1); % get sizes of the coordinates arrays
[m2,~] = size(coor2);

a = coor1(:,2:4);
b = coor2(:,2:4);
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));

Bfac_a = repcolvector(Bfac1,m2);
Bfac_b = repcolvector(Bfac2,m1).';
Bfac = reduce_B*(Bfac_a+Bfac_b)/2;
% Bfac = sqrt(Bfac_a.*Bfac_b);

Rmin2_i0 = uff.LJ_r(coor1(:,1));
eps_i = uff.LJ_D(coor1(:,1)); % no conversion, well depth is renormalized
Rmin2_j0 = uff.LJ_r(coor2(:,1));
eps_j = uff.LJ_D(coor2(:,1));

eps = sqrt(eps_i(:)*eps_j);
Rmin2_i = Rmin2_i0*forgive; % soften LJ with the forgive factor
Rmin2_j = Rmin2_j0*forgive; % soften LJ with the forgive factor
Rmin = sqrt(repmat(Rmin2_i(:),1,m2).*repmat(Rmin2_j,m1,1)); % UFF has geometric average combination rule
            % ATTENTION! Rmin2_i - is actually Rmin_i/2

mm = m1*m2;
r = reshape(pair_dist,1,mm)./reshape(Rmin,1,mm);
D = reshape(eps,1,mm);
std = sqrt(reshape(Bfac,1,mm)/(24*pi^2));
            
poti = interp2(avg_LJ_pot.rm,avg_LJ_pot.sigrm,avg_LJ_pot.pot,r,std,'linear',0);
Energy = gas_un*avg_LJ_pot.T*(D.*poti)/avg_LJ_pot.D;
Energy(pair_dist>10) = 0;
TotalEnergy = sum(Energy);

out.energy = TotalEnergy;
out.used = used_potential;
