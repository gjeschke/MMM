function out = get_energyLJ_charmm27(coor1,coor2,run_opt)

% out = get_energyLJ_charmm27(coor1,coor2,run_opt)

% accepts coor1 and coor2 - Nx4 arrays for the molecule 1 and molecule 2
% (:,1) - atom types (6 for carbon etc); (:,2:4) - xyz coordinates for each atom

% run_opt - various run options
% orig_ind - optinal. Here original index (better - numbers) of the original coordinate array

global LJ

conv_factor=(4.185*1000); % conversion to SI, - if working per mol

%--------------------------------------------------------------------------
used_potential='charmm27 parametrized Lennard-Jones';    % potential type used
%--------------------------------------------------------------------------
energy_stats=run_opt.energy_stats;  % choices: 0 or 1
% 0 - no statistics made (should be selected first)
% 1 - statistics for every atom pairs (additional information, but a huge file, time consuming)
forgive=run_opt.forgive; % get softening factor for the LJ potential
%--------------------------------------------------------------------------
[m1,n1] = size(coor1); % get sizes of the coordinates arrays
[m2,n2] = size(coor2);

a = coor1(:,2:4);
b = coor2(:,2:4);
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));

%{
% initiate statistics matrix if needed
if energy_stats==1
  stats=zeros(m1*m2,5);
  ind_coor1=run_opt.ind_coor1;
  ind_coor2=run_opt.ind_coor2;
  s=0; % initilize counter
end
%}
stats = [];

Rmin2_i0 = LJ.Rmin2(coor1(:,1));
eps_i = LJ.eps(coor1(:,1))*conv_factor;
Rmin2_j0 = LJ.Rmin2(coor2(:,1));
eps_j = LJ.eps(coor2(:,1))*conv_factor;

% [r_vdw_i,q_i,Rmin2_i0,eps_i] = get_simpleCH27forLJstst(coor1(:,1));
% [r_vdw_j,q_j,Rmin2_j0,eps_j] = get_simpleCH27forLJstst(coor2(:,1));
eps = sqrt(eps_i(:)*eps_j);
Rmin2_i = Rmin2_i0*forgive; % soften LJ with the forgive factor
Rmin2_j = Rmin2_j0*forgive; % soften LJ with the forgive factor
Rmin = repmat(Rmin2_i(:),1,m2) + repmat(Rmin2_j,m1,1);
            % ATTENTION! Rmin2_i - is actually Rmin_i/2
q = (Rmin./pair_dist).^2;
q = q.*q.*q;
Energy = eps.*(q.^2-2*q);
Energy(pair_dist>10) = 0;
TotalEnergy = sum(sum(Energy));

  %{
    if energy_stats==1
      s=s+1;
      stats(s,1)=ind_coor1(i1);
      stats(s,2)=ind_coor2(i2);
      stats(s,3)=Energy_;
      stats(s,4)=pair_dist;
      stats(s,5)=s;
    end
  %}

out.energy = TotalEnergy;
out.used = used_potential;
if energy_stats==1
  out.clash_stats=stats;
end
