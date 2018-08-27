function rms=rms_sites_homdimer(constraints,sites,v) 

% Compute r.m.s. deviation of site-to-site distances from given constraints
% version for a heterodimer dimer

% v             parameter vector, angles alpha v(1), beta v(2) and gama v(3)
%               (in degrees) Euler rotation of the (second) monomer molecule w.r.t. the
%               fixed (first) monomer molecule; 
%               x v(3), y v(4) and z v(3) (in Angstroem) determine translation of 
%               second molecule w.r.t. first molecule
%
% sites         coordinates of sites
% constraints   structure with constraints

% take care about homooligomer multiplicity
if isfield(constraints,'multiplicity')
    multi=constraints.multiplicity;
else
    multi = 2;
end

ph=v(1);
th=v(2);
x=v(3);
y=v(4);

sites0=sites.NOsite1;
pairs_exp=constraints.pairs_exp;
sig=pairs_exp(:,4);
[np,~]=size(pairs_exp);
    
% keyboard
%- up to 20% speed-up------------------------------------------------------
% define orientation of the c2 dimer axis wrt molecule one:
t0=[0,0,0]; % no translation of first molecule
eulang0=[ph,th,0]; % rotation of first molecule
eulang=[ph,th,360/multi*pi/180];
Rp1=get_rotmat(eulang0);
t=[x,y,0]; % translation of second moelcule (already in the "dimer frame")
% eulang=[360/multi,0,0]; % rotation of the second molecule (basically, a reflection about the c2 axis)
Rp2=get_rotmat(eulang);

sites1=rot_trans_vectorized(sites0,Rp1,t0); % protein coordinates of the first molecule in a dimer
sites2=rot_trans_vectorized(sites0,Rp2,t);

% % keep it for future cases for speed up (here it does not bring much)
% [m,~]=size(sites0);
% sites0=sites0';
% t0=t0';
% t0x = t0(:,ones(1, m));
% t=t';
% tx = t(:,ones(1, m));
% sites1=rot_trans_fast(sites0,Rp1,t0x); % protein coordinates of the first molecule in a dimer
% sites1=sites1';
% sites2=rot_trans_fast(sites0,Rp2,tx);
% sites2=sites2';

% sites1=rot_trans(sites0,Rp1,t0); % protein coordinates of the first molecule in a dimer
% sites2=rot_trans(sites0,Rp2,t);


% keyboard
diff=sites1-sites2;
dist=sqrt(diff(:,1).*diff(:,1)+diff(:,2).*diff(:,2)+diff(:,3).*diff(:,3));
% ### dist = sqrt(sum(diff.^2,2));
var=((dist-pairs_exp(:,3))./sig).*((dist-pairs_exp(:,3))./sig);
% var=((dist-pairs_exp(:,3))).*((dist-pairs_exp(:,3)));
rms=sqrt(sum(var/np));
%--------------------------------------------------------------------------


