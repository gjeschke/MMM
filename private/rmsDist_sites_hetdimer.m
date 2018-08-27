function [rms,dist]=rmsDist_sites_hetdimer(constraints,sites,v) 

% Compute r.m.s. deviation of site-to-site distances from given constraints
% version for a heterodimer dimer and returns also distances corresponding
% to the current transformation

% v             parameter vector, angles alpha v(1), beta v(2) and gama v(3)
%               (in degrees) Euler rotation of the (second) monomer molecule w.r.t. the
%               fixed (first) monomer molecule; 
%               x v(3), y v(4) and z v(3) (in Angstroem) determine translation of 
%               second molecule w.r.t. first molecule
%
% sites         coordinates of sites
% constrainst   structure with constraints

alpha=v(1);
beta=v(2);
gama=v(3);
x=v(4);
y=v(5);
z=v(6);
 
sites1=sites.NOsite1;
sites2=sites.NOsite2;
pairs_exp=constraints.pairs_exp;
sig=pairs_exp(:,4);
[np,~]=size(pairs_exp);
    

eulang=[alpha,beta,gama]; % Euler angles to rotate molecule 2
Rp2=get_rotmat(eulang); % rotation matrix to rotate molecule 2
t=[x,y,z]; % translation vector to translate molecule 2
sites2=rot_trans_vectorized(sites2,Rp2,t);

% keyboard
diff=sites1-sites2;
dist=sqrt(diff(:,1).*diff(:,1)+diff(:,2).*diff(:,2)+diff(:,3).*diff(:,3));
var=((dist-pairs_exp(:,3))./sig).*((dist-pairs_exp(:,3))./sig);
% var=((dist-pairs_exp(:,3))).*((dist-pairs_exp(:,3)));

rms=sqrt(sum(var/np));


