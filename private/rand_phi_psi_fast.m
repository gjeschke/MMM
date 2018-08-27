function [rphi,rpsi]=rand_phi_psi_fast(phi,psi,reduced)
% Pair of random backbone dihedral angles from a 
% given a population distribution of a Ramachandran plot and corresponding
% phi and psi axes, the distribution of the random pairs conforms to the
% given population distribution
%
% phi   axis of backbone dihedrals phi
% psi   axis of backbone dihedrals psi
% pop   population matrix with dimensions of phi, psi axes
%
% rphi  random backbone angle phi
% rpsi  random backbone angle psi
%
% G. Jeschke, 2013

[~,poi]=min(abs(reduced-rand));
rphi=phi(poi);
rpsi=psi(poi);
