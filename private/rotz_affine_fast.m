function transmat=rotz_affine_fast(angle,transmat)
% function transmat=rotz_affine_fast(angle,transmat)
%
% Creates a 4x4 matrix for an affine coordinate transformation in 3D space
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% use functions affine_trafo, affine_trafo_point, and affine_trafo_vector 
% for performing the actual transformations
%
% fast version for Monte Carlo trials in helix bundle assembler, created
% from affine.m
%
% angle   rotation angle phi (radians) for rotation about z axis
% transmat         must be an eye(4) matrix, providing it obviates
%                  expensive initialization
%
% G. Jeschke, 2012

cfi=cos(angle);
sfi=sin(angle);
transmat(1,1)=cfi;
transmat(1,2)=-sfi;
transmat(2,1)=sfi;
transmat(2,2)=cfi;
       