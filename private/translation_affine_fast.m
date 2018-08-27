function transmat=translation_affine_fast(transvec,transmat)
% function transmat=translation_affine_fast(transvec,transmat)
%
% Creates a 4x4 matrix for an affine coordinate transformation in 3D space
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% use functions affine_trafo, affine_trafo_point, and affine_trafo_vector 
% for performing the actual transformations
%
% fast version for Monte Carlo trials in helix bundle assembler, created
% from affine.m
% transvec         coordinates [x,y,z] of new origin
% transmat         must be an eye(4) matrix, providing it obviates
%                  expensive initialization
%
% G. Jeschke, 2012
%

transmat(1:3,4)=transvec';
       