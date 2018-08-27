function coor1=affine_trafo_point(coor0,matrices)
% function coor1=affine_trafo_point(coor0,matrices)
%
% Performs 3D affine transformation(s) on the three-element coordinates of
% a point
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% use function affine_trafo_vector for transformations of vectors
% use function affine to obtain affine transformation matrices
%
% no test for proper input format (vector & matrix dimensions)
% calls function affine_trafo
%
% coor0     original coordinates, can be row or column vector
% matrices  can be a single 4x4 affine transformation matrix of the form 
%           obtained with function affine or a one-dimensional cell array
%           of such matrices, in the latter case the transformations are
%           performed subsequently from left to right
%
% G. Jeschke, 2009

[m,n]=size(coor0);

if m<n,
    coor0=coor0';
end;

coor0=[coor0;1];

coor1=affine_trafo(coor0,matrices);
coor1=coor1(1:3);

if m<n,
    coor1=coor1';
end;