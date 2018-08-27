function coor1=affine_trafo(coor0,matrices)
% function coor1=affine_trafo(coor0,matrices)
%
% Performs 3D affine transformation(s) on a four-element coordinate vector
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% use functions affine_trafo_point or affine_trafo_vector 
% for transformations of three-element coordinate vectors
% use function affine to obtain affine transformation matrices
%
% no test for proper input format (vector & matrix dimensions)
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

if iscell(matrices),
    matrix=eye(4);
    num=length(matrices);
    for k=1:num,
        matrix=matrices{k}*matrix;
    end;
else
    matrix=matrices;
end;
        
coor1=matrix*coor0;

if m<n,
    coor1=coor1';
end;