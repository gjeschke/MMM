function coor1=affine_trafo_coor(coor0,matrices)
% function coor1=affine_trafo_coor(coor0,matrices)
%
% Performs 3D affine transformation(s) on a [m,3] coordinate array
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% use functions affine_trafo_point or affine_trafo_vector 
% for transformations of three-element coordinate vectors
% use function affine to obtain affine transformation matrices
%
% no test for proper input format (vector & matrix dimensions)
%
% coor0     [m,3] Cartesian coordinate array
% matrices  can be a single 4x4 affine transformation matrix of the form 
%           obtained with function affine or a one-dimensional cell array
%           of such matrices, in the latter case the transformations are
%           performed subsequently from left to right
%
% G. Jeschke, 2013

[mm,~]=size(coor0);

if iscell(matrices),
    matrix=eye(4);
    num=length(matrices);
    for k=1:num,
        matrix=matrices{k}*matrix;
    end;
else
    matrix=matrices;
end;

xyz=[coor0 ones(mm,1)];
xyz=matrix*xyz';
xyz=xyz';
coor1=xyz(:,1:3);
            
