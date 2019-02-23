function coor2 = trans_rot_trans(coor1,trans1,rot,trans2)
% coor2 = trans_rot_trans(coor1,trans1,rot,trans2)
% 
% Performs a translation, rotation, translation sequence on coordinates
% this function exists to avoid the large translations and rotations
% between very similar coordinate sets that may appear in an affine
% transformation matrix approach
%
% coor1     [m,3] coordinate array
% trans1    translation vector to frame where the rotation is performed
% rot       rotation matrix
% trans2    translation vector to final frame
%
% G. Jeschke, 6.2.2019

[m,~] = size(coor1);
tmat1 = repmat(trans1,m,1);
coor2 = coor1 + tmat1;
coor2 = rot*coor2';
coor2 = coor2';
tmat2 = repmat(trans2,m,1);
coor2 = coor2 + tmat2;