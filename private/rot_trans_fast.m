function newcoor=rot_trans_fast(coor,rotmat,t)
%
% newcoor=rot_trans(coor,rotmat,t),
%
% Rotation and subsequent translation of a coordinate array
% 
% coor      x,y,z coordinate array (3xN)
% rotmat    rotation matrix, use get_rotmat to create from Euler angles
% t         column vector for coordinate translation (xt;yt;zt)
%

newcoor = rotmat*coor+t;
