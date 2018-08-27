function newcoor=rot_trans_vectorized(coor,rotmat,t)
%
% newcoor=rot_trans(coor,rotmat,t),
%
% Rotation and subsequent translation of a coordinate array
% 
% coor      x,y,z coordinate array
% rotmat    rotation matrix, use get_rotmat to create from Euler angles
% t         vector for coordinate translation (xt,yt,zt)
%

[m,~] = size(coor);
coor=coor';
newcoor = rotmat*coor;
trans = t(ones(m, 1), :);
newcoor = newcoor + trans';
newcoor=newcoor';

% [m,~] = size(coor);
% coor=coor';
% newcoor = rotmat*coor;
% t=t';
% trans = t(:,ones(1,m));
% newcoor = newcoor + trans;
% newcoor=newcoor';