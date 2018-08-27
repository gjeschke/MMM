function newcoor=rot_trans(coor,rotmat,t)
%
% newcoor=rot_trans(coor,rotmat,t),
%
% Rotation and subsequent translation of a coordinate array
% 
% coor      x,y,z coordinate array
% rotmat    rotation matrix, use get_rotmat to create from Euler angles
% t         vector for coordinate translation (xt,yt,zt)
%

coor=coor';
[m,n]=size(coor);
newcoor=coor;
for k=1:n,
    newcoor(:,k)=rotmat*coor(:,k)+t';
end;
newcoor=newcoor';