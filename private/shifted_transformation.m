function coor=shifted_transformation(coor0,transmat,origin)
% function coor=shifted_transformation(coor0,transmat,origin)
%   Affine coordinate transformation with respect to coordinate frame with
%   shifted origin
%
% coor0     input coordinates
% transmat  affine transformation matrix 4x4
% origin    origin of the frame, to which transmat relates
%
% coor      output coordinates (in original frame)
%
% G. Jeschke, 2010

[m,n]=size(coor0);
ormat=repmat(origin,m,1);
coor0b=coor0-ormat;
coor1=affine_coor_set(coor0b,transmat);
coor=coor1+ormat;