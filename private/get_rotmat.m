function Rp=get_rotmat(euler)
%
% function Rp=get_rotmat(euler),
%
% Compute rotation matrix
%
% euler     Euler angles in radians
% Rp        rotation matrix
%

ca=cos(euler(1));
sa=sin(euler(1));
cb=cos(euler(2));
sb=sin(euler(2));
cg=cos(euler(3));
sg=sin(euler(3));

Rp=[ca*cb*cg-sa*sg,sa*cb*cg+ca*sg,-sb*cg;-ca*cb*sg-sa*cg,-sa*cb*sg+ca*cg,sb*sg;ca*sb,sa*sb,cb];

