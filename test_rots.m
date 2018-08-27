function test_rots

a=0;
b=40;
g=0;

euler1 = pi*[a,b,0]/180;
euler2 = pi*[a,b,180]/180;
trans1 = [0,0,0];
trans2 = [-10,30,0];

coor=[2,3,4];

Rp1=get_rotmat(euler1);
Rp2=get_rotmat(euler2);

coor1 = rot_trans(coor,Rp1,trans1),
coor2 = rot_trans(coor,Rp2,trans2),
dist1 = norm(coor1-coor2),

a=90;
b=120;
g=0;

euler1 = pi*[a,b,0]/180;
euler2 = pi*[a,b,180]/180;
trans1 = [0,0,0];
trans2 = [-30,10,0];

coor=[2,3,4];

Rp1=get_rotmat(euler1);
Rp2=get_rotmat(euler2);

coor3 = rot_trans(coor,Rp1,trans1),
coor4 = rot_trans(coor,Rp2,trans2)
dist1 = norm(coor3-coor4),
