function test_rots

a=56;
b=40;
g=137;

euler2 = pi*[a,b,180]/180;
trans2 = [-10,30,47];

% coor=[2,3,4;5,6,7;8,9,0;11,12,13];

aaa=[2 3 4];
N=100000;
coorN=aaa(ones(N,1),:);
coor=coorN;


Rp2=get_rotmat(euler2);

tic,
% for k=1:100000,
    coor1 = rot_trans(coor,Rp2,trans2);
% end;
toc,
tic,
% for k=1:100000,
    coor2 = rot_trans_vectorized(coor,Rp2,trans2);
% end;
toc,

[m,~] = size(coor);
coor=coor';
trans2 = trans2';
transx = trans2(:,ones(1, m));
tic,
% for k = 1:100000,
    coor3 = rot_trans_fast(coor,Rp2,transx);
% end;
toc,
coor3 = coor3';

keyboard

