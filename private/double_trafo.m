function [coorb,transmat] = double_trafo(coor,orig,Rp,iorig,iRp)

[ma,~] = size(coor);
transmat1 = eye(4);
transmat1(1:3,4) = -iorig';
transmat2 = eye(4);
transmat2(1:3,1:3) = iRp;
transmat3 = eye(4);
transmat3(1:3,1:3) = Rp';
transmat4 = eye(4);
transmat4(1:3,4) = orig';
transmat = transmat4*transmat3*transmat2*transmat1;
coor0 = [coor ones(ma,1)];
coor1=transmat*coor0';
coor1=coor1';
coorb=coor1(:,1:3);

