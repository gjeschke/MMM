function rotmat=opt_rot(coor1,coor2)
%
% function rotmat=opt_rot(coor1,coor2)
%
% returns the rotation matrix that superimposes coordinate set 2 onto 
% coordinate set 1 with least r.m.s.d.
% Least square fitting of two sets of atomic position vectors according to
% A. D. McLachlan, J. Mol. Biol. 1979, 128, 49-79 "Gene Duplications in the
% Structural Evolution of Chymotrypsin"
% 
% equal weights wi=1/m are assumed for all atomic positions
%
% coor1   master coordinate set, must be an (m x 3) array
% coor2   coordinate set to be fitted to coor1 by translation and rotation,
%         must be an (m x 3) array with the same m, otherwise rms=-1e6 is
%         output
%
% rotmat  rotation matrix for coordinate set 2, empty, if input is
%         inconsistent
%
% (c) G. Jeschke, 2010

coor1b=coor1;
coor2b=coor2;

[m1,n1]=size(coor1);
[m2,n2]=size(coor2);

% Check if input coordinate arrays have the correct dimension and are
% consistent
if m1~=m2 || n1~=3 || n2 ~=3,
    rotmat=[];
    return;
end;

cent1=sum(coor1)/m1;
cent2=sum(coor2)/m2;
for k=1:m1,
    coor1b(k,:)=coor1(k,:)-cent1;
    coor2b(k,:)=coor2(k,:)-cent2;
end;

U=coor1b'*coor2b/m1;

[h,D0,v]=svd(U);

if det(U)>0,
    rotmat=(v*h')';
else
    rotmat=zeros(3,3);
    for ii=1:3,
        for jj=1:3,
            for k=1:3,
                signum=1;
                if k==3, signum=-1; end;
                rotmat(jj,ii)=rotmat(jj,ii)+signum*v(ii,k)*h(jj,k);
            end;
        end;
    end;
end;

