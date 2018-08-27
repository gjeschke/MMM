function [rms,coor2b,trans,rot]=fit_structures(coor1,coor2)
%
% function [rms,coor2b,trans,rot]=fit_structures(coor1,coor2)
%
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
% rms     root mean square deviation of both coordinate sets
% coor2b  coordinate set 2 transformed to give the best fit with coordinate
%         set 1
% trans   translation vector for centroid of coordinate set 2
% rot     rotation matrix for coordinate set 2
%
% trans and rot are supplied so that after a backbone fit the whole
% structure 2 can be transformed with function rot_trans_structures
%
% (c) G. Jeschke, 2007

rms=-1e6; % values for failed fit
coor2b=coor2;
trans=zeros(1,3);
rot=eye(3);

[m1,n1]=size(coor1);
[m2,n2]=size(coor2);

% Check if input coordinate arrays have the correct dimension and are
% consistent
if m1~=m2 || n1~=3 || n2 ~=3,
    return;
end;

cent1=sum(coor1)/m1;
cent2=sum(coor2)/m2;
shift1=kron(ones(m1,1),cent1);
shift2=kron(ones(m2,1),cent2);
coor1b=coor1-shift1;
coor2b=coor2-shift2;
trans=cent1-cent2;
bas=0;
for k=1:m1,
    bas=bas+sum(coor1b(k,:).*coor1b(k,:))+sum(coor2b(k,:).*coor2b(k,:));
end;
bas=bas/(2*m1);
    

U=coor1b'*coor2b/m1;

[h,D0,v]=svd(U);

UUp=U*U';

D=sort(sqrt(eig(UUp)),1,'descend'); 


if det(U)>0,
    rms=bas-sum(D);
    rot=v*h';
else
    rms=bas-D(1)-D(2)+D(3);
    rot=zeros(3,3);
    for ii=1:3,
        for jj=1:3,
            for k=1:3,
                signum=1;
                if k==3, signum=-1; end;
                rot(ii,jj)=rot(ii,jj)+signum*v(ii,k)*h(jj,k);
            end;
        end;
    end;
end;

coor2c=coor2b;
for k=1:m1,
    coor2b(k,:)=coor2c(k,:)*rot+cent1;
end;

rms=sqrt(2*rms);
