function [rms,coor2b,transmat]=rmsd_superimpose0(coor1,coor2)
%
% function [rms,coor2b,transmat]=rmsd_superimpose(coor1,coor2)
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
% transmat   affine transformation matrix for coordinate set 2
%
% transmat is supplied so that after a fit of a substructure other parts of
% structure 2 can be transformed with transform_structure (or simple
% multiplication of a coordinate array with transmat)
%
% (c) G. Jeschke, 2009

rms=-1e6; % values for failed fit
coor1b=coor1;
coor2b=coor2;
transmat=eye(4);

[m1,n1]=size(coor1);
[m2,n2]=size(coor2);

% Check if input coordinate arrays have the correct dimension and are
% consistent
if m1~=m2 || n1~=3 || n2 ~=3,
    return;
end;

cent1=sum(coor1)/m1;
cent2=sum(coor2)/m2;
% shift1=repmat(cent1,m1,1);
% shift2=repmat(cent2,m2,1);
% coor1b=coor1-shift1;
% coor2b=coor2-shift2;
bas=0;
for k=1:m1,
    coor1b(k,:)=coor1(k,:)-cent1;
    coor2b(k,:)=coor2(k,:)-cent2;
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

transmat1=affine('translation',-cent2); % shift to centroid of structure to be fitted
transmat2(1:3,1:3)=rot'; % rotate to new model frame
transmat2(4,4)=1; 
transmat3=affine('translation',cent1); % shift to target centroid
transmat=transmat3*transmat2*transmat1;
