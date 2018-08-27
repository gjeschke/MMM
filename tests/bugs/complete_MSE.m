function complete_MSE

% N,CA,C,CG,SE,CE
coor307 = [41.237  16.737  53.130;...
    42.371  16.835  52.177;...
    43.568  15.856  52.378;...
    41.824  19.231  51.663;...
    41.936  20.998  52.566;...
    44.013  21.333  52.168];

% rms 0.3554 Å
coor52 = [14.766  55.032  33.921;...
    14.597  55.492  35.245;...
    15.731  55.146  36.118;...
    12.120  55.773  35.361;...
    10.395  55.367  36.241;
    10.745  54.451  37.561];


[rms,coor52b,transmat] = rmsd_superimpose(coor307,coor52);

coor52_N = [14.766  55.032  33.921 1];
coor52_O = [16.008  55.869  37.071 1];
coor52_CB = [13.319  54.939  35.819 1];

disp('Aber hallo!');


function [rms,coor2b,transmat]=rmsd_superimpose(coor1,coor2,weights)
%
% function [rms,coor2b,transmat]=rmsd_superimpose(coor1,coor2,weights)
%
% Least square fitting of two sets of atomic position vectors according to
% A. D. McLachlan, J. Mol. Biol. 1979, 128, 49-79 "Gene Duplications in the
% Structural Evolution of Chymotrypsin"
% 
% equal weights wi=1/m are assumed for all atomic positions, unless the
% third argument is given
%
% coor1   master coordinate set, must be an (m x 3) array
% coor2   coordinate set to be fitted to coor1 by translation and rotation,
%         must be an (m x 3) array with the same m, otherwise rms=-1e6 is
%         output
% weights vector of length m with weights for the atomic positions, can be
%         column or row vector, defaults to uniform weights
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

if nargin>2,
    [mw,nw]=size(weights);
    if nw>mw,
        mw=nw;
        weights=weights.';
    end;
else
    weights=ones(m1,1);
    mw=m1;
end;

% Check if input coordinate arrays have the correct dimension and are
% consistent
if m1~=m2 || n1~=3 || n2 ~=3 || mw~=m1,
    return;
end;

repweights=repmat(sqrt(weights),1,3);

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
    

U=(coor1b.*repweights)'*(coor2b.*repweights)/sum(weights);

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

diff=coor2b-coor1;
rms=sqrt(sum(sum(diff.^2))/m1);

transmat1=affine('translation',-cent2); % shift to centroid of structure to be fitted
transmat2(1:3,1:3)=rot'; % rotate to new model frame
transmat2(4,4)=1; 
transmat3=affine('translation',cent1); % shift to target centroid
transmat=transmat3*transmat2*transmat1;

function transmat=affine(mode,param1,param2,param3)
% function transmat=affine(mode,param1,param2,param3)
%
% Creates a 4x4 matrix for an affine coordinate transformation in 3D space
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% use functions affine_trafo, affine_trafo_point, and affine_trafo_vector 
% for performing the actual transformations
% an empty matrix is returned when a non-existing mode is requisted
% function assumes that the correct number of correctly dimensioned
% parameters is provided by the client (caller)
%
% mode      transformation mode, case sensitive string
% param     parameters
%
% Implemented combinations of modes and parameters
% mode              explanation             parameters
% 'identity'        identity transformation (none)
% 'translation'     translation             coordinates [x,y,z] of new
%                                           origin
% 'Euler'           rotation                Euler angles [alpha,beta,gamma]
%                                           for rotations about z, y', and
%                                           z'' axes, in radians
%                                           WARNING: the objects are
%                                           rotated, not the frame, use
%                                           -alpha, -beta, -gamma for the
%                                           frame rotation more usual in
%                                           magnetic resonance
% 'rotx'            rotation about x axis   rotation angle phi (radians)
% 'roty'            rotation about y axis   rotation angle phi (radians)
% 'rotz'            rotation about z axis   rotation angle phi (radians)
% 'rotn'            rotation about an       rotation angle phi (radians)
%                   arbitrary axis          and vector representation
%                                           [nx,ny,nz] of the rotation axis
% 'invert'          inversion at origin     (none)
% 'reflectxy'       reflection at xy plane  (none) 
% 'reflectxz'       reflection at xz plane  (none) 
% 'reflectyz'       reflection at yz plane  (none) 
% 'scale'           scaling                 scaling factors [ax,ay,az] for
%                                           the three coordinates or
%                                           uniform scaling factor a for
%                                           all coordinates
% 'shear'           shearing                vector
%                                           [sxy,sxz,syx,syz,szx,szy] of
%                                           shearing coefficients
% 'screw'           rotation about an       rotation angle phi (radians),
%                   arbitrary axis          vector representation
%                   followed by trans-      [nx,ny,nz] of the rotation axis
%                   lation along the same   and shift dz along the axis
%                   axis
%
% G. Jeschke, 2009

transmat=eye(4); % initialize identity transformation

switch mode
    case 'identity' % needed to prevent otherwise clause for identity mode
    case 'translation'
        transmat(1:3,4)=param1';
    case 'Euler'
        ca=cos(param1(1)); sa=sin(-param1(1));
        cb=cos(param1(2)); sb=sin(-param1(2));
        cg=cos(param1(3)); sg=sin(-param1(3));
        transmat(1,1)=cg*cb*ca-sg*sa;
        transmat(1,2)=cg*cb*sa+sg*ca;
        transmat(1,3)=-cg*sb;
        transmat(2,1)=-sg*cb*ca-cg*sa;
        transmat(2,2)=-sg*cb*sa+cg*ca;
        transmat(2,3)=sg*sb;
        transmat(3,1)=sb*ca;
        transmat(3,2)=sb*sa;
        transmat(3,3)=cb;
    case 'rotx'
        cfi=cos(param1(1));
        sfi=sin(param1(1));
        transmat(2,2)=cfi;
        transmat(2,3)=-sfi;
        transmat(3,2)=sfi;
        transmat(3,3)=cfi;
    case 'roty'
        cfi=cos(param1(1));
        sfi=sin(param1(1));
        transmat(1,1)=cfi;
        transmat(1,3)=sfi;
        transmat(3,1)=-sfi;
        transmat(3,3)=cfi;
    case 'rotz'
        cfi=cos(param1(1));
        sfi=sin(param1(1));
        transmat(1,1)=cfi;
        transmat(1,2)=-sfi;
        transmat(2,1)=sfi;
        transmat(2,2)=cfi;
    case 'rotn'
        c=cos(param1(1));
        s=sin(param1(1));
        t=1-c;
        n=param2/norm(param2);
        nx=n(1);
        ny=n(2);
        nz=n(3);
        transmat(1,1)=t*nx^2+c;
        transmat(1,2)=t*nx*ny-s*nz;
        transmat(1,3)=t*nx*nz+s*ny;
        transmat(2,1)=t*nx*ny+s*nz;
        transmat(2,2)=t*ny^2+c;
        transmat(2,3)=t*ny*nz-s*nx;
        transmat(3,1)=t*nx*nz-s*ny;
        transmat(3,2)=t*ny*nz+s*nx;
        transmat(3,3)=t*nz^2+c;
    case 'invert'
        transmat(1,1)=-1;
        transmat(2,2)=-1;
        transmat(3,3)=-1;
    case 'reflectxy'
        transmat(3,3)=-1;
    case 'reflectxz'
        transmat(2,2)=-1;
    case 'reflectyz'
        transmat(1,1)=-1;
    case 'scale'
        if length(param1)==1, param1=[1,1,1]*param1; end;
        transmat(1,1)=param1(1);
        transmat(2,2)=param1(2);
        transmat(3,3)=param1(3);
    case 'shear'
        transmat(2,1)=param1(1);
        transmat(3,1)=param1(2);
        transmat(1,2)=param1(3);
        transmat(3,2)=param1(4);
        transmat(1,3)=param1(5);
        transmat(2,3)=param1(6);
    case 'screw'
        c=cos(param1(1));
        s=sin(param1(1));
        t=1-c;
        n=param2/norm(param2);
        nx=n(1);
        ny=n(2);
        nz=n(3);
        transmat(1,1)=t*nx^2+c;
        transmat(1,2)=t*nx*ny-s*nz;
        transmat(1,3)=t*nx*nz+s*ny;
        transmat(2,1)=t*nx*ny+s*nz;
        transmat(2,2)=t*ny^2+c;
        transmat(2,3)=t*ny*nz-s*nx;
        transmat(3,1)=t*nx*nz-s*ny;
        transmat(3,2)=t*ny*nz+s*nx;
        transmat(3,3)=t*nz^2+c;
        transmat(1,4)=param3*param2(1);
        transmat(2,4)=param3*param2(2);
        transmat(3,4)=param3*param2(3);
    otherwise
        transmat=[];
end;
        
       