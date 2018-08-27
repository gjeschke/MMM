function [ origin , tmatrix , rmsd] = opt_transform( coor1, coor2, verbose )
%function [ origin , tmatrix ] = opt_transform( coor1, coor2, verbose )
%   finds optimum affine transformation for superposition of coordinate set
%   coor2 onto coordinate set coor1
%   origin of the frame is the point out of coor2 than is at minimum
%   distance from axis of a screw transform
%   corresponding affine transformation matrix and the rmsd are
%   returned
%
%   empty output values, if size of coordinate arrays does not match or
%   coordinate arrays are empty
%
%   coor1       nx3 set of Cartesian point coordinates
%   coor2       nx3 set of Cartesian point coordinates
%   verbose     optional, if true (1), verbose output is displayed
%
%   origin      origin of the optimum frame to perform transform
%   tmatrix     affine transformation matrix
%   rmsd        root mean square deviation between point sets after
%               transformation of coor2
%
% G. Jeschke, 2010

origin=[];
tmatrix=[];
rmsd=[];

[n,x]=size(coor1);
[n2,x]=size(coor2);

if n==0 || n~=n2, % empty input or size mismatch
    return;
end;

if nargin<3, % suppress verbose output, if not requested
    verbose=false;
end;

orig0=find_superposition(coor1,coor2);
ormat0=repmat(orig0,n,1);
[rmsd,coor2c,transmat]=rmsd_superimpose(coor1-ormat0,coor2-ormat0);
EV=affine2EV(transmat);
screw_axis=EV(1:3);
if verbose,
    fprintf(2,'RMSD of screw transformation is %5.2f\n',rmsd);
    translation_axis=transmat(1:3,4);
    translation_axis=translation_axis'/norm(translation_axis);
    deviat=sum(screw_axis.*translation_axis);
    fprintf(2,'Screw transformation has translation part of %5.4f\n',norm(transmat(1:3,4)));
    fprintf(2,'Cosine deviation of rotation and translation axis is %8.6f\n',1-abs(deviat));
end;
[origin, distance, number]=get_origin(coor2,orig0,screw_axis);
if verbose,
    fprintf(2,'Origin has distance of %5.2f from screw axis\n',distance);
end;
ormat=repmat(origin,n,1);
coor1b=coor1-ormat;
coor2b=coor2-ormat;
[rmsd,coor2c,tmatrix]=rmsd_superimpose(coor1b,coor2b);
if verbose,
    trans=tmatrix(1:3,4);
    fprintf(2,'Affine transformation has translation part of %5.4f\n',norm(trans));
    fprintf(2,'Remaining superposition RMSD is %5.2f\n',rmsd);
end;

function rms=rms_shift_superposition(v,coor1,coor2)

[mc,nc]=size(coor1);
for k=1:mc,
    coor1(k,:)=coor1(k,:)-v;
    coor2(k,:)=coor2(k,:)-v;
end;
[rms0,h102b,transmat]=rmsd_superimpose(coor1,coor2);
tvec=transmat(1:3,4);
rms=norm(tvec);
%fprintf(1,'rmsd: %6.3f\n',rms);

function v=find_superposition(coor1,coor2)

[rms0,h102b,transmat]=rmsd_superimpose(coor1,coor2);
tvec=transmat(1:3,4);
v0=tvec'/2;
v=fminsearch(@rms_shift_superposition,v0,[],coor1,coor2);

function [orig,distance,number]=get_origin(coor,orig0,screw_axis)
% finds point in coordinate set coor that is closest to a screw axis going
% through point orig0 and having a unit direction vector screw_axis
%
% orig      the new origin
% distance  distance from screw axis
% number    number of origin in coordinate set ccor1

p1=orig0;
p2=orig0+screw_axis;
[mc,nc]=size(coor);
dist=zeros(1,mc);
for k=1:mc,
    p0=coor(k,:);
    dist(k)=norm(cross((p0-p1),(p0-p2)))/norm(p2-p1);
end;
[distance,number]=min(dist);
orig=coor(number,:);



