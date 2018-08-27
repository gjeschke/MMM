function [err,helix_colors]=p_bundle_fit(coor,falpha,colors)
%
% Creates a cartoon plot of a helix bundle without spin labels
%
% helices are colored in a rainbow sequence
% each helix is defined by 4x3 Cartesian coordinates
%   1: cytoplasmic label
%   2: cytoplasmic end
%   3: periplasmic end
%   4: periplasmic label
%
% program returns with error code 1 if number of coordinate sets is not a
% multiple of four
%
% falpha    optional argument for transparency of the helix cylinder
%           surfaces, 0 fully transparent (invisible), 1 (default) fully
%           opaque
% colors    optional argument that specifies the colors of all helices as a
%           RGB triple, if coor has the dimension 4Nx3, colors must have
%           the dimension Nx3, otherwise interpolation is used
%
% (c) G. Jeschke, 2008

mask = zeros(1,13);
mask(1) = 1;
mask(12) = 1;
mask(13) = 1;


cytoplasmic=[25 41 91 128 179 190 253 276 349 371 423 430 470];
periplasmic=[5 62 80 149 164 210 233 294 326 391 400 446 450];
helices=[2,3,4,5,6,7,8,9,10,11];

bundle_coor = zeros(length(helices),3);
prot_coor = bundle_coor;
for k = 1:length(helices),
    bas1 = 4*(helices(k)-1);
    bas2 = 2*(k-1);
    bundle_coor(bas2+1,:) = coor(bas1+2,:);
    bundle_coor(bas2+2,:) = coor(bas1+3,:);
    adr1 = sprintf('%i.CA',cytoplasmic(helices(k)));
    [~,ccoor] = get_object(adr1,'coor');
    prot_coor(bas2+1,:) = ccoor;
    adr2 = sprintf('%i.CA',periplasmic(helices(k)));
    [~,pcoor] = get_object(adr2,'coor');
    prot_coor(bas2+2,:) = pcoor;
end;

[rms,~,transmat]=rmsd_superimpose(prot_coor,bundle_coor);
fprintf(1,'Helix end point r.m.s.d is %5.2f Å\n',rms);

[mm,~] = size(coor);
xyz=[coor ones(mm,1)];
xyz=transmat*xyz';
xyz=xyz';
coor=xyz(:,1:3);

if nargin<2,
    falpha=0.75;
end;

% Standard settings for plot appearance
rad_end_point=0.5;
rad_label_arm=0.3;
color_arm=[0.85 0.85 0.85];
rad_label=1.0;
color_end=[0.5 0.5 0.5];
color_cyto=[1.0 0.15 0.15];
color_peri=[0.15 0.15 1.0];


err=0;

[m,n]=size(coor);
if mod(m,4)~=0,
    disp('### Number of coordinate sets is not a multiple of four ###');
    return;
end;

helices=m/4;

if helices<=3,
    helix_colors=flipud(hsv(helices));
else
    helix_colors=jet(helices);
end;
if nargin>2,
    cmap=colors;
    helix_colors=zeros(helices,3);
    [m,n]=size(cmap);
    for k=1:helices,
        indy=1+round((m-1)*(k-1)/(helices-1));
        helix_colors(k,:)=cmap(indy,:);
    end;
end;

for hnum=1:helices,
    if mask(hnum),
        basnum=4*(hnum-1);
        [x,y,z]=helix_tube(coor(basnum+2,:),coor(basnum+3,:));
        han=surf(x,y,z);
        set(han,'EdgeColor','none','FaceColor',helix_colors(hnum,:),'FaceAlpha',falpha,'FaceLighting','phong');
        hold on;
        [x,y,z]=point2sphere(coor(basnum+2,:),rad_end_point);
        obj=surface(x,y,z);
        set(obj, 'FaceColor', color_cyto, 'EdgeColor', 'none');
        [x,y,z]=point2sphere(coor(basnum+3,:),rad_end_point);
        obj=surface(x,y,z);
        set(obj, 'FaceColor', color_peri, 'EdgeColor', 'none');
    end;
end;
% camlight
% lighting phong
% material dull
% set(gca,'FontSize',14);
% axis equal
% axis off
% view(0,90);

function [x,y,z]=helix_tube(acoor,ecoor)
% 
% Provides Cartesian coordinates for plotting a cylindrical tube with function surf
% when coordinates acoor of the starting point and ecoor of the end point
% are given

cent=(acoor+ecoor)/2;
hvec=ecoor-acoor;
[ofi,oth,hlen] = cart2sph(hvec(1),hvec(2),hvec(3));
oth=pi/2-oth;
rotmat=get_rotmat_old(180*[ofi,oth,0]/pi);
r=3.5;
[x,y,z]=cylinder(3.5);
z=hlen*(z-0.5*ones(size(z)));

[m,n]=size(x);
for k=1:m,
    for kk=1:n,
        vec=[x(k,kk) y(k,kk) z(k,kk)];
        vec=rot_trans_old(vec,rotmat',[0,0,0]);
        x(k,kk)=vec(1)+cent(1);
        y(k,kk)=vec(2)+cent(2);
        z(k,kk)=vec(3)+cent(3);
    end;
end;

function Rp=get_rotmat_old(euler)
%
% function Rp=get_rotmat(euler),
%
% Compute rotation matrix
%
% euler     Euler angles in degree
% Rp        rotation matrix
%

euler=pi*euler/180;
ca=cos(euler(1));
sa=sin(euler(1));
cb=cos(euler(2));
sb=sin(euler(2));
cg=cos(euler(3));
sg=sin(euler(3));

Rp=[ca*cb*cg-sa*sg,sa*cb*cg+ca*sg,-sb*cg;-ca*cb*sg-sa*cg,-sa*cb*sg+ca*cg,sb*sg;ca*sb,sa*sb,cb];

function newcoor=rot_trans_old(coor,rotmat,t)
%
% newcoor=rot_trans(coor,rotmat,t),
%
% Rotation and subsequent translation of a coordinate array
% 
% coor      x,y,z coordinate array
% rotmat    rotation matrix, use get_rotmat to create from Euler angles
% t         vector for coordinate translation (xt,yt,zt)
%

coor=coor';
[m,n]=size(coor);
newcoor=coor;
for k=1:n,
    newcoor(:,k)=rotmat*coor(:,k)+t';
end;
newcoor=newcoor'; 