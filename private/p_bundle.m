function handles=p_bundle(handles,coor,falpha,colors)
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

global hMain

if nargin<3,
    falpha=1;
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

helix_colors=zeros(helices,3);
colormap(jet);
cmap=flipud(colormap);
if nargin>3,
    cmap=colors;
end;
[m,n]=size(cmap);
for k=1:helices,
    indy=1+round((m-1)*(k-1)/(helices-1));
    helix_colors(k,:)=cmap(indy,:);
end;

axes(handles.axes_model);
cla;
set(gca,'Clipping','off');
hold on

% figure(1); hold on;

for hnum=1:helices,
    basnum=4*(hnum-1);
    [x,y,z]=helix_tube(coor(basnum+2,:),coor(basnum+3,:));
    han=surf(x,y,z);
    set(han,'EdgeColor','none','FaceColor',helix_colors(hnum,:),'FaceAlpha',falpha,'FaceLighting','phong','Clipping','off');
    hold on;
    [x,y,z]=point2sphere(coor(basnum+2,:),rad_end_point);
    obj=surface(x,y,z);
    set(obj, 'FaceColor', color_cyto, 'EdgeColor', 'none','Clipping','off');
    [x,y,z]=point2sphere(coor(basnum+3,:),rad_end_point);
    obj=surface(x,y,z);
    set(obj, 'FaceColor', color_peri, 'EdgeColor', 'none','Clipping','off');
end;
hMain.camlight=camlight;
guidata(handles.axes_model,hMain);


lighting phong
material shiny
set(gca,'FontSize',14);
axis equal
axis off
view(0,90);
set(handles.button_cam_up,'Enable','off');
set(handles.text_view,'String','view 0 90');

guidata(handles.axes_model,handles);

axes(handles.axes_frame);
cla;
hold on

plot3([0 1],[0 0],[0 0],'r','LineWidth',2); % x axis
plot3([0 0],[0 1],[0 0],'g','LineWidth',2); % y axis
plot3([0 0],[0 0],[0 1],'b','LineWidth',2); % z axis
view(0,90);
axis([-0.1,1.1,-0.1,1.1,-0.1,1.1]);

guidata(handles.axes_frame,handles);


function [x,y,z]=helix_tube(acoor,ecoor)
% 
% Provides Cartesian coordinates for plotting a cylindrical tube with function surf
% when coordinates acoor of the starting point and ecoor of the end point
% are given

cent=(acoor+ecoor)/2;
hvec=ecoor-acoor;
[ofi,oth,hlen] = cart2sph(hvec(1),hvec(2),hvec(3));
oth=pi/2-oth;
rotmat=get_rotmat([ofi,oth,0]);
r=3.5;
[x,y,z]=cylinder(3.5);
z=hlen*(z-0.5*ones(size(z)));

[m,n]=size(x);
for k=1:m,
    for kk=1:n,
        vec=[x(k,kk) y(k,kk) z(k,kk)];
        vec=rot_trans(vec,rotmat',[0,0,0]);
        x(k,kk)=vec(1)+cent(1);
        y(k,kk)=vec(2)+cent(2);
        z(k,kk)=vec(3)+cent(3);
    end;
end;

function [x,y,z] = point2sphere(coor,radius)
%VEC2TUBE	Create parametric spherical surface at point with coordinates coor
%	[x,y,z] = point2sphere(coor,radius).
%	

%	(c) 2007 by Gunnar Jeschke
%

[x,y,z]=sphere;
[m,n]=size(x);
x=radius*x;
y=radius*y;
z=radius*z;

x=x+coor(1)*ones(m,n);
y=y+coor(2)*ones(m,n);
z=z+coor(3)*ones(m,n);

