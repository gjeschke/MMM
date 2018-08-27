function gobjects=plot_label(indices,plot_mode)
% function gobjects=plot_label(indices,plot_mode)
%
% symbolic graphics for a spin label, residue type must be a defined spin
% label type, otherwise nothing is displayed and an empty vector gobjects
% returned
%
% indices       indices for addressing a residue (1x4 vector)
% plot_mode     0: sphere at spin label position, radius and color are set
%               in the label definition
%               1: sphere as above plus sticks for magnetic frame, with x:
%               red, y: green, z: blue
%               defaults to 0 (only sphere)
% gobjects      vector of graphics handles objects
%
% G. Jeschke, 2009-2013

global model
global label_defs
global hMain
global graph_settings

chi3_color = false; % separate coloring for different chi3 values for MTSL

rad0=graph_settings.label_radius;

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;

gobjects=[];

if length(indices)<4, return; end;

if nargin<2, plot_mode=0; end;

res_adr=mk_address(indices(1:4));

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
restag=info.name;
if ~strcmpi(restag,'R1A'),
    chi3_color = false;
end;
lab_id=tag2id(restag,label_defs.restags);
if isempty(lab_id), return; end;
frame_tags=label_defs.residues(lab_id).frame;
spin_density = label_defs.residues(lab_id).spin_density;
col=label_defs.residues(lab_id).color/255;
% rad=label_defs.residues(lab_id).radius; % old implementation with fixed radius and transparency

orig_adr=sprintf('%s.%s:',res_adr,id2tag(1,frame_tags));
x_adr=sprintf('%s.%s:',res_adr,id2tag(2,frame_tags));
y_adr=sprintf('%s.%s:',res_adr,id2tag(3,frame_tags));
if chi3_color,
    chi3=get_chi3(res_adr);
end;

[orig_indices,msg]=resolve_address(orig_adr);

[x_indices,msg]=resolve_address(x_adr);
if plot_mode,
    [y_indices,msg]=resolve_address(y_adr);
end;

[m,n]=size(orig_indices);

if m==0, return; end;

populations=zeros(1,m);
if m>1,
    for k=1:m,
        [msg,pop]=get_location(orig_indices(k,:),'population');
        populations(k)=pop;
    end;
    populations=populations.^(1/3)/max(populations.^(1/3));
    populations=populations/sum(populations)^(1/3);
else
    populations=1;
end;

for k=1:m,
    pop=populations(k);
    if chi3_color,
        % fprintf(1,'Rotamer %i has chi3= %4.2f degree\n',k,180*chi3(k)/pi);
        if chi3(k)>0,
            col=[255 0 0]/255;
        else
            col=[0 0 255]/255;
        end;
    end;
    if isempty(spin_density)
        [~,orig]=get_location(orig_indices(k,:),'xyz');
        [~,xdir]=get_location(x_indices(k,:),'xyz');
        coor=(orig+xdir)/2;
    else
        spindens = spin_density(:,2);
        spindens = spindens/sum(spindens);
        [msd,~] = size(spin_density);
        coor0 = zeros(msd,3);
        for kk = 1:msd,
            tag = id2tag(spin_density(kk,1),label_defs.residues(lab_id).atoms);
            [aindices,~]=resolve_address([res_adr '.' tag]);
            [~,ccoor] =  get_location([aindices k],'xyz');
            coor0(kk,:) = ccoor;
        end;
        coor = spindens'*coor0;
    end;
    rad=pop*rad0;
    [x,y,z,t]=point2trisphere(coor,rad);
    obj=trisurf(t,x,y,z);
    % reducepatch(obj);
    xyz=[coor(1),coor(2),coor(3),2*rad,2*rad,2*rad];
    record_object(obj,indices,xyz,k);
    gobjects=[gobjects obj];
    set(obj, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
    set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
    if plot_mode,
        [msg,ydir]=get_location(y_indices(k,:),'xyz');
        my_frame=frame(orig,xdir,ydir);
        cols={'r','g','b'};
        for kk=1:3,
            vec=1.4*rad*my_frame(kk,:);
            coor2=coor+vec;
            h=plot3([coor(1),coor2(1)],[coor(2),coor2(2)],[coor(3),coor2(3)],cols{kk},'LineWidth',2);
            gobjects=[gobjects h];
        end;
    end;
end;

info.label_graphics.objects=gobjects;
info.label_graphics.mode=plot_mode+1;
info.label_graphics.color=col;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4))=info;

function chi3=get_chi3(res_adr)

CB_adr=sprintf('%s.%s:',res_adr,'CB');
[CB_indices,msg]=resolve_address(CB_adr);
SG_adr=sprintf('%s.%s:',res_adr,'SG');
[SG_indices,msg]=resolve_address(SG_adr);
SD_adr=sprintf('%s.%s:',res_adr,'SD');
[SD_indices,msg]=resolve_address(SD_adr);
CE_adr=sprintf('%s.%s:',res_adr,'CE');
[CE_indices,msg]=resolve_address(CE_adr);
[m,~]=size(CB_indices);
chi3=zeros(1,m);
for k=1:m,
    [msg,CB]=get_location(CB_indices(k,:),'xyz');
    [msg,SG]=get_location(SG_indices(k,:),'xyz');
    [msg,SD]=get_location(SD_indices(k,:),'xyz');
    [msg,CE]=get_location(CE_indices(k,:),'xyz');
    chi3(k)=dihedral(CB,SG,SD,CE);
end;