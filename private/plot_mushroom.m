function gobjects=plot_mushroom(indices,plot_mode,label)
% function gobjects=plot_mushroom(indices,plot_mode)
%
% symbolic graphics for a spin label, residue type must be a defined spin
% label type, otherwise nothing is displayed and an empty vector gobjects
% returned
%
% indices       indices for addressing a residue (1x4 vector)
% plot_mode     0: point-spread ellipsoid
%               1: point-spread ellipsoid plus transparent ellipsoid for unrestrained
%               label
%               defaults to 0 (only sphere)
% gobjects      vector of graphics handles objects
%
% G. Jeschke, 2009-2013

global model
global hMain
global rotamer_libraries

col = [0.7,0.2,0];
alpha = 0.7;
col0 = [0.3,0.6,0];
alpha0 = 0.3;
scale = sqrt(2);

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;

gobjects=[];

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));

if length(indices)<4, return; end;

if ~exist('plot_mode','var') || isempty(plot_mode)
    plot_mode=0; 
end;

res_adr=mk_address(indices(1:4));
[msg,Cacoor] = get_object(sprintf('%s.CA',res_adr),'coor');

rot_lib_name = '';

for k = 1:length(rotamer_libraries),
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label),
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T),
            if rotamer_libraries(k).T(kk) == 298,
                id = kk;
            end;
        end;
        if id >0,
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
        end;
    end;
end;
if isempty(rot_lib_name)
    add_msg_board(sprintf('ERROR: Unknown label %s.',label));
    return
end

load(rot_lib_name);
label1 = rotamer_populations(indices,rot_lib);
[nr1,~] = size(label1.NOpos);
if nr1 < 2
    add_msg_board('Warning: At least two rotamers must exist. No display.',label);
    return
end
[pst0,center0] = point_spread_tensor(label1.NOpos(:,1:3),rot_lib.calibration.pop');
[dircos0,diag0] = eig(pst0);

[pst,center] = point_spread_tensor(label1.NOpos(:,1:3),label1.NOpos(:,4));
[dircos,diag1] = eig(pst);
if min(min(diag1)) <= 0 || min(min(diag0)) <= 0
    add_msg_board('Warning: No spread. Probably only one rotamer. No display.');
    return
end
[x,y,z] = ellipsoid(0,0,0,scale*sqrt(diag1(1,1)),scale*sqrt(diag1(2,2)),scale*sqrt(diag1(3,3)),20);
for k1 = 1:21,
    for k2 = 1:21,
        coor = [x(k1,k2) y(k1,k2) z(k1,k2)];
        coor = coor*dircos';
        x(k1,k2) = coor(1)+center(1);
        y(k1,k2) = coor(2)+center(2);
        z(k1,k2) = coor(3)+center(3);
    end;
end;

obj=surf(x,y,z);
xyz=[center(1),center(2),center(3),scale*sqrt(diag1(1,1)),scale*sqrt(diag1(2,2)),scale*sqrt(diag1(3,3))];
record_object(obj,indices,xyz,k);
gobjects=[gobjects obj];
set(obj, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha',alpha,'FaceLighting','gouraud','Clipping','off');
set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
obj2 = plot3([Cacoor(1),center(1)],[Cacoor(2),center(2)],[Cacoor(3),center(3)],'Color',col,'Linewidth',4);
gobjects=[gobjects obj2];

[short0,spoi0] = min(diag(diag0));
oblate0 = 1;
for k = 1:3
    if k ~= spoi0
        oblate0 = oblate0*diag0(k,k);
    end;
end;
oblate0 = sqrt(sqrt(oblate0)/short0);

[short,spoi] = min(diag(diag1));
shortdir = dircos(:,spoi);
% obj3 = plot3(center(1)+[0,3*shortdir(1)],center(2)+[0,3*shortdir(2)],center(3)+[0,3*shortdir(3)],'Color',[0,0,0.75],'Linewidth',4);
% gobjects=[gobjects obj3];

if plot_mode,
    [x,y,z] = ellipsoid(0,0,0,scale*sqrt(diag0(1,1)),scale*sqrt(diag0(2,2)),scale*sqrt(diag0(3,3)),20);
    for k1 = 1:21,
        for k2 = 1:21,
            coor = [x(k1,k2) y(k1,k2) z(k1,k2)];
            coor = coor*dircos0';
            x(k1,k2) = coor(1)+center0(1);
            y(k1,k2) = coor(2)+center0(2);
            z(k1,k2) = coor(3)+center0(3);
        end;
    end;
    obj=surf(x,y,z);
    xyz=[center0(1),center0(2),center0(3),scale*sqrt(diag0(1,1)),scale*sqrt(diag0(2,2)),scale*sqrt(diag0(3,3))];
    record_object(obj,indices,xyz,k);
    gobjects=[gobjects obj];
    set(obj, 'FaceColor', col0, 'EdgeColor', 'none', 'FaceAlpha',alpha0,'FaceLighting','gouraud','Clipping','off');
    set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
    obj2 = plot3([Cacoor(1),center0(1)],[Cacoor(2),center0(2)],[Cacoor(3),center0(3)],'Color',col0,'Linewidth',4);
    gobjects=[gobjects obj2];
end;

% Output diagnostics
spread = sqrt(prod(diag(diag1)));
oblate = 1;
for k = 1:3
    if k ~= spoi
        oblate = oblate*diag1(k,k);
    end;
end;
oblate = sqrt(sqrt(oblate)/short);

spread0 = sqrt(prod(diag(diag0)));
add_msg_board(sprintf('Label %s at residue %s.',label,res_adr));
add_msg_board(sprintf('Contraction is %5.3f',spread/spread0));
add_msg_board(sprintf('Oblaticity is %5.3f (for free label %5.3f)',oblate,oblate0));
stem0 = center0-Cacoor;
stem0e = stem0/norm(stem0);
stem = center-Cacoor;
steme = stem/norm(stem);
ang = 180*acos(sum(stem0e.*steme))/pi;
ang2 = 180*acos(sum(shortdir'.*steme))/pi;
if ang2 > 90
    ang2 = 180-ang2;
end;
add_msg_board(sprintf('Stalk reorientation is %5.1f°',ang));
add_msg_board(sprintf('Stalk vector is [%5.3f,%5.3f,%5.3f]',steme));
add_msg_board(sprintf('Cap normal vector is [%5.3f,%5.3f,%5.3f]',shortdir));
add_msg_board(sprintf('The two vectors deviate by %5.1f°.',ang2));



if isfield(info,'label_graphics') && isfield(info.label_graphics,'objects')
    info.label_graphics.objects=[info.label_graphics.objects gobjects];
else
    info.label_graphics.objects = gobjects;
end
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4))=info;