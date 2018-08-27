function handles=depth_cueing(handles)

global hMain
global model
global hModel

set(hMain.figure,'Pointer','watch');
if hMain.detached,
    set(hModel.figure,'Pointer','watch');
end;
if hMain.hierarchy_display,
%     if hMain.large,
%         set(hMain.hierarchy_window_large,'Pointer','watch');
%     else
%         set(hMain.hierarchy_window,'Pointer','watch');
%     end;
    set(hMain.hierarchy_window,'Pointer','watch');
end;

switch hMain.color,
    case 'white'
        bckg=[1,1,1];
    case 'grey'
        bckg=[0.941,0.941,0.941];
    case 'black'
        bckg=[0,0,0];
end;

m=model.graphics_lookup_pointer;
campos=get(hMain.axes_model,'CameraPosition');
camtarget=get(hMain.axes_model,'CameraTarget');
camvec=camtarget-campos;
camvec=camvec/norm(camvec);
z=zeros(1,m);
xyzmat=zeros(m,3);
for k=1:m,
    xyzext=model.graphics_xyz(k,:);
    xyz=(xyzext(1:3)+xyzext(4:6))/2;
    xyzmat(k,:)=xyz;
    z(k)=norm(xyz.*camvec);
end;

minz=min(z);
maxz=max(z);

selection_flag=0;
if isfield(model,'selections')
    if model.selections>=1,
        selection_flag=1;
        z2=zeros(1,model.selections);
        for k=1:model.selections,
            cindices=model.selected{k};
            cindices=cindices(cindices>0);
            [msg,xyz]=get_object(cindices,'xyz');
            z2(k)=norm(mean(xyz).*camvec);
        end;
        minz=min(z2)-5;
        maxz=max(z2);
    end;
end;

dz=(z-minz)/(maxz-minz);

backplane=1;
focusplane=0;
dz=(dz-focusplane)/(backplane-focusplane);
trans=ones(size(dz));
trans(find(dz>backplane))=0;
dz(find(dz<0))=0;
dz(find(dz>1))=1;
db=1-dz;
ambient=db;
if selection_flag,
    ambient=sqrt(db);
end;
for kk=1:m,
    indices=model.graphics_lookup(kk,2:n);
    indices=indices(indices>0);
    [msg,graphics]=get_object(indices,'graphics');
    if ~isempty(graphics),
        setcolor=db(kk)*bckg+dz(kk)*graphics.color(1,:);
        if graphics.mode==1,
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    set(graphics.objects(k),'Color',setcolor,'Alpha',trans(kk));
                end;
            end;
        elseif graphics.mode>1,
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    set(graphics.objects(k),'FaceColor',setcolor,'FaceAlpha',trans(kk),'AmbientStrength',ambient(kk));
                end;
            end;
        end;
    end;
end;

set(hMain.figure,'Pointer','arrow');
if hMain.detached,
    set(hModel.figure,'Pointer','arrow');
end;
if hMain.hierarchy_display,
%     if hMain.large,
%         set(hMain.hierarchy_window_large,'Pointer','arrow');
%     else
%         set(hMain.hierarchy_window,'Pointer','arrow');
%     end;
    set(hMain.hierarchy_window,'Pointer','arrow');
end;
