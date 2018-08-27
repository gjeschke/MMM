function unrecord_objects
% function unrecord_objects
%
% Removes graphics objects from the list of stored objects
% consolidation of the objects list, deleted objects are stored in field
% model.unrecord, which is emptied after removal
%
%
% G. Jeschke, 2009

global model

if ~isfield(model,'graphics_lookup'),
    return;
end;

[m,~]=size(model.graphics_lookup);
% keyboard
poi0=model.graphics_lookup_pointer;
poi1=0;
retained_objects=zeros(poi0,7);
retained_gobjects=gobjects(1,poi0);
retained_xyz=zeros(poi0,6);

for k=1:poi0,
    unrecord=find(model.unrecord==model.graphics_objects(k), 1);
    if isempty(unrecord)
        poi1=poi1+1;
        retained_gobjects(poi1) = model.graphics_objects(k);
        retained_objects(poi1,2:7)=model.graphics_lookup(k,2:7);
        retained_xyz(poi1,:)=model.graphics_xyz(k,:);
    end;
end;
model.graphics_objects(1:poi1) = retained_gobjects(1:poi1);
model.graphics_lookup(1:poi1,2:7)=retained_objects(1:poi1,2:7);
model.graphics_xyz(1:poi1,:)=retained_xyz(1:poi1,:);
if poi1<m,
    model.graphics_lookup(poi1+1:m,:)=0;
    model.graphics_xyz(poi1+1:m,:)=0;
end;
model.unrecord=[];
model.graphics_lookup_pointer=poi1;