function record_cartoon_object(obj,indices,xyz,num)
% function record_cartoon_object(obj,indices,xyz,num)
%
% Stores graphics object handles and corresponding indices of structure
% elements in a central lookup table and sets the function that is executed
% on clicking an object, the extensions of the object are also
% stored
%
% gobjects  vector of object handles
% indices   corresponding index vector
% xyz       midpoint and extension of graphics object, midpoint values
%           x,y,z, extensions max(x)-min(x), max(y)-min(y), max(z)-min(z)
%           as 1*6 vector, storage is as minima and maxima of coordinates
% num       auxiliary index for label location
%
% G. Jeschke, 2009

global model

maxx=xyz(1)+xyz(4)/2;
minx=xyz(1)-xyz(4)/2;
maxy=xyz(2)+xyz(5)/2;
miny=xyz(2)-xyz(5)/2;
maxz=xyz(3)+xyz(6)/2;
minz=xyz(3)-xyz(6)/2;

set(obj, 'ButtonDownFcn',@cartoon_object_clicked);

if isfield(model,'graphics_lookup')
    poi=model.graphics_lookup_pointer;
else
    poi=0;
    model.graphics_objects=gobjects(1,50000);
    model.graphics_lookup=zeros(50000,7);
    model.graphics_lookup_pointer=0;
end;

if ~isfield(model,'graphics_xyz')
    model.graphics_xyz=zeros(50000,6);
end;

if obj~=0 % discard empty objects
    poi=poi+1;
    model.graphics_objects(poi)=obj;
    model.graphics_lookup(poi,2:length(indices)+1)=indices;
    if nargin>3,
        model.graphics_lookup(poi,7)=num;
    end;
    model.graphics_xyz(poi,:)=[minx miny minz maxx maxy maxz];
end;

model.graphics_lookup_pointer=poi;

    