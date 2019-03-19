function record_cartoon_object(obj,indices,xyz,num)
% function record_cartoon_object(obj,indices,xyz,num)
%
% Stores indices of structure elements in the object's UserData and
% sets the callback for clicking an object, the extensions of the
% object are also stored
%
% obj       graphics object handle
% indices   corresponding index vector
% xyz       midpoint and extension of graphics object, midpoint values
%           x,y,z, extensions max(x)-min(x), max(y)-min(y), max(z)-min(z)
%           as 1*6 vector, storage is as minima and maxima of coordinates
% num       auxiliary index for label location
%
% G. Jeschke, 2009

if obj==0, return;  end

obj.ButtonDownFcn = @cartoon_object_clicked;
obj.UserData.lookup(2:length(indices)+1) = indices;
if nargin>3
    obj.UserData.lookup(7)=num;
end

% Calculate graphics object extent and store it
minx=xyz(1)-xyz(4)/2;
miny=xyz(2)-xyz(5)/2;
minz=xyz(3)-xyz(6)/2;
maxx=xyz(1)+xyz(4)/2;
maxy=xyz(2)+xyz(5)/2;
maxz=xyz(3)+xyz(6)/2;
xyzrange = [minx miny minz maxx maxy maxz];
obj.UserData.xyz = xyzrange;
