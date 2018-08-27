function objects = plot_bilabel(popcoor,indices,col)
% Plots the spin label cloud for a label specified by Cartesian coordinates
% (first three columns of popcoor) and populations (fourth column)
% indices of the current chain must be supplied
% a color can be specified by col (RGB triple, default royalblue)
%
% handles to the graphics objects are returned
%
% G. Jeschke, 27.3.2018

global graph_settings

rad0 = graph_settings.label_radius;

if ~exist('col','var') || isempty(col)
    col = [65,105,225]/255;
end

[m,~] = size(popcoor);

objects = gobjects(1,m);

for k = 1:m
    if popcoor(k,4) > 0.002,
        rad = popcoor(k,4)^(1/3)*rad0;
        [x,y,z,t]=point2trisphere(popcoor(k,1:3),rad);
        obj=trisurf(t,x,y,z);
        % reducepatch(obj);
        xyz=[popcoor(k,1),popcoor(k,2),popcoor(k,3),2*rad,2*rad,2*rad];
        record_object(obj,indices,xyz,k);
        objects(k) = obj; 
        set(obj, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
        set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
    end
end