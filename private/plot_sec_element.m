function [objects,x,y,z,CA]=plot_sec_element(indices,range,sec_type)
% Plots a secondary structure element addressed by structure, chain and
% coordinate set indices, range of residues in the chain, and secondary
% structure type and returns object handles
%
% indices   first three indices into structure
% range     range of residues in chain (internal indexing)
% sec_type  0 loop, 1 helix, 2,-2 strand (sheet), -1 helix cartoon (cylinder)

global model
global graph_settings

objects=[];

range(2) = range(2)+1;
n=range(end)-range(1)+1;
CA=zeros(n,3);
O=zeros(n,3);
% Extract Calpha and O coordinates and residue indices
rindices = zeros(n,4);
rindices(:,1:3) = repmat(indices(1:3),n,1);
cmadr = mk_address(indices(1:3));
for k=range(1):range(2),
    ind=k-range(1)+1;
    rindices(ind,4) = k;
    resadr = sprintf('%s%i',cmadr,k);
    [msg,CA_coor] = get_object([resadr '.CA'],'coor');
    CA(ind,:)=CA_coor;
    [msg,O_coor] = get_object([resadr '.O'],'coor');
    O(ind,:)=O_coor;
end;
CA=CA(1:ind,:);
O=O(1:ind,:);
if ind>1,
    CA(ind,:)=(CA(ind-1,:)+CA(ind,:))/2;
    O(ind,:)=(O(ind-1,:)+O(ind,:))/2;
end;
if sec_type==2,
    sheetflag=1;
else
    sheetflag=0;
end;
[backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,sheetflag);
if ~isempty(backbone)
    switch sec_type
        case 0
            [x,y,z] = tubeplot(backbone',graph_settings.coil_radius);
            objects = surf(x,y,z);
            set(objects,'EdgeColor','none','FaceColor',graph_settings.coil_color,'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
            set(objects, 'CDataMapping','direct','AlphaDataMapping','none');
        case -1
            [objects,x,y,z,hends]=plot_helix_cartoon(CA); %graph_settings.helix_color
            CA = hends;
        case 1
            [objects,x,y,z]=plot_helix(backbone,rung,normal,[1,0,0]); %graph_settings.helix_color
        case {-2,2}
            [x,y,z]=generalized_tubeplot(backbone',rung',normal',abs(sec_type),graph_settings.spr);
            objects = surf(x,y,z);
            set(objects,'EdgeColor','none','FaceColor',[1,0,0],'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
            set(objects, 'CDataMapping','direct','AlphaDataMapping','none');
    end;
end;