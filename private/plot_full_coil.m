function [gobjects,x,y,z]=plot_full_coil(indices)
% Plots a smooth coil addressed by structure, chain and
% coordinate set indices
%
% indices   first three indices into structure
% gobjects  graphics objects
% x,y,z     extensions in x,y,z

global model
global graph_settings
global hMain

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;

falpha=1;

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'secondary_graphics'),
    allgraphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics;
    if ~isempty(allgraphics),
        for k=1:length(allgraphics),
            graphics=allgraphics(k);
            if ~isempty(graphics)
                if ~isempty(graphics.objects),
                    for kk=1:length(graphics.objects),
                        if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                            delete(graphics.objects(kk));
                            model.unrecord=[model.unrecord graphics.objects(kk)];
                        end;
                    end;
                end;
                if isfield(model,'secondary_selected') && ~isempty(model.secondary_selected)
                    poi=0;
                    for kk=1:length(model.secondary_selected),
                        if isempty(find(graphics.objects==model.secondary_selected{kk}.objects, 1)),
                            poi=poi+1;
                            model.secondary_selected{poi}=model.secondary_selected{kk};
                        end;
                    end;
                    model.secondary_selected=model.secondary_selected{1:poi};
                end;
            end;
        end;
    end;
end;

model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics=[];

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'chain_model_graphics'),
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics;
    if ~isempty(graphics)
        if ~isempty(graphics.objects),
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    delete(graphics.objects(k));
                    model.unrecord=[model.unrecord graphics.objects(k)];
                end;
            end;
        end;
    end;
end;

gobjects=[];

n=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);
CA=zeros(n,3);
O=zeros(n,3);
poi=0;
% Extract Calpha and O coordinates
for k=1:n,
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(k);
    if isfield(info,'terminal'),
        if info.terminal, terminate=1; else terminate=0;
        end;
    else
        terminate=0;
    end;
    CA_id=tag2id('CA',info.atom_tags);
    O_id=tag2id('O',info.atom_tags);
    if ~isempty(CA_id) && ~isempty(O_id),
        poi=poi+1;
        CA_num=info.atom_numbers{CA_id};
        [mm,nn]=size(CA_num);
        if nn==1, CA_num=[CA_num 1]; end;
        CA_coor=zeros(1,3);
        pop=0;
        for kk=1:mm,
            CA_coor=CA_coor+CA_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(CA_num(kk,1),:);
            pop=pop+CA_num(kk,2);
        end;
        CA_coor=CA_coor/pop;
        CA(poi,:)=CA_coor;
        O_num=info.atom_numbers{O_id};
        [mm,nn]=size(O_num);
        if nn==1, O_num=[O_num 1]; end;
        O_coor=zeros(1,3);
        pop=0;
        for kk=1:mm,
            O_coor=O_coor+O_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(O_num(kk,1),:);
            pop=pop+O_num(kk,2);
        end;
        O_coor=O_coor/pop;
        O(poi,:)=O_coor;
    end;
    if terminate, break; end;
end;
if poi>0,
    CA=CA(1:poi,:);
    O=O(1:poi,:);
    [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,0);
    if ~isempty(backbone)
        [x,y,z]=tubeplot(backbone',graph_settings.coil_radius);
        han=surf(x,y,z);
        set(han,'EdgeColor','none','FaceColor',graph_settings.coil_color,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
        set(han, 'CDataMapping','direct','AlphaDataMapping','none');
        gobjects=[gobjects han];
        minx=min([min(x),minx]);
        maxx=max([max(x),maxx]);
        miny=min([min(y),miny]);
        maxy=max([max(y),maxy]);
        minz=min([min(z),minz]);
        maxz=max([max(z),maxz]);        
    end;
end;

graphics.mode=3;
graphics.objects=gobjects;    
graphics.color=graph_settings.coil_color;
graphics.opaque=falpha;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics=graphics;

for k=1:length(gobjects),
    if ~isempty(gobjects(k)),
        if gobjects(k)~=0,
            record_object(gobjects(k),indices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
        end;
    end;
end;

x=[minx,maxx];
y=[miny,maxy];
z=[minz,maxz];