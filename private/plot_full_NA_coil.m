function [gobjects,x,y,z]=plot_full_NA_coil(indices)
% Plots a smooth coil for DNA/RNA addressed by structure, chain and
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
                if ~isempty(model.secondary_selected)
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
% Extract C1', C4' and C3' coordinates
for k=1:n,
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(k);
    if isfield(info,'terminal'),
        if info.terminal, terminate=1; else terminate=0;
        end;
    else
        terminate=0;
    end;
    C1_id=tag2id('C1''',info.atom_tags);
    if isempty(C1_id),
        C1_id=tag2id('C1*',info.atom_tags);
    end;
    C3_id=tag2id('C3''',info.atom_tags);
    if isempty(C4_id),
        C3_id=tag2id('C3*',info.atom_tags);
    end;
    C4_id=tag2id('C4''',info.atom_tags);
    if isempty(C4_id),
        C4_id=tag2id('C4*',info.atom_tags);
    end;
    if ~isempty(C1_id) && ~isempty(C3_id) && ~isempty(C4_id),
        poi=poi+1;
        resnum(poi)=k;
        C1_num=info.atom_numbers{C1_id};
        [mm,nn]=size(C1_num);
        if nn==1, C1_num=[C1_num 1]; end;
        C1_coor=zeros(1,3);
        pop=0;
        for kk=1:mm,
            C1_coor=C1_coor+C1_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(C1_num(kk,1),:);
            pop=pop+C1_num(kk,2);
        end;
        C1_coor=C1_coor/pop;
        C1(poi,:)=C1_coor;
        C3_num=info.atom_numbers{C3_id};
        [mm,nn]=size(C3_num);
        if nn==1, C3_num=[C3_num 1]; end;
        C3_coor=zeros(1,3);
        pop=0;
        for kk=1:mm,
            C3_coor=C3_coor+C3_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(C3_num(kk,1),:);
            pop=pop+C3_num(kk,2);
        end;
        C3_coor=C3_coor/pop;
        C3(poi,:)=C3_coor;
        C4_num=info.atom_numbers{C4_id};
        [mm,nn]=size(C4_num);
        if nn==1, C4_num=[C4_num 1]; end;
        C4_coor=zeros(1,3);
        pop=0;
        for kk=1:mm,
            C4_coor=C4_coor+C4_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(C4_num(kk,1),:);
            pop=pop+C4_num(kk,2);
        end;
        C4_coor=C4_coor/pop;
        C4(poi,:)=C4_coor;
    end;
    if terminate, break; end;
end;
col=graph_settings.NA_color;
if poi>0 && sum(sum(isnan(C3)))==0 && sum(sum(isnan(C4)))==0,
    C1=C1(1:poi,:);
    C3=C3(1:poi,:);
    C4=C4(1:poi,:);
    [backbone,rung,normal]=mk_ribbon(graph_settings.spr,(C3+C4)/2,C1);
    if ~isempty(backbone)
        [x,y,z]=generalized_tubeplot(backbone',rung',normal',3,graph_settings.NA_spr);
        han=surf(x,y,z);
        set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
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