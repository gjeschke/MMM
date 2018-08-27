function gobjects=plot_full_ribbon(indices)
% Plots a smooth coil addressed by structure, chain and
% coordinate set indices
%
% indices   first three indices into structure
% objects   graphics objects
% x,y,z     extensions in x,y,z

global model
global graph_settings

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;

falpha=1;

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
            end;
        end;
    end;
end;

model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics=[];

gobjects=[];
colors=zeros(1000,3);
secondary=zeros(1,1000);
range=zeros(1000,2);
minmax=zeros(1000,6);

n=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);
CA=zeros(n,3);
O=zeros(n,3);
poi=0;
cpoi=0;
starter=1;
% Extract Calpha and O coordinates
sec=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(1).secondary;
for k=1:n,
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(k);
    CA_id=tag2id('CA',info.atom_tags);
    O_id=tag2id('O',info.atom_tags);
    secend=1;
    if ~isempty(CA_id) && ~isempty(O_id),
        poi=poi+1;
        secend=0;
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
    if info.secondary~=sec && poi~=0,
        switch sec
            case 0
                col=graph_settings.coil_color;
            case 1
                col=graph_settings.helix_color;
            case 2
                col=graph_settings.sheet_color;
        end;
        CA=CA(1:poi,:);
        O=O(1:poi,:);
        if poi>1,
            CA(poi,:)=(CA(poi-1,:)+CA(poi,:))/2;
            O(poi,:)=(O(poi-1,:)+O(poi,:))/2;
        end;
        if sec==2, sheetflag=1; else sheetflag=0; end;
        [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,sheetflag);
        if ~isempty(backbone)
            switch sec
                case 0
                    col=graph_settings.loop_color;
                    [x,y,z]=tubeplot(backbone',graph_settings.coil_radius/2);
                case 1
                    col=graph_settings.helix_color;
                    [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec);
               case 2
                    col=graph_settings.sheet_color;
                    [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec,graph_settings.spr);
            end;
            han=surf(x,y,z);
            set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
            gobjects=[gobjects han];
            cpoi=cpoi+1;
            colors(cpoi,:)=col;
            range(cpoi,1)=starter;
            range(cpoi,2)=k-1;
            secondary(k)=sec;
            minx=min(min(x));
            maxx=max(max(x));
            miny=min(min(y));
            maxy=max(max(y));
            minz=min(min(z));
            maxz=max(max(z));
            minmax(cpoi,:)=[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz];
            starter=k;
        end;
        if secend,
            poi=0;
        else
            if poi>1,
                CA(1,:)=(CA(poi-1,:)+CA(poi,:))/2;
                O(1,:)=(O(poi-1,:)+O(poi,:))/2;
            else
                CA(1,:)=CA(poi,:);
                O(1,:)=O(poi,:);
            end;
            poi=1;
        end;
        sec=info.secondary;
    end;
end;
if poi>0,
    CA=CA(1:poi,:);
    O=O(1:poi,:);
    if sec==2, sheetflag=1; else sheetflag=0; end;
    [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,sheetflag);
    if ~isempty(backbone)
        switch sec
            case 0
                col=graph_settings.loop_color;
                [x,y,z]=tubeplot(backbone',graph_settings.coil_radius/2);
            case 1
                col=graph_settings.helix_color;
                [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec);
           case 2
                col=graph_settings.sheet_color;
                [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec,graph_settings.spr);
        end;
        han=surf(x,y,z);
        set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
        gobjects=[gobjects han];
        cpoi=cpoi+1;
        colors(cpoi,:)=col;
        range(cpoi,1)=starter;
        range(cpoi,2)=k-1;
        secondary(k)=sec;
        minx=min(min(x));
        maxx=max(max(x));
        miny=min(min(y));
        maxy=max(max(y));
        minz=min(min(z));
        maxz=max(max(z));
        minmax(cpoi,:)=[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz];
    end;
end;

graphics.mode=4;
graphics.opaque=falpha;

allgraphics=[];

selflag=is_selected(indices(1:3));
selflag=0;

for k=1:length(gobjects),
    if ~isempty(gobjects(k)),
        if gobjects(k)~=0,
            pindices=[indices -k];
            record_object(gobjects(k),pindices,minmax(k,:));
            graphics.objects=gobjects(k);    
            graphics.color=colors(k,:);
            graphics.sec=secondary(k);
            graphics.range=range(k,:);
            allgraphics=[allgraphics graphics];
            if selflag,
                if isempty(model.secondary_selected),
                    model.secondary_selected{1}=graphics;
                    model.secondary_indices{1}=pindices;
                else
                    n=length(model.secondary_selected);
                    model.secondary_selected{n+1}=graphics;
                    model.secondary_indices{n+1}=pindices;
                end;                                                
            end;
        end;
    end;
end;

model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics=allgraphics;
