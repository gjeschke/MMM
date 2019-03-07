function gobjects=plot_tiled_ribbon(indices)
% Plots a tiled ribbon addressed by structure, chain and
% coordinate set indices and assigns the graphics objects to residues
%
% indices   first three indices into structure
% objects   graphics objects
% x,y,z     extensions in x,y,z

global model
global graph_settings
global hMain

if gca~=hMain.axes_model
    axes(hMain.axes_model);
end

maxx=-1e6;
maxy=-1e6;
maxz=-1e6;
minx=1e6;
miny=1e6;
minz=1e6;

falpha=1;

if isfield(model.structures{indices(1)}(indices(2)),'helix_defs') && ~isempty(model.structures{indices(1)}(indices(2)).helix_defs),
    helix_defs=model.structures{indices(1)}(indices(2)).helix_defs;
else
    helix_defs=[];
end

if isfield(model.structures{indices(1)}(indices(2)),'sheet_defs') && ~isempty(model.structures{indices(1)}(indices(2)).sheet_defs),
    sheet_defs=model.structures{indices(1)}(indices(2)).sheet_defs;
else
    sheet_defs=[];
end


if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'chain_model_graphics')
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics;
    if ~isempty(graphics)
        if ~isempty(graphics.objects)
            for k=1:length(graphics.objects)
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0
                    delete(graphics.objects(k));
                    model.unrecord=[model.unrecord graphics.objects(k)];
                end
            end
        end
    end
end

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'secondary_graphics')
    allgraphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics;
    if ~isempty(allgraphics)
        for k=1:length(allgraphics)
            graphics=allgraphics(k);
            if ~isempty(graphics)
                if ~isempty(graphics.objects)
                    for kk=1:length(graphics.objects)
                        if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0
                            delete(graphics.objects(kk));
                            model.unrecord=[model.unrecord graphics.objects(kk)];
                        end
                    end
                end
            end
        end
    end
end

graphics.objects=[];

gobjects=[];
colors=zeros(1000,3);
secondary=zeros(1,1000);
range=zeros(1000,2);
minmax=zeros(1000,6);

if isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info)
    return;
end

n=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);

TM_flags=zeros(1,n);
for k=1:n
    mynum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(k).number;
    if ~isempty(mynum) && ~isnan(mynum)
        if ~isempty(helix_defs)
            for kk=1:length(helix_defs)
                if isfield(helix_defs{kk},'TM') && helix_defs{kk}.TM
                    if mynum>=helix_defs{kk}.range(1) && mynum<=helix_defs{kk}.range(2)
                        TM_flags(k)=1;
                    end
                end
            end
        end
        if ~isempty(sheet_defs)
            for kk=1:length(sheet_defs)
                if isfield(sheet_defs{kk},'TM') && sheet_defs{kk}.TM
                    if mynum>=sheet_defs{kk}.range(1) && mynum<=sheet_defs{kk}.range(2)
                        TM_flags(k)=1;
                    end
                end
            end
        end
    end
end

% Set up TM flags

CA=zeros(n,3);
O=zeros(n,3);
resnum=zeros(1,n);
poi=0;
cpoi=0;
starter=1;
% Extract Calpha and O coordinates
sec=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(1).secondary;
for k=1:n
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(k);
    if isfield(info,'terminal')
        if info.terminal, terminate=1; else terminate=0;
        end
    else
        terminate=0;
    end
    CA_id=tag2id('CA',info.atom_tags);
    O_id=tag2id('O',info.atom_tags);
    secend=1;
    if ~isempty(CA_id) && ~isempty(O_id)
        poi=poi+1;
        secend=0;
        CA_num=info.atom_numbers{CA_id};
        [mm,nn]=size(CA_num);
        if nn==1, CA_num=[CA_num 1]; end
        CA_coor=zeros(1,3);
        pop=0;
        for kk=1:mm
            CA_coor=CA_coor+CA_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(CA_num(kk,1),:);
            pop=pop+CA_num(kk,2);
        end
        CA_coor=CA_coor/pop;
        CA(poi,:)=CA_coor;
        resnum(poi)=k;
        O_num=info.atom_numbers{O_id};
        [mm,nn]=size(O_num);
        if nn==1, O_num=[O_num 1]; end
        O_coor=zeros(1,3);
        pop=0;
        for kk=1:mm
            O_coor=O_coor+O_num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(O_num(kk,1),:);
            pop=pop+O_num(kk,2);
        end
        O_coor=O_coor/pop;
        O(poi,:)=O_coor;
    end
    if info.secondary~=sec && poi~=0
        switch sec
            case 0
                col=graph_settings.coil_color;
            case 1
                col=get_color(sec,TM_flags(k-1));
            case 2
                col=get_color(sec,TM_flags(k-1));
        end
        CA=CA(1:poi,:);
        O=O(1:poi,:);
        if poi>1
            CA(poi,:)=(CA(poi-1,:)+CA(poi,:))/2;
            O(poi,:)=(O(poi-1,:)+O(poi,:))/2;
        end
        if sec==2, sheetflag=1; else, sheetflag=0; end
        [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,sheetflag);
        if ~isempty(backbone)
            switch sec
                case 0
                    col=graph_settings.loop_color;
                    [x,y,z]=tubeplot(backbone',graph_settings.coil_radius/2);
                case 1
                    col=get_color(sec,TM_flags(k-1));
                    [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec);
               case 2
                    col=get_color(sec,TM_flags(k-1));
                    [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec,graph_settings.spr);
            end
            [mm,nn]=size(x);
            for kk=1:length(CA)-1
                rnum=resnum(kk);
                aa=(kk-1)*graph_settings.spr;
                ee=aa+graph_settings.spr+2;
                if aa<1, aa=1; end
                if ee>nn, ee=nn; end
                if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum),'residue_graphics')
                    if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics,'objects')
                        if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects)
                            for kkk=1:length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects)
                                if ishandle(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects(kkk))
                                    delete(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects(kkk));
                                end
                            end
                        end
                    end
                end
                if ee-aa>=2
                    han=surface(x(:,aa:ee),y(:,aa:ee),z(:,aa:ee));
                    set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
                    set(han, 'CDataMapping','direct');
                    minx=min(min(x(:,aa:ee)));
                    maxx=max(max(x(:,aa:ee)));
                    miny=min(min(y(:,aa:ee)));
                    maxy=max(max(y(:,aa:ee)));
                    minz=min(min(z(:,aa:ee)));
                    maxz=max(max(z(:,aa:ee)));
                    rindices=[indices(1:3) rnum];
                    record_object(han,rindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
                    residue_graphics.mode=4;
                    residue_graphics.objects=han;    
                    residue_graphics.color=col;
                    residue_graphics.opaque=falpha;
                    model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics=residue_graphics;
                end
            end
            starter=k;
        end
        if secend
            poi=0;
        else
            CA(1,:)=CA(poi,:);
            O(1,:)=O(poi,:);
            poi=1;
            resnum(1)=k;
        end
        sec=info.secondary;
    end
    if terminate, break; end
end
if poi>0 && sum(sum(isnan(CA)))==0 && sum(sum(isnan(O)))==0
    CA=CA(1:poi,:);
    O=O(1:poi,:);
    if sec==2, sheetflag=1; else sheetflag=0; end
    [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,sheetflag);
    if ~isempty(backbone)
        switch sec
            case 0
                col=graph_settings.loop_color;
                [x,y,z]=tubeplot(backbone',graph_settings.coil_radius/2);
            case 1
                col=get_color(sec,TM_flags(n));
                [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec);
           case 2
                col=get_color(sec,TM_flags(n));
                [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec,graph_settings.spr);
        end
        [mm,nn]=size(x);
        for kk=1:length(CA)-1
            rnum=resnum(kk);
            aa=(kk-1)*graph_settings.spr;
            ee=aa+graph_settings.spr;
            if aa<1, aa=1; end
            if ee>nn, ee=nn; end
            if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum),'residue_graphics')
                if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics,'objects')
                    if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects)
                        for kkk=1:length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects)
                            if ishandle(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects(kkk))
                                delete(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects(kkk));
                            end
                        end
                    end
                end
            end
            han=surf(x(:,aa:ee),y(:,aa:ee),z(:,aa:ee));
            set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
            set(han, 'CDataMapping','direct');
            minx=min(min(x(:,aa:ee)));
            maxx=max(max(x(:,aa:ee)));
            miny=min(min(y(:,aa:ee)));
            maxy=max(max(y(:,aa:ee)));
            minz=min(min(z(:,aa:ee)));
            maxz=max(max(z(:,aa:ee)));
            rindices=[indices(1:3) rnum];
            record_object(han,rindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            residue_graphics.mode=4;
            residue_graphics.objects=han;    
            residue_graphics.color=col;
            residue_graphics.opaque=falpha;
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics=residue_graphics;
        end
    end
end

graphics.mode=4;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics=graphics;

function col=get_color(sec,TM_flag)

global graph_settings

switch sec
    case 0
        col=graph_settings.coil_color;
    case 1
        if TM_flag
            col=graph_settings.TM_helix_color;
        else
            col=graph_settings.helix_color;
        end
    case 2
        if TM_flag
            col=graph_settings.TM_sheet_color;
        else
            col=graph_settings.sheet_color;
        end
end
