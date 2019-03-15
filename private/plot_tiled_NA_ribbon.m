function gobjects=plot_tiled_NA_ribbon(indices,coil)
% Plots a tiled ribbon for DNA or RNA addressed by structure, chain and
% coordinate set indices and assigns the graphics objects to residues
%
% indices   first three indices into structure
% coil      optional flag, if 1 a coil (only backbone ribbon) is displayed,
%           otherwise tubes to base and base polygons are also shown
%
% gobjects  graphics objects
%
% G. Jeschke, 2010

global model
global graph_settings
global hMain

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;

if nargin<2,
    coil=false;
end;

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
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

graphics.objects=[];

gobjects=[];
colors=zeros(1000,3);
secondary=zeros(1,1000);
range=zeros(1000,2);
minmax=zeros(1000,6);

if isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info),
    return;
end;

n=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);

C1=zeros(n,3);
C3=zeros(n,3);
C4=zeros(n,3);
resnum=zeros(1,n);
poi=0;
cpoi=0;
starter=1;
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
    if isempty(C3_id),
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
        [mm,nn]=size(x);
        for kk=1:length(C3),
            rnum=resnum(kk);
            cindices=[indices rnum];
            aa=(kk-1)*graph_settings.spr;
            ee=aa+graph_settings.spr;
            if aa<1, aa=1; end;
            if ee>nn, ee=nn; end;
            if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum),'residue_graphics'),
                if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics,'objects'),
                    if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects)
                        for kkk=1:length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects),
                            if ishandle(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects(kkk)),
                                delete(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics.objects(kkk));
                            end;
                        end;
                    end;
                end;
            end;
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
            if ~coil,
                [nuc_coor,ncol]=get_nucleotide_coor(cindices);
                if ~isempty(nuc_coor),
                    back_point=sum(nuc_coor(1:2,:),1)/2;
                    base_edge=nuc_coor(3,:);
                    set(han,'FaceColor',ncol);
                    han2=capped_tube(back_point,base_edge,ncol,graph_settings.NA_width/4,falpha);
                    set(han2,'EdgeColor','none','FaceLighting','gouraud','Clipping','off','CDataMapping','direct');
                    for khan=1:length(han2),
                        record_object(han2(khan),cindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
                    end;
                else
                    han2=[];
                    ncol=col;
                end;
                [mnc,nnc]=size(nuc_coor);
                if mnc>=5,
                    polygon=nuc_coor(3:end,:);
                    han3=fill3(polygon(:,1),polygon(:,2),polygon(:,3),ncol);
                    set(han3,'EdgeColor','none','FaceColor',ncol,'FaceLighting','gouraud','Clipping','off','CDataMapping','direct');
                    for khan=1:length(han3),
                        record_object(han3(khan),cindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
                    end;
                else
                    han3=[];
                end;
            else
                han2=[];
                han3=[];
                ncol=col;
            end;
            residue_graphics.mode=4;
            residue_graphics.objects=[han han2 han3];    
            residue_graphics.color=ncol;
            residue_graphics.opaque=falpha;
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics=residue_graphics;
        end;
    end;
end;

graphics.mode=4;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics=graphics;
