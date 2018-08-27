function plot_model
%

global hMain
global model

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;


snum=length(model.structures);
for ks=1:snum,
    cnum=length(model.structures{ks});
    for kc=1:cnum,
        modnum=length(model.structures{ks}(kc).residues);
        for km=1:modnum,
            if isfield(model.structures{ks}(kc).residues{km},'chain_model_graphics')
                graphics=model.structures{ks}(kc).residues{km}.chain_model_graphics;
                if ~isempty(graphics),
                    if ~isempty(graphics.objects),
                        for kk=1:length(graphics.objects),
                            if graphics.objects(kk)~=0,
                                model.unrecord=[model.unrecord graphics.objects(kk)];
                            end;
                        end;
                    end;
                    switch graphics.mode
                        case 3
                            set_object([ks,kc,km],'show',{'string'});
                            [mcol,ncol]=size(graphics.color);
                            if mcol>1,
                                set_object([ks,kc,km],'color',{graphics.color(1,:),0});
                            end;
                    end;
                end;
            end;
            if isfield(model.structures{ks}(kc).residues{km},'secondary_graphics')
                allgraphics=model.structures{ks}(kc).residues{km}.secondary_graphics;
                if ~isempty(allgraphics),
                    for k=1:length(allgraphics),
                        graphics=allgraphics(k);
                        if ~isempty(graphics),
                            if ~isempty(graphics.objects),
                                for kk=1:length(graphics.objects),
                                    if graphics.objects(kk)~=0,
                                        model.unrecord=[model.unrecord graphics.objects(kk)];
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
            rnum=length(model.structures{ks}(kc).residues{km}.info);
            for kr=1:rnum, % residue graphics
                graphics=model.structures{ks}(kc).residues{km}.info(kr).residue_graphics;
                model.structures{ks}(kc).residues{km}.info(kr).residue_graphics.objects={};
                if ~isempty(graphics),
                    if ~isempty(graphics.objects),
                        for k=1:length(graphics.objects),
                            if graphics.objects(k)~=0,
                                model.unrecord=[model.unrecord graphics.objects(k)];
                            end;
                        end;
                    end;
                    if isfield(graphics,'mode'),
                        switch graphics.mode
                            case 1
                                set_object([ks,kc,km,kr],'show',{'CaWire'});
                            case 2
                                set_object([ks,kc,km,kr],'show',{'CaStick'});
                            case 3
                                set_object([ks,kc,km,kr],'show',{'coil'});
                            case 4
                                set_object([ks,kc,km,kr],'show',{'ribbon'});
                            case 5
                                set_object([ks,kc,km,kr],'show',{'cartoon'});
                        end;
                        if graphics.mode>0,
                            [mcol,ncol]=size(graphics.color);
                            if mcol>1,
                                set_object([ks,kc,km,kr],'color',{graphics.color(1,:)});
                            end;
                        end;
                    end;
                end;
            end;
            for kr=1:rnum, % residue label graphics
                graphics=model.structures{ks}(kc).residues{km}.info(kr).label_graphics;
                model.structures{ks}(kc).residues{km}.info(kr).label_graphics.objects=[];
                if ~isempty(graphics)
                    if ~isempty(graphics.objects),
                        for k=1:length(graphics.objects),
                            if graphics.objects(k)~=0,
                                model.unrecord=[model.unrecord graphics.objects(k)];
                            end;
                        end;
                        switch graphics.mode
                            case 1
                                set_object([ks,kc,km,kr],'show',{'label'});
                            case 2
                                set_object([ks,kc,km,kr],'show',{'label_frame'});
                        end;
                        if graphics.mode>0,
                            [mcol,ncol]=size(graphics.color);
                            if mcol>1,
                                set_object([ks,kc,km,kr],'color',{graphics.color(1,:)});
                            end;
                        end;
                    end;
                end;
            end;
            for kr=1:rnum, % atom graphics
                anum=length(model.structures{ks}(kc).residues{km}.info(kr).atom_numbers);
                for ka=1:anum,
                    atoms=model.structures{ks}(kc).residues{km}.info(kr).atom_numbers{ka};
                    [lnum,n]=size(atoms);
                    if isfield(model.structures{ks}(kc),'atom_graphics')
                        if length(model.structures{ks}(kc).atom_graphics)>=km
                            mmm=length(model.structures{ks}(kc).atom_graphics{km}(:));
                            for kl=1:lnum,
                                if atoms(kl,1)<=mmm,
                                    graphics=model.structures{ks}(kc).atom_graphics{km}(atoms(kl,1));
                                    if ~isempty(graphics),
                                        if ~isempty(graphics.objects),
                                            for kk=1:length(graphics.objects),
                                                if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                                                    model.unrecord=[model.unrecord graphics.objects(kk)];
                                                end;
                                            end;      
                                        end;
                                    end;
                                end;
                            end;
                        end;
                        if ~isempty(model.structures{ks}(kc).atom_graphics),
                            mmm=length(model.structures{ks}(kc).atom_graphics{km}(:));
                            for kl=1:lnum,
                                if atoms(kl,1)<=mmm,
                                    atom_graphics=model.structures{ks}(kc).atom_graphics{km}(atoms(kl,1));
                                    if ~isempty(atom_graphics.mode)
                                        switch atom_graphics.mode
                                            case 1
                                                set_object([ks,kc,km,kr,ka,kl],'show',{'wire'});
                                            case 2
                                                set_object([ks,kc,km,kr,ka,kl],'show',{'stick'});
                                            case 3
                                                set_object([ks,kc,km,kr,ka,kl],'show',{'ball&stick'});
                                            case 4
                                                set_object([ks,kc,km,kr,ka,kl],'show',{'space-filling'});
                                        end;
                                        if atom_graphics.mode>1,
                                            [mcol,ncol]=size(atom_graphics.color);
                                            if mcol>1,
                                                set_object([ks,kc,km,kr,ka,kl],'color',{atom_graphics.color(1,:)});
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    if isfield(model,'info') && isfield(model.info{snum},'bilayer'),
        if ~isempty(model.info{snum}.bilayer),
            model.info{snum}.bilayer=make_planes(model.info{snum}.bilayer);
        end;
    end;
end;

set(gcf,'Renderer','opengl');
lighting gouraud 
material shiny

if isfield(model,'surfaces'),
    for k=1:length(model.surfaces),
        [mode,tag]=strtok(model.surfaces(k).tag,':');
        switch mode
            case 'density'
                id=tag2id(tag(2:end),model.density_tags);
                plot_density(k,id,model.surfaces(k).active);
            case 'SAS'
                plot_SAS(k,tag(2:end),model.surfaces(k).active);
        end;
    end;
end;

if isfield(model,'motion'),
    for k=1:length(model.motion),
        graphics=[];
        color=model.motion(k).color;
        alpha=model.motion(k).transparency;
        radius=model.motion(k).radius;
        coor0=model.motion(k).template;
        coor1=model.motion(k).target;
        if model.motion(k).active,
            [m0,n0]=size(coor0);
            [m1,n1]=size(coor1);
            if ~isempty(coor0) && m0==m1,
                for kca=1:m0,
                    obj=arrow(coor0(kca,:),coor1(kca,:),color,radius,alpha(kca));
                    graphics=[graphics obj];
                end;
            end;
        end;
        model.motion(k).gobjects=graphics;
    end;
end;

highlight_selection;

function plot_density(k,id,active)

global model

dg=model.surfaces(k);

level=dg.level;

cube=model.densities{id}.cube;
cube=cube/max(max(max(cube)));
x=model.densities{id}.x;
y=model.densities{id}.y;
z=model.densities{id}.z;

[xg,yg,zg]=meshgrid(x,y,z);
p = patch(isosurface(xg,yg,zg,cube,level));
set(p, 'FaceColor', dg.color, 'EdgeColor', 'none','FaceAlpha',dg.transparency,'FaceLighting','gouraud','Clipping','off');
set(p, 'CDataMapping','direct','AlphaDataMapping','none');
if ~active,
    set(p,'Visible','off');
end;
dg.gobjects=p;
model.surfaces(k)=dg;

function plot_SAS(k,tag,active)

global model

dg=model.surfaces(k);
for kk=1:length(model.SAS),
    if strcmp(tag,model.SAS(kk).tag),
        obj = trisurf(model.SAS(kk).tri,model.SAS(kk).x,model.SAS(kk).y,model.SAS(kk).z);
        set(obj, 'FaceColor', dg.color, 'EdgeColor', 'none', 'FaceAlpha',dg.transparency,'FaceLighting','gouraud','Clipping','off');
        set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
        dg.gobjects=obj;
        if ~active,
            set(obj,'Visible','off');
        end;
    end;
end;
model.surfaces(k)=dg;

function bilayer=make_planes(bilayer)

global graph_settings

x=[-bilayer.width/2 bilayer.width/2 bilayer.width/2 -bilayer.width/2];
y=[-bilayer.width/2 -bilayer.width/2 bilayer.width/2 bilayer.width/2];
z=zeros(1,4);

hc=patch(x,y,z,graph_settings.bilayer_midplane_color);
hb=patch(x,y,z-bilayer.thickness/2*ones(1,4),graph_settings.bilayer_outer_color);
ht=patch(x,y,z+bilayer.thickness/2*ones(1,4),graph_settings.bilayer_outer_color);

set(hc,'FaceAlpha',0.25);
set(hb,'FaceAlpha',0.25);
set(ht,'FaceAlpha',0.25);
bilayer.graphics=[hc hb ht];

if bilayer.show_center,
    set(hc,'Visible','on');
else
    set(hc,'Visible','off');
end;

if bilayer.show_outer,
    set(hb,'Visible','on');
    set(ht,'Visible','on');
else
    set(hb,'Visible','off');
    set(ht,'Visible','off');
end;
