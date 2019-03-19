function gobjects=plot_backbone_segment(indices,mode)
% function gobjects=plot_backbone_segment(indices,mode)
%
% Plots a backbone segment for a particular residue and returns a
% vector of graphics objects handles
%
% the handle graphics objects and graphics settings are stored in a
% structure variable that is indexed in the same way as the atom
% coordinates, this structure residue_graphics has the fields
% .objects  handles to Matlab graphics objects
% .mode     residue graphics mode, 0 none, 1 CaWire, 2 CaStick, 3 coil, 4:
%           ribbon
% .color    RGB color vector
% .opaque   opaqueness (1-transparency) 
%
% indices   indices to identify structure, model, and chain
% atom      atom identifier in the model
% elnum     element number
% mode      plot mode
%
% gobjects  vector of handles to graphics objects of all alternate locations
%
% G. Jeschke, 2009

global graph_settings % default graphics settings
global model
global hMain

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;

gobjects=[];
residue_graphics={};

maxx=-1e6;
minx=1e6;
maxy=-1e6;
miny=1e6;
maxz=-1e6;
minz=1e6;

recorded=0;

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'residue_graphics'),
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;
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

opaque=1;
switch mode
    case 'CaWire'
        obj1=[];
        obj2=[];
        my_color=graph_settings.CaWire_color;
        [previous,this,next]=segment_CA_coor(indices);
        if ~isempty(this) && ~isempty(previous),
            link=this+(previous-this)/2;
            x=[this(1) link(1)];
            y=[this(2) link(2)];
            z=[this(3) link(3)];
            maxx=max([maxx,max(x)]);
            minx=min([minx,min(x)]);
            maxy=max([maxy,max(y)]);
            miny=min([miny,min(y)]);
            maxz=max([maxz,max(z)]);
            minz=min([minz,min(z)]);
            obj1=plot3(x,y,z,'LineWidth',graph_settings.CaWire_width,'Color',my_color);
            set(obj1,'Clipping','off');
        end;
        if ~isempty(this) && ~isempty(next),
            link=this+(next-this)/2;
            x1=[this(1) link(1)];
            y1=[this(2) link(2)];
            z1=[this(3) link(3)];
            obj2=plot3(x1,y1,z1,'LineWidth',graph_settings.CaWire_width,'Color',my_color);
            maxx=max([maxx,max(x1)]);
            minx=min([minx,min(x1)]);
            maxy=max([maxy,max(y1)]);
            miny=min([miny,min(y1)]);
            maxz=max([maxz,max(z1)]);
            minz=min([minz,min(z1)]);
            set(obj2,'Clipping','off');
        end;
        gobjects=[obj1 obj2];
        residue_graphics.mode=1;
    case 'CaStick'
        obj1=[]; obj2=[]; obj3=[];
        my_color=graph_settings.CaStick_color;
        [previous,this,next]=segment_CA_coor(indices);
        if ~isempty(this) && ~isempty(previous),
            link=this+(previous-this)/2;
            [x,y,z] = vec2tube(this,link,graph_settings.CaStick_radius,20);
            maxx=max([maxx,max(x)]);
            minx=min([minx,min(x)]);
            maxy=max([maxy,max(y)]);
            miny=min([miny,min(y)]);
            maxz=max([maxz,max(z)]);
            minz=min([minz,min(z)]);
            obj1=surf(x,y,z);
            set(obj1, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj1, 'CDataMapping','direct','AlphaDataMapping','none');
        end;
        if ~isempty(this) && ~isempty(next),
            link=this+(next-this)/2;
            [x,y,z] = vec2tube(this,link,graph_settings.CaStick_radius,20);
            maxx=max([maxx,max(x)]);
            minx=min([minx,min(x)]);
            maxy=max([maxy,max(y)]);
            miny=min([miny,min(y)]);
            maxz=max([maxz,max(z)]);
            minz=min([minz,min(z)]);
            obj2=surf(x,y,z);
            set(obj2, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj2, 'CDataMapping','direct','AlphaDataMapping','none');
        end;
        if ~isempty(this),
            [x,y,z]=point2sphere(this,graph_settings.CaStick_radius);
            maxx=max([maxx,max(x)]);
            minx=min([minx,min(x)]);
            maxy=max([maxy,max(y)]);
            miny=min([miny,min(y)]);
            maxz=max([maxz,max(z)]);
            minz=min([minz,min(z)]);
            obj3=surf(x,y,z);
            set(obj3, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj3, 'CDataMapping','direct','AlphaDataMapping','none');
        end;
        gobjects=[obj1 obj2 obj3];
        residue_graphics.mode=2;
    case 'coil'
        type=1;
        if isfield(model.structures{indices(1)}(indices(2)),'seqtype'), % determine sequence type, information is present
            type=model.structures{indices(1)}(indices(2)).seqtype;
        end;
        switch type
            case 1
                mk_secondary_ribbon(indices,1);
                recorded=1;
            case 2
                mk_NA_ribbon(indices,1);
                recorded=1;
        end;
        residue_graphics.mode=3;
     case 'ribbon'
        type=1;
        if isfield(model.structures{indices(1)}(indices(2)),'seqtype'), % determine sequence type, information is present
            type=model.structures{indices(1)}(indices(2)).seqtype;
        end;
        switch type
            case 1
                mk_secondary_ribbon(indices);
                recorded=1;
            case 2
                mk_NA_ribbon(indices);
                recorded=1;
            otherwise
                recorded=1;
        end;
        residue_graphics.mode=4;
    otherwise % default is CaWire
        my_color=graph_settings.CaWire_color;
        [previous,this,next]=segment_CA_coor(indices);
        if ~isempty(previous) && ~isempty(this) && ~isempty(next),
            x=[this(1) previous(1)];
            y=[this(2) previous(2)];
            z=[this(3) previous(3)];
            obj1=plot3(x,y,z,'LineWidth',graph_settings.CaWire_width,'Color',my_color);
            x1=[this(1) next(1)];
            y1=[this(2) next(2)];
            z1=[this(3) next(3)];
            maxx=max([max(x),max(x1)]);
            minx=min([min(x),min(x1)]);
            maxy=max([max(y),max(y1)]);
            miny=min([min(y),min(y1)]);
            maxz=max([max(z),max(z1)]);
            minz=min([min(z),min(z1)]);
            obj2=plot3(x,y,z,'LineWidth',graph_settings.CaWire_width,'Color',my_color);
            gobjects=[obj1 obj2];
        end;
        residue_graphics.mode=1;
end;
if ~recorded,
    for k=1:length(gobjects),
        if ~isempty(gobjects(k)),
            if gobjects(k)~=0,
                record_object(gobjects(k),indices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            end;
        end;
    end;
    residue_graphics.objects=gobjects;    
    residue_graphics.color=my_color;
    residue_graphics.opaque=opaque;
    model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=residue_graphics;
end;


function [previous,this,next]=segment_CA_coor(indices)
% Returns the Cartesian coordinates of the C_alpha atoms of the previous
% residue, the current (this) residue, and the next residue
% the current residue is addressed by indices
% if the previous or next residue do not exist, the missing coordinate is
% the same as for the current residue

global model

monotype=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).type;

this=[];
previous=[];
next=[];

if monotype==1,
    this=atom_coor(indices,'CA');

    if indices(4)>1 
        pindices=indices;
        pindices(4)=pindices(4)-1;
        previous=atom_coor(pindices,'CA');
    end;
    if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'terminal')
        if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).terminal), return; end;
        pindices=indices;
        pindices(4)=pindices(4)+1;
        next=atom_coor(pindices,'CA');
    end;
end;
if monotype==2,
    this1=atom_coor(indices,'C4''');
    this2=atom_coor(indices,'C3''');
    if ~isempty(this1) && ~isempty(this2),
        this=(this1+this2)/2;
    else
        this=[];
    end;
    if indices(4)>1 
        pindices=indices;
        pindices(4)=pindices(4)-1;
        previous1=atom_coor(pindices,'C4''');
        previous2=atom_coor(pindices,'C3''');
        if ~isempty(previous1) && ~isempty(previous2),
            previous=(previous1+previous2)/2;
        else
            previous=[];
        end;
    end;
    if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'terminal')
        if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).terminal), return; end;
        pindices=indices;
        pindices(4)=pindices(4)+1;
        next1=atom_coor(pindices,'C4''');
        next2=atom_coor(pindices,'C3''');
        if ~isempty(next1) && ~isempty(next2),
            next=(next1+next2)/2;
        else
            next=[];
        end;
    end;
end;

function error=mk_secondary_ribbon(indices,coil)
% 
% Plots and records a backbone segment of a tiled ribbon or coil for a
% single residue
%
% indices   residue indices (1...4)
% coil      optional flag that requests a coil segment (if 1) rather than 
%           a ribbon segment (if 0), defaults to 0
% error     flag that indicates when nothing was plotted (1), 0 if
%           successful

global model
global graph_settings

if nargin<2,
    coil=0;
end;

error=1;

falpha=1;

% Determine range of residues with the same secondary structure
sec=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).secondary; % secondary structure type
range=indices(4);
secn=sec;
poi=indices(4);
while poi>1 && secn==sec,
    poi=poi-1;
    secn=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(poi).secondary;
    if secn==sec,
        range=[poi range];
    end;
end;
if indices(4)>1,
    range=[poi range];
    extend=1;
else
    extend=0;
end;
secn=sec;
poi=indices(4);
stopflag=0;
while poi<length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info) && secn==sec,
    poi=poi+1;
    secn=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(poi).secondary;
    if secn==sec,
        range=[range poi];
    end;
    if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'terminal')
        if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).terminal),
            stopflag=1;
            break; 
        end; 
    end;
end;

if ~stopflag,
    range=[range poi];
    endextend=1;
else
    endextend=0;
end;

n=length(range);
pindices=indices;

mynum=[];

limit=4;
if sec==2, 
    limit=6;
end;
CA=zeros(2*(limit-1)+1,3);
O=zeros(2*(limit-1)+1,3);

poi=0;
for k=1:n,
    if abs(range(k)-indices(4))<limit,
        info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(range(k));
        CA_id=tag2id('CA',info.atom_tags);
        O_id=tag2id('O',info.atom_tags);
        secend=1;
        if ~isempty(CA_id) && ~isempty(O_id),
            poi=poi+1;
            if range(k)==indices(4), mynum=poi; end;
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
            resnum(poi)=range(k);
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
    end;
end;
if extend,
    CA(2,:)=(CA(1,:)+CA(2,:))/2;
    O(2,:)=(O(1,:)+O(2,:))/2;
    CA=CA(2:poi,:);
    O=O(2:poi,:);
    poi=poi-1;
    mynum=mynum-1;
else
    CA=CA(1:poi,:);
    O=O(1:poi,:);
end;

if poi>1,
    CA(poi,:)=(CA(poi-1,:)+CA(poi,:))/2;
    O(poi,:)=(O(poi-1,:)+O(poi,:))/2;
end;

if coil, sec=4; end;
if sec==2, sheetflag=1; else sheetflag=0; end;

if sec>0 && sec<4,
    [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA,O,sheetflag);
else
    [backbone,rung,normal]=mk_ribbon(graph_settings.spr,CA);
end;

if ~isempty(backbone)
    switch sec
        case 0
            col=graph_settings.loop_color;
            [x,y,z]=tubeplot(backbone',graph_settings.coil_radius/2);
            residue_graphics.mode=4;
        case 1
            col=graph_settings.helix_color;
            [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec);
            residue_graphics.mode=4;
        case 2
            col=graph_settings.sheet_color;
            [x,y,z]=generalized_tubeplot(backbone',rung',normal',sec,graph_settings.spr);
            residue_graphics.mode=4;
        case 4
            col=graph_settings.coil_color;
            [x,y,z]=tubeplot(backbone',graph_settings.coil_radius);
            residue_graphics.mode=3;
    end;
    [mm,nn]=size(x);
    rnum=indices(4);
    if ~isempty(mynum),
        aa=(mynum-1)*graph_settings.spr;
        ee=aa+graph_settings.spr+2;
        if aa<1, aa=1; end;
        if ee>nn, ee=nn; end;
        if ee>aa,
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
            han=surf(x(:,aa:ee),y(:,aa:ee),z(:,aa:ee));
            set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
            set(han, 'CDataMapping','direct','AlphaDataMapping','none');
            minx=min(min(x(:,aa:ee)));
            maxx=max(max(x(:,aa:ee)));
            miny=min(min(y(:,aa:ee)));
            maxy=max(max(y(:,aa:ee)));
            minz=min(min(z(:,aa:ee)));
            maxz=max(max(z(:,aa:ee)));
            record_object(han,indices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            residue_graphics.objects=han;    
            residue_graphics.color=col;
            residue_graphics.opaque=falpha;
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(rnum).residue_graphics=residue_graphics;
            error=0;
        end;
    end;
end;

function error=mk_NA_ribbon(indices,coil)
% 
% Plots and records a backbone segment of a tiled ribbon or coil for a
% single nucleotide in a DNA/RNA
%
% indices   residue indices (1...4)
% coil      optional flag that requests a coil segment (if 1) rather than 
%           a ribbon segment (if 0), defaults to 0
% error     flag that indicates when nothing was plotted (1), 0 if
%           successful

global model
global graph_settings

if nargin<2,
    coil=0;
end;

error=1;

falpha=1;

% Determine range of residues with the same secondary structure
range=indices(4);
poi=indices(4);
while poi>1,
    poi=poi-1;
    range=[poi range];
end;
poi=indices(4);
stopflag=0;
while poi<length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info),
    poi=poi+1;
    range=[range poi];
end;

n=length(range);
pindices=indices;

mynum=[];

limit=6;

C1=zeros(2*(limit-1)+1,3);
C3=zeros(2*(limit-1)+1,3);
C4=zeros(2*(limit-1)+1,3);

poi=0;
for k=1:n,
    if abs(range(k)-indices(4))<limit,
        info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(range(k));
        C4_id=tag2id('C4''',info.atom_tags);
        C3_id=tag2id('C3''',info.atom_tags);
        C1_id=tag2id('C1''',info.atom_tags);
        if ~isempty(C4_id) && ~isempty(C3_id) && ~isempty(C1_id),
            poi=poi+1;
            resnum(poi)=range(k);
            if range(k)==indices(4), mynum=poi; end;
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
        end;
    end;
end;

C4=C4(1:poi,:);
C3=C3(1:poi,:);
C1=C1(1:poi,:);

[backbone,rung,normal]=mk_ribbon(graph_settings.spr,(C4+C3)/2,C1);

col=graph_settings.NA_color;

if ~isempty(backbone)
    [x,y,z]=generalized_tubeplot(backbone',rung',normal',3,graph_settings.NA_spr);
    if coil,
        residue_graphics.mode=3;
    else
        residue_graphics.mode=4;
    end;
    [mm,nn]=size(x);
    rnum=indices(4);
    if ~isempty(mynum),
        aa=(mynum-1)*graph_settings.spr;
        ee=aa+graph_settings.spr+2;
        if aa<1, aa=1; end;
        if ee>nn, ee=nn; end;
        if ee>aa,
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
            han=surf(x(:,aa:ee),y(:,aa:ee),z(:,aa:ee));
            set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
            set(han, 'CDataMapping','direct','AlphaDataMapping','none');
            minx=min(min(x(:,aa:ee)));
            maxx=max(max(x(:,aa:ee)));
            miny=min(min(y(:,aa:ee)));
            maxy=max(max(y(:,aa:ee)));
            minz=min(min(z(:,aa:ee)));
            maxz=max(max(z(:,aa:ee)));
            record_object(han,indices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            if ~coil,
                [nuc_coor,ncol]=get_nucleotide_coor(indices);
                if ~isempty(nuc_coor),
                    back_point=sum(nuc_coor(1:2,:),1)/2;
                    base_edge=nuc_coor(3,:);
                    set(han,'FaceColor',ncol);
                    han2=capped_tube(back_point,base_edge,ncol,graph_settings.NA_width/4,falpha);
                    set(han2,'EdgeColor','none','FaceLighting','gouraud','Clipping','off','CDataMapping','direct');
                    for khan=1:length(han2),
                        record_object(han2(khan),indices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
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
                        record_object(han3(khan),indices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
                    end;
                else
                    han3=[];
                end;
            else
                han2=[];
                han3=[];
                ncol=col;
            end;
            residue_graphics.objects=[han han2 han3];    
            residue_graphics.color=ncol;
            residue_graphics.opaque=falpha;
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=residue_graphics;
            error=0;
        end;
    end;
end;