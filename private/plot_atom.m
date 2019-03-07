function gobjects=plot_atom(indices,elnum,mode)
% function gobjects=plot_atom(indices,elnum,mode)
%
% Plots an atom and the bonds originating from that atom and returns a
% vector of graphics objects handles
%
% the handle graphics objects and graphics settings are stored in a
% structure variable that is indexed in the same way as the atom
% coordinates, this structure atom_graphics has the fields
% .objects  handles to Matlab graphics objects
% .mode     atom graphics mode, 0 none, 1 wire, 2 stick, 3 ball&stick, 4:
%           space-filling
% .color    RGB color vector
% .opaque   opaqueness (1-transparency) 
%
% indices   indices to identify structure, chain, model, residue, and atom
% elnum     element number
% mode      plot mode
%
% gobjects  vector of handles to graphics objects of all alternate locations
%
% G. Jeschke, 2009

global graph_settings % default graphics settings
global chemistry % pse information (colors, vdW radii)
global model
global hMain

if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;

if ~isfield(hMain,'atom_graphics') || ~ishandle(hMain.atom_graphics),
    hMain.atom_graphics=hgtransform;
end;
if ~isfield(hMain,'atom_graphics_reduced') || ~ishandle(hMain.atom_graphics_reduced),    
    hMain.atom_graphics_reduced=hgtransform;
end;

switch hMain.atom_graphics_mode,
    case 0
        full_vis='off';
        red_vis='off';
    case 1
        full_vis='on';
        red_vis='off';
    case 2
        full_vis='off';
        red_vis='on';
end;
maxbonds=graph_settings.maxbonds; 

gobjects=[];

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};


[m,n]=size(atom); % determine whether there are alternate locations (n>1)
if n==1,
    atom=[atom 1.0]; % default transparency is unity
end;

coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
my_color=chemistry.pse(elnum).rgb;
my_vdW=chemistry.pse(elnum).vdW;

if length(indices)>5,
    kstart=indices(6);
    kend=indices(6);
else
    kstart=1;
    kend=m;
end;

for k=kstart:kend, % loop over all alternate locations
    indices(6)=k;
    if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
        if length(model.structures{indices(1)}(indices(2)).atom_graphics)>=indices(3)
            mmm=length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(:));
            if mmm>=atom(k,1),
                graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
                if ~isempty(graphics),
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
    coor=coors(atom(k,1),:);
    opaque=atom(k,2);
    obj=plot3(coor(1),coor(2),coor(3),'.','Color',my_color,'Visible',red_vis,'Parent',hMain.atom_graphics_reduced);
    aobjects=obj;
    switch mode
        case 'wire'
            x=coor(1); y=coor(2); z=coor(3);
            obj=plot3(x,y,z,'.','Color',my_color,'Parent',hMain.atom_graphics,'Visible',full_vis);
            set(obj,'Clipping','off');
            xyz=[x,y,z,0.1,0.1,0.1];
            record_object(obj,indices,xyz);
            aobjects=[aobjects obj];
            atom_graphics.mode=1;
        case {'stick','cofactor'}
%             [x,y,z]=point2sphere(coor,graph_settings.stick_radius);
%             obj=surface(x,y,z);
            [x,y,z,t]=point2trisphere(coor,graph_settings.stick_radius);
            obj=trisurf(t,x,y,z);
            % reducepatch(obj);
            xyz=[coor(1),coor(2),coor(3),2*graph_settings.stick_radius,2*graph_settings.stick_radius,2*graph_settings.stick_radius];
            record_object(obj,indices,xyz);
            aobjects=[aobjects obj];
            set(obj, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            atom_graphics.mode=2;
        case 'ball&stick'
%             [x,y,z]=point2sphere(coor,graph_settings.ball_radius);
%             obj=surface(x,y,z);
            [x,y,z,t]=point2trisphere(coor,graph_settings.ball_radius);
            obj=trisurf(t,x,y,z);
            xyz=[coor(1),coor(2),coor(3),2*graph_settings.ball_radius,2*graph_settings.ball_radius,2*graph_settings.ball_radius];
            record_object(obj,indices,xyz);
            aobjects=[aobjects obj];
            set(obj, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            atom_graphics.mode=3;
        case 'space-filling'
%             [x,y,z]=point2sphere(coor,my_vdW);
%             obj=surface(x,y,z);
            [x,y,z,t]=point2trisphere(coor,my_vdW,2);
            obj=trisurf(t,x,y,z);
            aobjects=[aobjects obj];
            set(obj, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            xyz=[mean(x),mean(y),mean(z),2*my_vdW,2*my_vdW,2*my_vdW];
            record_object(obj,indices,xyz);
            atom_graphics.mode=4;
        case 'marker'            
            marker_color = [0.25,0.25,0.25];
            [x,y,z] = vec2tube(coor-[1,0,0],coor+[1,0,0],graph_settings.stick_radius,10);
            obj=surface(x,y,z);
            set(obj, 'FaceColor', marker_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            xyz=[min(min(x)),min(min(y)),min(min(z)),max(max(x)),max(max(y)),max(max(z))];
            record_object(obj,indices,xyz);
            [x,y,z] = vec2tube(coor-[0,0,1],coor+[0,0,1],graph_settings.stick_radius,10);
            obj=surface(x,y,z);
            set(obj, 'FaceColor', marker_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            xyz=[min(min(x)),min(min(y)),min(min(z)),max(max(x)),max(max(y)),max(max(z))];
            record_object(obj,indices,xyz);
            [x,y,z] = vec2tube(coor-[0,1,0],coor+[0,1,0],graph_settings.stick_radius,10);
            obj=surface(x,y,z);
            set(obj, 'FaceColor', marker_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            xyz=[min(min(x)),min(min(y)),min(min(z)),max(max(x)),max(max(y)),max(max(z))];
            record_object(obj,indices,xyz);
            atom_graphics.mode=5;
        otherwise % default is stick
%             [x,y,z]=point2sphere(coor,graph_settings.stick_radius);
%             obj=surface(x,y,z);
            [x,y,z,t]=point2trisphere(coor,graph_settings.stick_radius);
            obj=trisurf(t,x,y,z);
            xyz=[mean(mean(x)),mean(mean(y)),mean(mean(z)),max(max(x))-min(min(x)),max(max(y))-min(min(y)),max(max(z))-min(min(z))];
            record_object(obj,indices,xyz);
            aobjects=[aobjects obj];
            set(obj, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
            atom_graphics.mode=2;
    end;
    bonded=model.structures{indices(1)}(indices(2)).conn(atom(k,1),:);
    bonds=zeros(maxbonds,3);
    poi=0;
    for kk=1:length(bonded),
        if bonded(kk)>0, % there is a bonded atom
            poi=poi+1;
            bond_end=coor+(coors(bonded(kk),:)-coor)/2;
            xm=(coor(1)+bond_end(1))/2;
            ym=(coor(2)+bond_end(2))/2;
            zm=(coor(3)+bond_end(3))/2;
            xd=max([coor(1) bond_end(1)])-min([coor(1) bond_end(1)]);
            yd=max([coor(2) bond_end(2)])-min([coor(2) bond_end(2)]);
            zd=max([coor(3) bond_end(3)])-min([coor(3) bond_end(3)]);
            switch  mode
                case 'stick'
                    atom_graphics.mode=2;
                case 'cofactor'
                    atom_graphics.mode=2;
                case 'ball&stick'
                    atom_graphics.mode=3;
                case 'space-filling'
                    atom_graphics.mode=4;
            end;
            switch mode
                case 'wire'
                    x=[coor(1),bond_end(1)];
                    y=[coor(2),bond_end(2)];
                    z=[coor(3),bond_end(3)];
                    obj=plot3(x,y,z,'-','LineWidth',graph_settings.wire_width,'Color',my_color,'Parent',hMain.atom_graphics,'Visible',full_vis);
                    set(obj,'Clipping','off');
                    xyz=[xm,ym,zm,xd,yd,zd];
                    record_object(obj,indices,xyz);
                    aobjects=[aobjects obj];
                    atom_graphics.mode=1;
                case {'stick','ball&stick','cofactor'}
                    [x,y,z] = vec2tube(coor,bond_end,graph_settings.stick_radius,10);
                    obj=surface(x,y,z);
                    xyz=[xm,ym,zm,xd,yd,zd];
                    record_object(obj,indices,xyz);
                    aobjects=[aobjects obj];
                    set(obj, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
                    set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
                otherwise
                    [x,y,z] = vec2tube(coor,bond_end,graph_settings.stick_radius,10);
                    obj=surface(x,y,z);
                    xyz=[xm,ym,zm,xd,yd,zd];
                    record_object(obj,indices,xyz);
                    aobjects=[aobjects obj];
                    set(obj, 'FaceColor', my_color, 'EdgeColor', 'none', 'FaceAlpha',opaque,'FaceLighting','gouraud','Clipping','off');
                    set(obj, 'CDataMapping','direct','AlphaDataMapping','none','Parent',hMain.atom_graphics,'Visible',full_vis);
                    atom_graphics.mode=2;
            end;
        end;
    end;
    atom_graphics.objects=aobjects;    
    atom_graphics.color=my_color;
    atom_graphics.opaque=opaque;
    gobjects=[gobjects aobjects]; % collect handles of all graphics objects
    model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=atom_graphics;
end;