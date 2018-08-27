function msg=set_surface(address,property,argin)
% Set a property of a surface specified by its address
% the surface address must begin with $, otherwise an empty object array is
% returned
% an empty object array is also returned, if a surface with the suppliedd
% address does not exist
% surface addresses have the format $type:tag
% where type is the surface type (currently density or SAS) and tag the
% surface tag
%
% G. Jeschke, 2009

global model

msg.error=0;
msg.text='No error.';

if (~isfield(model,'surfaces') || isempty(model.surfaces)) && (~isfield(model,'motion') || isempty(model.motion)), 
    msg.error=1;
    msg.text='ERROR: No surfaces or motion arrows defined in this model.';
    return; 
end; % return immediately if there are no surfaces

if ~strcmp(address(1),'$') || length(address)<2, % return immediately if address is invalid
    msg.error=2;
    msg.text='ERROR: Surface or motion arrow address must begin with $.';
    return
else
    mytag=address(2:end);
    found=0;
    if isfield(model,'surfaces') && ~isempty(model.surfaces)
        for k=1:length(model.surfaces),
            if strcmp(mytag,model.surfaces(k).tag),
                found=1;
                id=k;
                type='surface';
                item=model.surfaces(k);
                break
            end;
        end;
    end;
    if isfield(model,'motion') && ~isempty(model.motion)
        for k=1:length(model.motion),
            if strcmp(mytag,model.motion(k).tag),
                found=1;
                id=k;
                type='motion';
                item=model.motion(k);
                break
            end;
        end;
    end;
end;

if ~found,
    msg.error=3;
    msg.text=sprintf('ERROR: Surface or set of motion arrows with address "%s" does not exist.',address);
    return
end;

switch property
    case 'all'
        item=argin;
    case 'graphics'
        item.gobjects=argin;
    case 'color'
        % check argument
        proper=0;
        if ischar(argin) && numel(argin)==1,
            proper=1;
        elseif isa(argin,'double') && numel(argin)==3 && min(argin)>=0 && max(argin)<=1,
            proper=1;
        end;
        if proper,
            item.color=argin;
            set(item.gobjects,'FaceColor',argin);
        else
            msg.error=6;
            msg.text{1}='ERROR: Color argument must be either a single letter';
            msg.text{2}='       or an RGB vector with three numbers between 0 and 1.';
        end;
    case 'delete'
        obj=item.gobjects;
        if ~isempty(obj),
            for k=1:length(obj),
                if ishandle(obj),
                    delete(obj);
                end;
            end;
        end;
        % delete the underlying data
        msg=delete_data(item);
        % and the surface itself
        if isfield(model,'surfaces') && ~isempty(model.surfaces)
            surf=[];
            for k=1:length(model.surfaces),
                if ~strcmp(model.surfaces(k).tag,item.tag),
                    surf=[surf model.surfaces(k)];
                end;
            end;
            model.surfaces=surf;
        end;
        if isfield(model,'motion') && ~isempty(model.motion)
            motion=[];
            for k=1:length(model.motion),
                if ~strcmp(model.motion(k).tag,item.tag),
                    motion=[motion model.motion(k)];
                end;
            end;
            model.motion=motion;
        end;
        % the following would require a redraw method
        %     case 'level'
        %         model.surfaces(id).level=argin;
    case 'transparency'
        if isa(argin,'double') && numel(argin)==1 && argin>=0 && argin<=1,
            item.transparency=argin;
            set(item.gobjects,'FaceAlpha',argin);
        else
            msg.error=6;
            msg.text='ERROR: Transparency argument must be between 0 and 1.';
        end;
    case 'visible'
        if strcmp(argin,'off'),
            item.active=0;
        elseif strcmp(argin,'on')
            item.active=1;
        end;
        if strcmp(argin,'off') || strcmpi(argin,'on')
            set(item.gobjects,'Visible',argin);
        else
            msg.error=7;
            msg.text='ERROR: Visibility must be "on" or "off".';
        end;
    case 'show'
        if strcmpi(type,'surface'),
            item=show_surface(item);
        elseif strcmpi(type,'motion'),
            item=show_motion(item);
        end;
    otherwise
        msg.error=4;
        msg.text=sprintf('ERROR: Specified property "%s" does not exist.',property);
end;
if ~strcmp(property,'delete'),
    switch type
        case 'surface'
            model.surfaces(id)=item;
        case 'motion'
            model.motion(id)=item;
    end;
end;

function item=show_surface(item)
mode=strtok(item.tag,':');
switch mode
    case 'density'
        item=plot_density(item);
    case 'SAS'
        item=plot_SAS(item);
end;

function item=show_motion(item)
graphics=[];
color=item.color;
alpha=item.transparency;
radius=item.radius;
coor0=item.template;
coor1=item.target;
item.active=1;
[m0,n0]=size(coor0);
[m1,n1]=size(coor1);
if ~isempty(coor0),
    for kca=2:min([m0,m1])-1
        obj=arrow(coor0(kca,:),coor1(kca,:),color,radius,alpha(kca));
        graphics=[graphics obj];
    end;
end;
item.gobjects=graphics;

function msg=delete_data(surface)

global model

msg.error=0;
msg.text='No error.';

[type,tag]=strtok(surface.tag,':');
tag=tag(2:end);
switch type
    case 'density'
        if isfield(model,'densities') && ~isempty(model.densities),
            found=0;
            for k=1:length(model.densities),
                if strcmp(tag,model.densities{k}.tag),
                    found=1;
                    id=k;
                    break;
                end;
            end;
            if ~found,
                msg.error=12;
                msg.text='Warning: Underlying data of surface does not exist.';
            else
                tags=':';
                for k=1:length(model.densities),
                    ktag=id2tag(k,model.density_tags);
                    if k~=id,
                        tags=[tags ktag ':'];
                    end;
                end;
                model.density_tags=tags;
                if id==length(model.densities),
                    if id==1,
                        model.densities=[];
                    else
                        model.densities=model.densities{1:id-1};
                    end;
                else
                    if id==1,
                        model.densities=model.densities{2:end};
                    else
                        model.densities=[model.densities{1:id-1} model.densities{id+1:end}];
                    end;
                end;
            end;
        else
            msg.error=11;
            msg.text='Warning: No density data for this model.';
        end;
    case 'SAS'
        if isfield(model,'SAS') && ~isempty(model.SAS),
            found=0;
            for k=1:length(model.SAS),
                if strcmp(tag,model.SAS(k).tag),
                    found=1;
                    id=k;
                    break;
                end;
            end;
            if ~found,
                msg.error=12;
                msg.text='Warning: Underlying data of surface does not exist.';
            else
                if id==length(model.SAS),
                    if id==1,
                        model.SAS=[];
                    else
                        model.SAS=model.SAS(1:id-1);
                    end;
                else
                    if id==1,
                        model.SAS=model.SAS(2:end);
                    else
                        model.SAS=[model.SAS(1:id-1) model.SAS(id+1:end)];
                    end;
                end;
            end;
        else
            msg.error=13;
            msg.text='Warning: No solvent accessible surface data for this model.';
        end;
    case 'motion'
    otherwise
        msg.error=10;
        msg.text='Warning: Unknown surface type.';
end;

function dg=plot_density(dg)

global model

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
item.active=true;
dg.gobjects=p;

function dg=plot_SAS(dg)

global model

for kk=1:length(model.SAS),
    if strcmp(tag,model.SAS(kk).tag),
        obj = trisurf(model.SAS(kk).tri,model.SAS(kk).x,model.SAS(kk).y,model.SAS(kk).z);
        set(obj, 'FaceColor', dg.color, 'EdgeColor', 'none', 'FaceAlpha',dg.transparency,'FaceLighting','gouraud','Clipping','off');
        set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
        dg.gobjects=obj;
        item.active=true;
    end;
end
