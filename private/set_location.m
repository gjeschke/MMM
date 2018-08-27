function [message,argout]=set_location(indices,property,argin)
% function [message,argout]=set_location(indices,property,argin)
%
% Sets properties of an atom location in MMM
%
% indices   index vector that identifies the atom location in the model
% property  property to be set, given as a string, available properties and
%           corresponding set functions are defined at the beginning of the
%           function source code in the switch board
% argin     further input arguments given as a cell array argin{:}
%
% message   error message structure with fields .error and .text,
%           .error=0 indicates no error
% argout    further output arguments given as structure argout{n}
%
% G. Jeschke, 2009

argout={};

message.error=0;
message.text='';

if nargin<3,
    argin{1}='';
end;

switch property
    case 'annotations'
	message=set_annotations(indices,argin{1});
    case 'color'
        message=color_location(indices,argin{1});
    case 'colorscheme'
        message=colorscheme_location(indices,argin{1});
    case 'hide'
        message=hide_location(indices);
    case 'show'
        message=show_location(indices,argin{1});
    case 'transform'
        message=transform_location(indices,argin);
    case 'transparency'
        message=transparency_location(indices,argin{1});
    case 'uncolor'
        message=uncolor_location(indices);
    case 'untransparency'
        message=untransparency_location(indices);
    case 'xyz'
        message=xyz_location(indices,argin{1});
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function message=color_location(indices,rgb)
% Sets color of atom graphics
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(indices(6),1);
graphics=[];
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics'),
    if length(model.structures{indices(1)}(indices(2)).atom_graphics)>=indices(3),
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom,
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
        end;
    end;
end;

if isempty(graphics), return; end;
if isempty(graphics.objects), return; end;


if ~isempty(graphics)
    if ~isempty(graphics.objects),
        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                if isprop(graphics.objects(k),'Color'),
                    set(graphics.objects(k),'Color',rgb);
                elseif isprop(graphics.objects(k),'FaceColor'),
                    set(graphics.objects(k),'FaceColor',rgb);
                end;
            end;
        end;
        graphics.color=[rgb;graphics.color];
    end;
end;

model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom)=graphics;

function message=uncolor_location(indices)
% Resets color of an atom location to previous value
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(indices(6),1);
graphics=[];
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics'),
    if length(model.structures{indices(1)}(indices(2)).atom_graphics)>=indices(3),
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom,
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
        end;
    end;
end;
if isempty(graphics), return; end;
if isempty(graphics.objects), return; end;

if graphics.mode==1,
    prop='Color';
else
    prop='FaceColor';
end;

[m,n]=size(graphics.color);
if m<2,
    message.error=3;
    message.text='No previous color available';
    return
end;

graphics.color=graphics.color(2:m,:);
rgb=graphics.color(1,:);

for k=1:length(graphics.objects),
    if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
        set(graphics.objects(k),prop,rgb);
    end;
end;

model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom)=graphics;

function message=colorscheme_location(indices,scheme)
% Sets color of atom graphics for a location according to color scheme
%

global model
global graph_settings

message.error=0;
message.text='';

known=0; % flag for known scheme

switch scheme
    case 'chain'
        chains=length(model.structures{indices(1)});
        rgb=color_grade(indices(2),chains);
        known=1;
    case 'secondary'
        sec_type=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).secondary;
        rgb=graph_settings.loop_color;
        if sec_type==1,
            rgb=graph_settings.helix_color;
        end;
        if sec_type==2,
            rgb=graph_settings.sheet_color;
        end;
        known=1;
    case 'sequence'
        residues=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);
        rgb=color_grade(indices(4),residues);
        known=1;
    case 'Bfactor'
        [message,Bfactors]=get_location(indices,'Bfactor');
        Bfactor=sqrt(mean(Bfactors));
        rgb=color_grade(Bfactor-graph_settings.Bmin+1,graph_settings.Bmax-graph_settings.Bmin+1);
        known=1;
    case 'Bfactor_tight'
        [message,Bfactors]=get_location(indices,'Bfactor');
        Bfactor=sqrt(mean(Bfactors));
        Bmin=round(sqrt(model.info{indices(1)}.B_range(1)));
        Bmax=round(sqrt(model.info{indices(1)}.B_range(2)));
        rgb=color_grade(Bfactor-Bmin+1,Bmax-Bmin+1);
        known=1;
    case 'charge'
        num=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
        if num>0 && num<=length(model.structures{indices(1)}(indices(2)).sequence),
            slc=model.structures{indices(1)}(indices(2)).sequence(num);
        else
            slc='?';
        end;
        rgb=color_by_charge(slc);
        known=1;
end;

if known,
    color_location(indices,rgb);
end;


function message=transform_location(indices,matrices)
% Coordinate transformation defined by an affine 4x4 transformation matrix
% or a cell array of such matrices

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
xyz0=coors(atoms(indices(6),1),:);

xyz1=affine_trafo_point(xyz0,matrices);
coors(atoms(indices(6),1),:)=xyz1;
model.structures{indices(1)}(indices(2)).xyz{indices(3)}=coors;


function message=transparency_location(indices,alpha)
% Sets transparency of atom location graphics
%

global model

message.error=0;
message.text='';

atom0=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom0(indices(6),1);

if alpha<0, % alpha<0 selects transparency by location
    if length(atom0(indices(6),:))>1,
        alpha=atom0(indices(6),2);
    else
        alpha=1;
    end;
end;

graphics=[];
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics') ...
        && ~isempty(model.structures{indices(1)}(indices(2)).atom_graphics)...
        && length(model.structures{indices(1)}(indices(2)).atom_graphics)>=indices(3),
    if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom,
        graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
    end;
end;
if isempty(graphics), return; end;
if isempty(graphics.objects), return; end;

if graphics.mode==1,
    message.error=3;
    message.text='No transparency for line objects.';
    return
end;

if ~isempty(graphics)
    if ~isempty(graphics.objects),
        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                if isprop(graphics.objects(k),'FaceAlpha'),
                    set(graphics.objects(k),'FaceAlpha',alpha);
                end;
            end;
        end;
        graphics.opaque=[alpha;graphics.opaque];
    end;
end;

model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom)=graphics;

function message=untransparency_location(indices)
% Resets transparency of an atom location to previous value
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(indices(6),1);
graphics=[];
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics'),
    if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom,
        graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
    end;
end;

if isempty(graphics), return; end;
if isempty(graphics.objects), return; end;

if graphics.mode==1,
    message.error=3;
    message.text='No transparency for line objects.';
    return
end;

m=length(graphics.opaque);
if m<2,
    message.error=3;
    message.text='No previous transparency value available';
    return
end;

graphics.opaque=graphics.opaque(2:m,:);
alpha=graphics.opaque(1);

for k=1:length(graphics.objects),
    if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
        set(graphics.objects(k),'FaceAlpha',alpha);
    end;
end;

model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom)=graphics;

function message=hide_location(indices)
% Hides an atom location by deleting graphics objects
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(indices(6),1);

graphics=[];
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
    if length(model.structures{indices(1)}(indices(2)).atom_graphics)>=indices(3),
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
            model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom).objects=[];
            model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom).mode=0;
        end;
    end;
end;

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

    
function message=show_location(indices,mode)
% Plots an atom location by calling plot routine
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

if isempty(mode),
    mode='space-filling';
end;

elnum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).elements(indices(5));
gobjects=plot_atom(indices,elnum,mode);

if isempty(gobjects),
    message.error=1;
    message.text='Nothing to plot.';
end;
    
function message = xyz_location(indices,ncoor)
% sets new xyz coordinates for a location
%
% indices  index vector, length at least 6


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
coors(atoms(indices(6),1),:) = ncoor;
model.structures{indices(1)}(indices(2)).xyz{indices(3)} = coors;