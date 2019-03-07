function [message,argout]=set_atom(indices,property,argin)
% function [message,argout]=set_atom(indices,property,argin)
%
% Sets properties of a residue in MMM
%
% indices   index vector that identifies the atom in the model
% property  property to be set, given as a string, available properties and
%           corresponding set functions are defined at the beginning of the
%           function source code in the switch board
% argin     further input arguments given as a structure argin{:}
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
        message=color_atom(indices,argin{1});
    case 'colorscheme'
        message=colorscheme_atom(indices,argin{1});
    case 'hide'
        message=hide_atom(indices);
    case 'show'
        message=show_atom(indices,argin{1});
    case 'transform'
        message=transform_atom(indices,argin);
    case 'transparency'
        message=transparency_atom(indices,argin{1});
    case 'uncolor'
        message=uncolor_atom(indices);
    case 'untransparency'
        message=untransparency_atom(indices);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function message=color_atom(indices,rgb)
% Sets color of atom graphics
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
[m,n]=size(atom); % determine whether there are alternate locations (n>1)
if n==1,
    atom=[atom 1.0]; % default transparency is unity
end;

kstart=1;
kend=m;

for k=kstart:kend,
    graphics=[];
    if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
        end;
    end;
    if isempty(graphics), continue; end;
    if isempty(graphics.objects), continue; end;
    for kk=1:length(graphics.objects),
        if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
            if isprop(graphics.objects(kk),'Color'),
                set(graphics.objects(kk),'Color',rgb);
            elseif isprop(graphics.objects(kk),'FaceColor'),
                set(graphics.objects(kk),'FaceColor',rgb);
            end;
        end;
    end;
    graphics.color=[rgb;graphics.color];
    model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=graphics;
end;

function message=uncolor_atom(indices)
% Resets color of atom graphics to previous one
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
[m,n]=size(atom); % determine whether there are alternate locations (n>1)
if n==1,
    atom=[atom 1.0]; % default transparency is unity
end;

kstart=1;
kend=m;

for k=kstart:kend,
    graphics=[];
    if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
        end;
    end;
    if isempty(graphics), continue; end;
    if isempty(graphics.objects), continue; end;

    if graphics.mode==1,
        prop='Color';
    else
        prop='FaceColor';
    end;

    [m,n]=size(graphics.color);
    if m<2,
        message.error=3;
        message.text='No previous color available';
        continue
    end;
    
    graphics.color=graphics.color(2:m,:);
    rgb=graphics.color(1,:);

    for kk=1:length(graphics.objects),
        if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
            set(graphics.objects(kk),prop,rgb);
        end;
    end;
    model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=graphics;
end;

function message=colorscheme_atom(indices,scheme)
% Sets color of atom graphics according to color scheme
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
        [message,Bfactors]=get_atom(indices,'Bfactor');
        Bfactor=sqrt(mean(Bfactors));
        rgb=color_grade(Bfactor-graph_settings.Bmin+1,graph_settings.Bmax-graph_settings.Bmin+1);
        known=1;
    case 'Bfactor_tight'
        [message,Bfactors]=get_atom(indices,'Bfactor');
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
    color_atom(indices,rgb);
end;



function message=transform_atom(indices,matrices)
% Coordinate transformation defined by an affine 4x4 transformation matrix
% or a cell array of such matrices

[message,newindices]=get_atom(indices,'children');
[m,n]=size(newindices);
for k=1:m,
    message=set_location(newindices(k,:),'transform',matrices);
end;

function message=transparency_atom(indices,alpha)
% Sets transparency of atom graphics
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
[m,n]=size(atom); % determine whether there are alternate locations (n>1)
if n==1,
    atom=[atom 1.0]; % default transparency is unity
end;

kstart=1;
kend=m;

for k=kstart:kend,
    if alpha<0, % alpha<0 selects transparency by location
        ralpha=atom(k,2);
    else
        ralpha=alpha;
    end;
    graphics=[];
    if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
        end;
    end;
    if isempty(graphics), continue; end;
    if isempty(graphics.objects), continue; end;

    if graphics.mode==1,
        message.error=3;
        message.text='No transparency for line objects.';
        continue
    end;

    for kk=1:length(graphics.objects),
        if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
            try
                set(graphics.objects(kk),'FaceAlpha',ralpha);
            catch
            end;
        end;
    end;
    graphics.opaque=[ralpha;graphics.opaque];
    model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=graphics;
end;


function message=untransparency_atom(indices)
% Resets transparency of atom graphics to previous one
%

global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
[m,n]=size(atom); % determine whether there are alternate locations (n>1)
if n==1,
    atom=[atom 1.0]; % default transparency is unity
end;

kstart=1;
kend=m;

for k=kstart:kend,
    graphics=[];
    if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
        end;
    end;
    if isempty(graphics), continue; end;
    if isempty(graphics.objects), continue; end;

    if graphics.mode==1,
        message.error=3;
        message.text='No transparency for line objects.';
        continue
    end;

    m=length(graphics.opaque);
    if m<2,
        message.error=3;
        message.text='No previous transparency value available';
        continue
    end;
    
    graphics.opaque=graphics.opaque(2:m);
    alpha=graphics.opaque(1);

    for kk=1:length(graphics.objects),
        if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
            set(graphics.objects(kk),'FaceAlpha',alpha);
        end;
    end;
    model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=graphics;
end;

function message=hide_atom(indices)
% Hides an atom by deleting graphics objects
%

global model
global hMain

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(1,1);
graphics=[];
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
    if ~isempty(model.structures{indices(1)}(indices(2)).atom_graphics)
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
            end;
        end;
    end;
end;

function message=show_atom(indices,mode)
% Plots an atom by calling plot routines for all locations
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
    
