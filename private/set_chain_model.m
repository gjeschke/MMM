function [message,argout]=set_chain_model(indices,property,argin)
% function [message,argout]=set_chain_model(indices,property,argin)
%
% Sets properties of one model (coordinate set) of a chain in MMM
%
% indices   index vector that identifies the residue in the model
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

if nargin<3
    argin{1}='';
end

if length(argin) < 2
    argin{2}=[];
end


switch property
    case 'annotations'
	message=set_annotations(indices,argin{1});
    case 'color'
        message=color_chain_model(indices,argin);
    case 'colorscheme'
        message=colorscheme_chain_model(indices,argin);
    case 'hide'
        message=hide_chain_model(indices);
    case 'show'
        message=show_chain_model(indices,argin{1},argin{2});
    case 'transform'
        message=transform_chain_model(indices,argin);
    case 'transparency'
        message=transparency_chain_model(indices,argin{1});
    case 'update'
        message=update_chain_model(indices);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end


function message=update_chain_model(indices)
% Updates a ribbon model of a chain after sceondary structure changes

global model

message.error=0;
message.text='';


if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'secondary_graphics'),
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics;
    if isfield(graphics,'mode'),
        if graphics.mode==4;
            show_chain_model(indices,'ribbon');
        end;
    end;
end;

function message=color_chain_model(indices,argin)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

color=argin{1};

if length(argin)>1 && ~isempty(argin{2})
    child_flag=argin{2};
else
    child_flag=1;
end;

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'chain_model_graphics')
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics;
    if ~isempty(graphics),
        if ~isempty(graphics.objects),
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    set(graphics.objects(k),'FaceColor',color);
                end;
            end;
            graphics.color=[color;graphics.color];
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics=graphics;
        end;
    end;
end;

if child_flag, % residues should have same color
    [msg,info_text]=get_chain_model(indices,'info');
    residue_address=[info_text{1} ':'];
    set_object(residue_address,'color',{color});
end;

function message=colorscheme_chain_model(indices,argin)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

scheme=argin{1};

if length(argin)>2,
    child_flag=argin{3};
else
    child_flag=1;
end;

switch scheme
    case 'chain'
        chains=length(model.structures{indices(1)});
        color=color_grade(indices(2),chains);
        if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'chain_model_graphics')
            graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics;
            if ~isempty(graphics),
                if ~isempty(graphics.objects),
                    for k=1:length(graphics.objects),
                        if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                            set(graphics.objects(k),'FaceColor',color);
                        end;
                    end;
                    graphics.color=[color;graphics.color];
                    model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics=graphics;
                end;
            end;
        end;
end;

if child_flag, % residues should have same color
    [msg,info_text]=get_chain_model(indices,'info');
    residue_address=[info_text{1} ':'];
    set_object(residue_address,'colorscheme',argin);
end;

function message=transparency_chain_model(indices,alpha)

message.error=0;
message.text='';

if alpha<0, alpha=1; end; % default transparency for chain models is opaque

[msg,info_text]=get_chain_model(indices,'info');
residue_address=[info_text{1} ':'];
set_object(residue_address,'transparency',{alpha});

function message=hide_chain_model(indices)
% Deletes residue graphics of residue children and atom graphics of atom children 
%

global model

message.error=0;
message.text='';

message.error=0;
message.text='';

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'chain_model_graphics')
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics;
    if ~isempty(graphics),
        if ~isempty(graphics.objects),
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    delete(graphics.objects(k));
                end;
            end;
        end;
    end;
    model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics={};
end;

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'secondary_graphics')
    allgraphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics;
    if ~isempty(allgraphics),
        for k=1:length(allgraphics),
            graphics=allgraphics(k);
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
    model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics={};
end;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;

if ~isempty(info) && length(info)>=1,
    % Loop over all residues
    for k=1:length(info),
        pindices=[indices(1:3) k];
        [message,argout]=set_residue(pindices,'hide');
    end;
end;

function message=show_chain_model(indices,mode,radius)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model
global graph_settings

if exist('radius','var') && ~isempty(radius)
    default_coil_radius = graph_settings.coil_radius;
    graph_settings.coil_radius = str2double(radius);
end

message.error=0;
message.text='';

plotted=0;

if isempty(mode),
    mode='wire';
end;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;

if strcmp(mode,'string')
    type=1;
    if isfield(model.structures{indices(1)}(indices(2)),'seqtype'), % determine sequence type, information is present
        type=model.structures{indices(1)}(indices(2)).seqtype;
    end;
    switch type
        case 1
            objects=plot_full_coil(indices);
            plotted=1;
        case 2
            objects=plot_full_NA_coil(indices);
            plotted=1;
    end;
elseif strcmp(mode,'ribbon')
    type=1;
    if isfield(model.structures{indices(1)}(indices(2)),'seqtype'), % determine sequence type, information is present
        type=model.structures{indices(1)}(indices(2)).seqtype;
    end;
    switch type
        case 1
            objects=plot_tiled_ribbon(indices);
            plotted=1;
        case 2
            objects=plot_tiled_NA_ribbon(indices);
            plotted=1;
    end;
elseif strcmp(mode,'coil')
    type=1;
    if isfield(model.structures{indices(1)}(indices(2)),'seqtype'), % determine sequence type, information is present
        type=model.structures{indices(1)}(indices(2)).seqtype;
    end;
    switch type
        case 1
            objects=plot_tiled_coil(indices);
            plotted=1;
        case 2
            objects=plot_tiled_NA_ribbon(indices,true);
            plotted=1;
    end;
end;        
if ~isempty(info) && length(info)>=1 && ~plotted,
    % Loop over all residues
    for k=1:length(info),
        my_mode{1}='wire';
        pindices=[indices(1:3) k];
        switch info(k).type
            case 1 % amino acid
                switch mode
                    case {'CaStick','CaWire','wire','stick','ball&stick','spacefilling','label'} % cofactor to be removed on implementation of proper nucleic acid plotting
                        my_mode{1}=mode;
                    otherwise
                        my_mode{1}=[];
                end;
            case 2 % nuclei acid
                switch mode
                    case {'CaStick','wire','CaWire'} % cofactor to be removed on implementation of proper nucleic acid plotting
                        my_mode{1}=mode;
                    case {'ribbon','coil'} % cofactorwire to be removed on implementation of proper nucleic acid plotting
                        my_mode{1}='wire';
                    otherwise
                        my_mode{1}=mode;
                end;
            otherwise % atomistic plot
                switch mode
                    case {'CaStick','cofactor'}
                        my_mode{1}='stick';
                    case {'wire','CaWire','cofactorwire'}
                        my_mode{1}='wire';
                    otherwise
                        my_mode{1}=mode;
                end;
        end;
        if ~isempty(my_mode{1}),
            [message,argout]=set_residue(pindices,'show',my_mode);
        end;
        if message.error==0,
            plotted=1;
        end;
    end;
end;

if plotted==0,
    message.error=2;
    message.text='Nothing to plot.';
end;

if exist('radius','var') && ~isempty(radius)
    graph_settings.coil_radius = default_coil_radius;
end


function message=transform_chain_model(indices,matrices)
% Coordinate transformation defined by an affine 4x4 transformation matrix
% or a cell array of such matrices

[message,newindices]=get_chain_model(indices,'children');
[m,n]=size(newindices);
for k=1:m,
    message=set_residue(newindices(k,:),'transform',matrices);
end;
