function [message,argout]=set_structure(indices,property,argin)
% function [message,argout]=set_structure(indices,property,argin)
%
% Sets properties of a structure (all chains/molecules) in MMM
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

if nargin<3,
    argin{1}='';
end;

switch property
    case 'annotations'
	message=set_annotations(indices,argin{1});
    case 'color'
        message=color_structure(indices,argin{1});
    case 'colorscheme'
        message=colorscheme_structure(indices,argin);
    case 'hide'
        message=hide_structure(indices);
    case 'show'
        message=show_structure(indices,argin{1});
    case 'transform'
        message=transform_structure(indices,argin);
    case 'transparency'
        message=transparency_structure(indices,argin{1});
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function message=color_structure(indices,color)

message.error=0;
message.text='';

[msg,info_text]=get_structure(indices,'info');
chain_address=[info_text{1} '(:)'];
set_object(chain_address,'color',{color});

function message=colorscheme_structure(indices,argin)

message.error=0;
message.text='';

[msg,info_text]=get_structure(indices,'info');
chain_address=[info_text{1} '(:)'];
set_object(chain_address,'colorscheme',argin);

function message=transform_structure(indices,matrices)
% Coordinate transformation defined by an affine 4x4 transformation matrix
% or a cell array of such matrices

[message,newindices]=get_structure(indices,'children');
[m,n]=size(newindices);
for k=1:m,
    message=set_chain(newindices(k,:),'transform',matrices);
end;

function message=transparency_structure(indices,alpha)

message.error=0;
message.text='';

[msg,info_text]=get_structure(indices,'info');
chain_address=[info_text{1} '(:)'];
set_object(chain_address,'transparency',{alpha});


function message=hide_structure(indices)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

num_chains=length(model.structures{indices(1)});

if ~isempty(num_chains) && num_chains>=1,
    for k=1:num_chains,
        pindices=[indices(1) k];
        [message,argout]=set_chain(pindices,'hide');
    end;
    if isfield(model,'motion') && ~isempty(model.motion),
        stag=id2tag(indices(1),model.structure_tags);
        for k=1:length(model.motion),
            if length(model.motion(k).tag)>7,
                tag=model.motion(k).tag(8:end);
                address=['$motion:' tag];
                stagc=strtok(tag,'_');
                if strcmp(stagc,stag),
                    message=set_surface(address,'delete');
                end;
            end;
        end;
    end;
end;

function message=show_structure(indices,mode)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

plotted=0;

if isempty(mode),
    mode='wire';
end;

my_mode{1}=mode;

num_chains=length(model.structures{indices(1)});

if ~isempty(num_chains) && num_chains>=1,
    for k=1:num_chains,
        pindices=[indices(1) k];
        [message,argout]=set_chain(pindices,'show',my_mode);
        if message.error==0,
            plotted=1;
        end;
    end;
end;

if plotted==0,
    message.error=2;
    message.text='Nothing to plot.';
end;
