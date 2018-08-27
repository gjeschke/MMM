function [message,argout]=set_chain(indices,property,argin)
% function [message,argout]=set_chain(indices,property,argin)
%
% Sets properties of a chain (all models of that) in MMM
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
        message=color_chain(indices,argin{1});
    case 'colorscheme'
        message=colorscheme_chain(indices,argin);
    case 'hide'
        message=hide_chain(indices);
    case 'show'
        message=show_chain(indices,argin{1});
    case 'transform'
        message=transform_chain(indices,argin);
    case 'transparency'
        message=transparency_chain(indices,argin{1});
    case 'undefine'
        undefine_secondary(indices(1),indices(2));
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function message=color_chain(indices,color)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

message.error=0;
message.text='';

[msg,info_text]=get_chain(indices,'info');
chain_model_address=[info_text{1} '{:}'];
set_object(chain_model_address,'color',{color});

function message=colorscheme_chain(indices,scheme)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

message.error=0;
message.text='';

[msg,info_text]=get_chain(indices,'info');
chain_model_address=[info_text{1} '{:}'];
set_object(chain_model_address,'colorscheme',scheme);

function message=transform_chain(indices,matrices)
% Coordinate transformation defined by an affine 4x4 transformation matrix
% or a cell array of such matrices

[message,newindices]=get_chain(indices,'children');
[m,n]=size(newindices);
for k=1:m,
    message=set_chain_model(newindices(k,:),'transform',matrices);
end;

function message=transparency_chain(indices,alpha)

message.error=0;
message.text='';

[msg,info_text]=get_chain(indices,'info');
chain_model_address=[info_text{1} '{:}'];
set_object(chain_model_address,'transparency',{alpha});

function message=hide_chain(indices)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

num_models=length(model.structures{indices(1)}(indices(2)).residues);

if ~isempty(num_models) && num_models>=1,
    for k=1:num_models,
        pindices=[indices(1:2) k];
        [message,argout]=set_chain_model(pindices,'hide');
    end;
end;

function message=show_chain(indices,mode)
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

num_models=length(model.structures{indices(1)}(indices(2)).residues);

if ~isempty(num_models) && num_models>=1,
    for k=1:num_models,
        pindices=[indices(1:2) k];
        [message,argout]=set_chain_model(pindices,'show',my_mode);
        if message.error==0,
            plotted=1;
        end;
    end;
end;

if plotted==0,
    message.error=2;
    message.text='Nothing to plot.';
end;

function undefine_secondary(snum,cnum)

global model

model.structures{snum}(cnum).loop_defs={};
model.structures{snum}(cnum).helix_defs={};
model.structures{snum}(cnum).sheet_defs={};
