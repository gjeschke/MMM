function [message,argout]=get_location(indices,property,varargin)
% function [message,argout]=get_location(indices,property,varargin)
%
% Gets properties of an atom location in MMM
%
% indices   index vector that identifies the atom location in the model
% property  property to be set, given as a string, available properties and
%           corresponding set functions are defined at the beginning of the
%           function source code in the switch board
% varargin  further input arguments given as a cell array varargin{:}
%
% message   error message structure with fields .error and .text,
%           .error=0 indicates no error
% argout    further output arguments given as structure argout{n}
%
% G. Jeschke, 2009

argout=[];
if nargin>=3,
    argin=varargin{1};
else
    argin='';
end;

message.error=0;
message.text='';

switch property
    case 'annotations'
        [message,argout]=get_annotations(indices);
    case 'Bfactor'
        [message,argout]=Bfactor_location(indices);
    case {'children','descendants'}
        argout=indices;
    case 'graphics'
        [message,argout]=graphics_location(indices);
    case 'info'
        [message,argout]=info_location(indices);
    case 'element'
        [message,argout]=element_location(indices);
    case 'mass'
        [message,argout]=mass_location(indices);
    case 'population'
        [message,argout]=population_location(indices);
    case 'xyz'
        [message,argout]=xyz_location(indices);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function [message,info_text]=info_location(indices)
% Plots an atom location by calling plot routine
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

stag=id2tag(indices(1),model.structure_tags,model.structure_ids);
ctag=id2tag(indices(2),model.chain_tags{indices(1)},model.chain_ids{indices(1)});
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
rtag=id2tag(indices(4),model.structures{indices(1)}(indices(2)).residues{indices(3)}.residue_tags);
atag=id2tag(indices(5),info.atom_tags);
adr_string=sprintf('[%s](%s){%i}%s.%s',stag,ctag,indices(3),rtag,atag);
if isfield(info,'location_tags')
    ltag=id2tag(indices(6),info.location_tags);
    if isempty(ltag), 
        ltag='-'; 
    else
        adr_string=sprintf('%s:%s',adr_string,ltag);
    end;
else
    ltag='-';
end;
info_text={adr_string sprintf('Structure: %s',stag) sprintf('Chain      : %s',ctag) sprintf('Model      : %i',indices(3))...
    sprintf('Residue  : %s %s',info.name,rtag) sprintf('Atom       : %s',atag) sprintf('Location  : %s',ltag)};

function [message,graphics]=graphics_location(indices)
% returns graphics handles for a location
%
% indices  index vector, length at least 6


global model

message.error=0;
message.text='';

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(indices(6),1);

graphics=[];

if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics'),
    if length((model.structures{indices(1)}(indices(2)).atom_graphics))>=indices(3),
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom,
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
        end;
    end;
end;

function [message,element]=element_location(indices)
% returns element (atomic number) for a location
%
% indices  index vector, length at least 6


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
element=isotopes(atoms(indices(6),1),1);

function [message,mass]=mass_location(indices)
% returns mass for a location
%
% indices  index vector, length at least 6


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
mass=isotopes(atoms(indices(6),1),2);

function [message,xyz]=xyz_location(indices)
% returns xyz coordinates for a location
%
% indices  index vector, length at least 6


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
xyz=coors(atoms(indices(6),1),:);

function [message,Bfactor]=Bfactor_location(indices)
% returns Bfactor for a location
%
% indices  index vector, length at least 6

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
Bfactors=model.structures{indices(1)}(indices(2)).Bfactor{indices(3)};
Bfactor=Bfactors(atoms(indices(6),1));

function [message,pop]=population_location(indices)
% returns population for a location
%
% indices  index vector, length at least 6


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
pop=atoms(indices(6),2);

