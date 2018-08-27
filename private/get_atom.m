function [message,argout]=get_atom(indices,property,varargin)
% function [message,argout]=get_atom(indices,property,varargin)
%
% Gets properties of an atom in MMM
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
        [message,argout]=Bfactor_atom(indices);
    case 'numbers'
        [message,argout]=number_atom(indices);
    case 'coor'
        [message,argout]=coor_atom(indices);
    case {'children','descendants'}
        [message,argout]=children_atom(indices);
    case 'graphics'
        [message,argout]=graphics_atom(indices);
    case 'info'
        [message,argout]=info_atom(indices);
    case 'element'
        [message,argout]=element_atom(indices);
    case 'mass'
        [message,argout]=mass_atom(indices);
    case 'populations'
        [message,argout]=pop_atom(indices);
    case 'xyz'
        [message,argout]=xyz_atom(indices);
end;

function [message,newindices]=children_atom(indices)
% returns indices of atom children
%
% indices  index vector, length must be 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};

[m,n]=size(atoms); % determine whether there are alternate locations (n>1)
if m>1,
    newindices=zeros(m,6);
    for k=1:m,
        newindices(k,1:5)=indices;
        newindices(k,6)=k;
    end;
else
    newindices=[indices 1];
end;

function [message,info_text]=info_atom(indices)
% returns info string for an atom
%
% indices  index vector, length at least 5

global model

message.error=0;
message.text='';

stag=id2tag(indices(1),model.structure_tags,model.structure_ids);
ctag=id2tag(indices(2),model.chain_tags{indices(1)},model.chain_ids{indices(1)});
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
rtag=id2tag(indices(4),model.structures{indices(1)}(indices(2)).residues{indices(3)}.residue_tags);
codelist=zeros(1,length(info.atom_numbers));
for k=1:length(codelist), codelist(k)=info.atom_numbers{k}(1,1); end;
atag=id2tag(indices(5),info.atom_tags);
adr_string=sprintf('[%s](%s){%i}%s.%s',stag,ctag,indices(3),rtag,atag);
info_text={adr_string sprintf('Structure: %s',stag) sprintf('Chain      : %s',ctag) sprintf('Model      : %i',indices(3))...
    sprintf('Residue  : %s %s',info.name,rtag) sprintf('Atom       : %s',atag)};


function [message,graphics]=graphics_atom(indices)
% returns graphics handles for an atom
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

graphics=[];

atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
atom=atom(1,1);
if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics'),
    if length((model.structures{indices(1)}(indices(2)).atom_graphics))>=indices(3),
        if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom,
            graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom);
        end;
    end;
end;

function [message,element]=element_atom(indices)
% returns element (atomic number) for an atom
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
element=isotopes(atoms(1,1),1);

function [message,mass]=mass_atom(indices)
% returns mass for an atom
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
mass=isotopes(atoms(1,1),2);


function [message,pop]=pop_atom(indices)
% returns populations of alternate atom locations, 1 is returned if the
% atom location is unique
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};

[m,n]=size(atoms); % determine whether there are alternate locations (n>1)
if n>1,
    pop=atoms(:,2);
else
    pop=1;
end;


function [message,xyz]=xyz_atom(indices)
% returns xyz coordinates for an atom
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};

[m,n]=size(atoms); % determine whether there are alternate locations (n>1)
xyz=zeros(m,3);

for k=1:m,
    xyz(k,:)=coors(atoms(k,1),:);
end;

function [message,coor]=coor_atom(indices)
% returns mean coordinate for an atom, averaged over all locations
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};

[m,n]=size(atoms); % determine whether there are alternate locations (n>1)
if n==1, atoms=[atoms 1]; end;
xyz=zeros(m,3);

popsum=0;
for k=1:m,
    xyz(k,:)=coors(atoms(k,1),:)*atoms(k,2);
    popsum=popsum+atoms(k,2);
end;

coor=sum(xyz,1)/popsum;

function [message,atom_numbers]=number_atom(indices)
% returns atom numbers (indices in the coordinate and isotope arrays)
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atom_numbers=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};

function [message,Bfactor]=Bfactor_atom(indices)
% returns xyz coordinates for an atom
%
% indices  index vector, length at least 5


global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{indices(5)};
Bfactors=model.structures{indices(1)}(indices(2)).Bfactor{indices(3)};

[m,n]=size(atoms); % determine whether there are alternate locations (n>1)
Bfactor=zeros(1,m);

for k=1:m,
    Bfactor(k)=Bfactors(atoms(k,1));
end;