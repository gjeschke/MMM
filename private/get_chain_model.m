function [message,argout]=get_chain_model(indices,property,varargin)
% function [message,argout]=get_chain_model(indices,property,varargin)
%
% Gets properties of a chain model in MMM
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


if nargin<3,
    argin={};
else
    argin=varargin{1};
end;
message.error=0;
message.text='';

switch property
    case 'annotations'
	[message,argout]=get_annotations(indices);
    case 'children'
        [message,argout]=children_chain_model(indices);
    case 'descendants'
        [message,argout]=descendants_chain_model(indices);
    case 'graphics'
        [message,argout]=graphics_chain_model(indices);
    case 'info'
        [message,argout]=info_chain_model(indices);
    case 'mass'
        [message,argout]=mass_chain_model(indices,argin);
    case 'secondary_graphics'
        [message,argout]=secondary_graphics_chain_model(indices);
    case 'xyz'
        [message,argout]=xyz_chain_model(indices,argin);
    case 'rcoor'
        [message,argout]=rcoor_chain_model(indices,argin);
    case 'xyz_heavy'
        [message,argout]=xyz_heavy_chain_model(indices,argin);
    case 'populations'
        [message,argout]=populations_chain_model(indices,argin);
    case 'elements'
        [message,argout]=elements_chain_model(indices,argin);
    case 'elements_heavy'
        [message,argout]=elements_heavy_chain_model(indices,argin);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function [message,newindices]=children_chain_model(indices)
% returns children indices for a chain model
%
% indices  index vector, length 3

global model

message.error=0;
message.text='';

residues=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
if isempty(residues)
    newindices=indices;
else
    newindices=zeros(length(residues),4);
    for k=1:length(residues)
        newindices(k,1:3)=indices;
        newindices(k,4)=k;
    end
end

function [message,newindices]=descendants_chain_model(indices)
% returns descendants indices for a chain model
%
% indices  index vector, length 3

message.error=0;
message.text='';

[message,mindices]=children_chain_model(indices);
[m,n]=size(mindices);
if n>3, % there are residues children, we have to check for children of their children
    newindices=zeros(10000,6);
    newindices(1:m,1:n)=mindices;
    poi=m;
    maxdim=4;
    for k=1:m,
        [message,nindices]=get_residue(mindices(k,:),'descendants');
        [mm,nn]=size(nindices);
        if nn>maxdim,
            maxdim=nn;
        end;
        newindices(poi+1:poi+mm,1:nn)=nindices;
        poi=poi+mm;
    end;
    newindices=newindices(1:poi,1:maxdim);
else
    newindices=indices;
end;

function [message,graphics]=graphics_chain_model(indices)
% returns graphics handles for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';
if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'chain_model_graphics'),
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.chain_model_graphics;
else
    graphics={};
end;

function [message,graphics]=secondary_graphics_chain_model(indices)
% returns graphics handles for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';
if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)},'secondary_graphics'),
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics;
else
    graphics={};
end;

function [message,info_text]=info_chain_model(indices)
% Plots an atom location by calling plot routine
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

stag=id2tag(indices(1),model.structure_tags,model.structure_ids);
ctag=id2tag(indices(2),model.chain_tags{indices(1)},model.chain_ids{indices(1)});
adr_string=sprintf('[%s](%s){%i}',stag,ctag,indices(3));
info_text={adr_string sprintf('Structure: %s',stag) sprintf('Chain      : %s',ctag) sprintf('Model      : %i',indices(3))};

function [message,mass]=mass_chain_model(indices,varargin)
% returns mass for a chain model
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';

if nargin>=2,
    argin=varargin{1};
else
    argin='';
end;

% mass=sum(model.structures{indices(1)}(indices(2)).isotopes(:,2)); % old
% version

if isempty(argin),
    waterflag=true;
else
    if strcmp(argin,'-nowater')
        waterflag=false;
    else
        waterflag=true;
    end;
end;

max_atoms=100000;
masses=zeros(max_atoms,1);

poi=0;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
ch_masses=model.structures{indices(1)}(indices(2)).isotopes(:,2);
for kk=1:length(info), % loop over residues
    if waterflag || ~strcmpi(info(kk).name,'HOH'),
        for kkk=1:length(info(kk).atom_numbers),
            anum=info(kk).atom_numbers{kkk};
            [ma,na]=size(anum);
            for k4=1:ma,
                poi=poi+1;
                masses(poi)=ch_masses(anum(k4,1));
            end;
        end;
    end;
end;
masses=masses(1:poi);
mass=sum(masses);


function [message,xyz]=xyz_chain_model(indices,varargin)
% returns xyz coordinates for a chain model
% water can be excluded by setting varargin{1}='-nowater';
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';

if nargin>=2,
    argin=varargin{1};
else
    argin='';
end;

% xyz=model.structures{indices(1)}(indices(2)).xyz{indices(3)}; % old
% version

if isempty(argin),
    waterflag=true;
else
    if strcmp(argin,'-nowater')
        waterflag=false;
    else
        waterflag=true;
    end;
end;

max_atoms=100000;
xyz=zeros(max_atoms,3);

poi=0;

if isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)})
    xyz = zeros(0,3);
    return
else
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
end
xyz_ch=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
for kk=1:length(info), % loop over residues
    if waterflag || ~strcmpi(info(kk).name,'HOH'),
        for kkk=1:length(info(kk).atom_numbers),
            anum=info(kk).atom_numbers{kkk};
            [ma,na]=size(anum);
            for k4=1:ma,
                poi=poi+1;
                xyz(poi,:)=xyz_ch(anum(k4,1),:);
            end;
        end;
    end;
end;
xyz=xyz(1:poi,:);

function [message,rcoor]=rcoor_chain_model(indices,varargin)
% returns xyz coordinates for a chain model
% water can be excluded by setting varargin{1}='-nowater';
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';

if nargin>=2,
    argin=varargin{1};
else
    argin='';
end;

% xyz=model.structures{indices(1)}(indices(2)).xyz{indices(3)}; % old
% version

if isempty(argin),
    waterflag=true;
else
    if strcmp(argin,'-nowater')
        waterflag=false;
    else
        waterflag=true;
    end;
end;

max_atoms=100000;
rcoor=zeros(max_atoms,4);

poi=0;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
xyz_ch=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
for kk=1:length(info), % loop over residues
    if waterflag || ~strcmpi(info(kk).name,'HOH'),
        for kkk=1:length(info(kk).atom_numbers),
            anum=info(kk).atom_numbers{kkk};
            [ma,na]=size(anum);
            for k4=1:ma,
                poi=poi+1;
                rcoor(poi,1)=kk;
                rcoor(poi,2:4)=xyz_ch(anum(k4,1),:);
            end;
        end;
    end;
end;
rcoor=rcoor(1:poi,:);

function [message,xyz]=xyz_heavy_chain_model(indices,varargin)
% returns xyz coordinates of only heavy (non hydrogen) atoms for a chain model
% only the first location of each atom is returned
% water is excluded
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';

if nargin>=2,
    argin=varargin{1};
else
    argin='';
end;

max_atoms=100000;
xyz=zeros(max_atoms,3);

poi=0;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
for kk=1:length(info), % loop over residues
    if ~strcmpi(info(kk).name,'HOH'),
        [message,coor]=get_residue([indices,kk],'xyz_heavy');
        if message.error ~= 0,
            xyz = [];
            return
        end;
        [mc,~] = size(coor);
        xyz(poi+1:poi+mc,:) = coor;
        poi = poi + mc;
    end;
end;
xyz=xyz(1:poi,:);

function [message,elements]=elements_chain_model(indices,varargin)
% returns atomic numbers of elements for a chain model
% water can be excluded by setting varargin{1}='-nowater';
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';

if nargin>=2,
    argin=varargin{1};
else
    argin='';
end;

% elements=model.structures{indices(1)}(indices(2)).isotopes(:,1); % old
% version

if isempty(argin),
    waterflag=true;
else
    if strcmp(argin,'-nowater')
        waterflag=false;
    else
        waterflag=true;
    end;
end;

max_atoms=100000;
elements=zeros(max_atoms,1);

poi=0;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
ch_elements=model.structures{indices(1)}(indices(2)).isotopes(:,1);
for kk=1:length(info), % loop over residues
    if waterflag || ~strcmpi(info(kk).name,'HOH'),
        for kkk=1:length(info(kk).atom_numbers),
            anum=info(kk).atom_numbers{kkk};
            [ma,na]=size(anum);
            for k4=1:ma,
                poi=poi+1;
                elements(poi)=ch_elements(anum(k4,1));
            end;
        end;
    end;
end;
elements=elements(1:poi);

function [message,elements]=elements_heavy_chain_model(indices,varargin)
% returns atomic numbers of elements for a chain model
% water can be excluded by setting varargin{1}='-nowater';
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';


max_atoms=100000;
elements=zeros(max_atoms,1);

poi=0;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
for kk=1:length(info), % loop over residues
    if ~strcmpi(info(kk).name,'HOH'),
        [message,celements]=get_residue([indices,kk],'elements_heavy');
        if message.error ~= 0,
            elements = [];
            return
        end;
        elements(poi+1:poi+length(celements)) = celements;
        poi=poi+length(celements);
    end
end
elements=elements(1:poi);

function [message,populations]=populations_chain_model(indices,varargin)
% returns populations of atom locations for a chain model
% water can be excluded by setting varargin{1}='-nowater';
%
% indices  index vector, length at least 3

global model

message.error=0;
message.text='';

if nargin>=2,
    argin=varargin{1};
else
    argin='';
end;

if isempty(argin),
    waterflag=true;
else
    if strcmp(argin,'-nowater')
        waterflag=false;
    else
        waterflag=true;
    end;
end;

max_atoms=100000;
populations=zeros(max_atoms,1);

poi=0;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info;
for kk=1:length(info), % loop over residues
    if waterflag || ~strcmpi(info(kk).name,'HOH'),
        for kkk=1:length(info(kk).atom_numbers),
            anum=info(kk).atom_numbers{kkk};
            [ma,na]=size(anum);
            if na==1,
                anum=[anum 1];
            end;
            for k4=1:ma,
                poi=poi+1;
                populations(poi)=anum(k4,2);
            end;
        end;
    end;
end;
populations=populations(1:poi);

