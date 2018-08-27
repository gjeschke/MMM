function [message,argout]=get_structure(indices,property,varargin)
% function [message,argout]=get_structure(indices,property,varargin)
%
% Gets properties of a structure in MMM
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
    argin={};
end;

message.error=0;
message.text='';

switch property
    case 'annotations'
	[message,argout]=get_annotations(indices);
    case 'children'
        [message,argout]=children_structure(indices);
    case 'descendants'
        [message,argout]=descendants_structure(indices);
    case 'info'
        [message,argout]=info_structure(indices);
    case 'mass'
        [message,argout]=mass_structure(indices,argin);
    case 'xyz'
        [message,argout]=xyz_structure(indices,argin);
    case 'xyz_heavy'
        [message,argout]=xyz_heavy_structure(indices,argin);
    case 'populations'
        [message,argout]=xyz_structure(indices,argin);
    case 'elements'
        [message,argout]=elements_structure(indices,argin);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function [message,newindices]=children_structure(indices)
% returns children indices for a structure
%
% indices  index vector, length 1

global model

message.error=0;
message.text='';

chains=model.structures{indices(1)};
if isempty(chains)
    newindices=indices;
else
    newindices=zeros(length(chains),2);
    for k=1:length(chains),
        newindices(k,1)=indices;
        newindices(k,2)=k;
    end;
end;

function [message,newindices]=descendants_structure(indices)
% returns descendants indices for a structure
%
% indices  index vector, length 1

message.error=0;
message.text='';

[message,mindices]=children_structure(indices);
[m,n]=size(mindices);
if n>1, % there are chain children, we have to check for children of their children
    newindices=zeros(50000,6);
    newindices(1:m,1:n)=mindices; % store the existing indices
    poi=m;
    maxdim=2;
    for k=1:m,
        if mindices(k,2)>0,
            [message,nindices]=get_chain(mindices(k,:),'descendants');
            [mm,nn]=size(nindices);
            if nn>maxdim,
                maxdim=nn;
            end;
            newindices(poi+1:poi+mm,1:nn)=nindices;
            poi=poi+mm;
        end;
    end;
    newindices=newindices(1:poi,1:maxdim);
else
    newindices=indices;
end;

function [message,info_text]=info_structure(indices)
% Plots an atom location by calling plot routine
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

stag=id2tag(indices(1),model.structure_tags,model.structure_ids);
adrstring=sprintf('[%s]',stag);
info_text={adrstring sprintf('Structure: %s',stag)};

function [message,mass]=mass_structure(indices,varargin)
% returns mass for a chain
%
% indices  index vector, length at least 1

global model

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

message.error=0;
message.text='';

nc=length(model.structures{indices(1)}); % number of chains
mass=0;
for k=1:nc,
    chindices=[indices(1),k];
    [message,mass_ch]=get_chain(chindices,'mass',argin);
    mass=mass+mass_ch;
end;

function [message,xyz]=xyz_structure(indices,varargin)
% returns xyz coordinates for a structure 
% if there are several models, the mean coordinates are returned
%
% indices  index vector, length at least 1

global model

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

message.error=0;
message.text='';

max_atoms=100000;
xyz=zeros(max_atoms,3);

poi=0;

nc=length(model.structures{indices(1)}); % number of chains
for k=1:nc,
    chindices=[indices(1),k];
    [message,xyz_ch]=get_chain(chindices,'xyz',argin);
    [m,n]=size(xyz_ch);
    xyz(poi+1:poi+m,:)=xyz_ch;
    poi=poi+m;
end;
xyz=xyz(1:poi,:);

function [message,xyz]=xyz_heavy_structure(indices,varargin)
% returns xyz coordinates for a structure 
% if there are several models, the mean coordinates are returned
%
% indices  index vector, length at least 1

global model

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

message.error=0;
message.text='';

max_atoms=100000;
xyz=zeros(max_atoms,3);

poi=0;

nc=length(model.structures{indices(1)}); % number of chains
for k=1:nc,
    chindices=[indices(1),k];
    [message,xyz_ch]=get_chain(chindices,'xyz_heavy',argin);
    [m,n]=size(xyz_ch);
    xyz(poi+1:poi+m,:)=xyz_ch;
    poi=poi+m;
end;
xyz=xyz(1:poi,:);


function [message,elements]=elements_structure(indices,varargin)
% returns atomic numbers of elements for a structure 
%
% indices  index vector, length at least 1

global model

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

message.error=0;
message.text='';

max_atoms=100000;
elements=zeros(max_atoms,1);

poi=0;

nc=length(model.structures{indices(1)}); % number of chains
for k=1:nc,
    chindices=[indices(1),k];
    [message,elm]=get_chain(chindices,'elements',argin);
    n=length(elm);
    elements(poi+1:poi+n,:)=elm;
    poi=poi+n;
end;
elements=elements(1:poi,:);

function [message,populations]=populations_structure(indices,varargin)
% returns populations of all atom locations for a structure 
%
% indices  index vector, length at least 1

global model

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

message.error=0;
message.text='';

max_atoms=100000;
populations=zeros(max_atoms,1);

poi=0;

nc=length(model.structures{indices(1)}); % number of chains
for k=1:nc,
    chindices=[indices(1),k];
    [message,elm]=get_chain(chindices,'populations',argin);
    n=length(elm);
    populations(poi+1:poi+n,:)=elm;
    poi=poi+n;
end;
populations=populations(1:poi,:);

