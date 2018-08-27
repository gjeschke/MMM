function [message,argout]=get_chain(indices,property,varargin)
% function [message,argout]=get_chain(indices,property,varargin)
%
% Gets properties of a chain in MMM
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

message.error=0;
message.text='';

if nargin>=3,
    argin=varargin{1};
else
    argin={};
end;

switch property
    case 'annotations'
	[message,argout]=get_annotations(indices);
    case 'children'
        [message,argout]=children_chain(indices);
    case 'descendants'
        [message,argout]=descendants_chain(indices);
    case 'info'
        [message,argout]=info_chain(indices);
    case 'mass'
        [message,argout]=mass_chain(indices,argin);
    case 'xyz'
        [message,argout]=xyz_chain(indices,argin);
    case 'xyz_heavy'
        [message,argout]=xyz_heavy_chain(indices,argin);
    case 'populations'
        [message,argout]=populations_chain(indices,argin);
    case 'elements'
        [message,argout]=elements_chain(indices,argin);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function [message,newindices]=children_chain(indices)
% returns children indices for a chain
%
% indices  index vector, length 2

global model

message.error=0;
message.text='';

chain_models=model.structures{indices(1)}(indices(2)).residues;
% keyboard;
if isempty(chain_models)
    newindices=indices;
else
    newindices=zeros(length(chain_models),3);
    for k=1:length(chain_models),
        newindices(k,1:2)=indices;
        newindices(k,3)=k;
    end;
end;

function [message,newindices]=descendants_chain(indices)
% returns descendants indices for a chain
%
% indices  index vector, length 2

[message,mindices]=children_chain(indices);
[m,n]=size(mindices);
if n>2, % there are chain model children, we have to check for children of their children
    newindices=zeros(10000,6);
    newindices(1:m,1:n)=mindices;
    poi=m;
    maxdim=3;
    for k=1:m,
        [message,nindices]=get_chain_model(mindices(k,:),'descendants');
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

function [message,info_text]=info_chain(indices)
% Plots an atom location by calling plot routine
%
% mode  string that determines the appearance of the plot

global model

message.error=0;
message.text='';

stag=id2tag(indices(1),model.structure_tags,model.structure_ids);
ctag=id2tag(indices(2),model.chain_tags{indices(1)},model.chain_ids{indices(1)});
adr_string=sprintf('[%s](%s)',stag,ctag);
info_text={adr_string sprintf('Structure: %s',stag) sprintf('Chain      : %s',ctag)};


function [message,mass]=mass_chain(indices,varargin)
% returns mass for a chain
%
% indices  index vector, length at least 2

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

% mass=sum(model.structures{indices(1)}(indices(2)).isotopes(:,2)); % old
% version

mindices=[indices(1:2) 1];
[message,mass]=get_chain_model(mindices,'mass',argin);

function [message,xyz]=xyz_chain(indices,varargin)
% returns xyz coordinates for a chain 
% if there are several models, the mean coordinates are returned
%
% indices  index vector, length at least 2

global model


if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

coor=model.structures{indices(1)}(indices(2)).xyz;
n=length(coor);
mindices=[indices(1:2) 1];
[message,xyz]=get_chain_model(mindices,'xyz',argin);
if n>1,
    for k=2:n,
        mindices(3)=k;
        [message,xyz0]=get_chain_model(mindices,'xyz',argin);
        xyz=xyz+xyz0;
    end;
end;
xyz=xyz/n;
% xyz=coor{1}; % old version
% if n>1,
%     for k=2:n,
%         xyz=xyz+coor{k};
%     end;
%     xyz=xyz/n;
% end;

function [message,xyz]=xyz_heavy_chain(indices,varargin)
% returns xyz coordinates for a chain 
% if there are several models, the mean coordinates are returned
%
% indices  index vector, length at least 2

global model

if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

coor=model.structures{indices(1)}(indices(2)).xyz;
n=length(coor);
mindices=[indices(1:2) 1];
[message,xyz]=get_chain_model(mindices,'xyz_heavy',argin);
if n>1,
    for k=2:n,
        mindices(3)=k;
        [message,xyz0]=get_chain_model(mindices,'xyz',argin);
        xyz=xyz+xyz0;
    end;
end;
xyz=xyz/n;

function [message,elements]=elements_chain(indices,varargin)
% returns atomic numbers of elements for a chain 
%
% indices  index vector, length at least 2


if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

mindices=[indices(1:2) 1];
[message,elements]=get_chain_model(mindices,'elements',argin);

% elements=model.structures{indices(1)}(indices(2)).isotopes(:,1); % old
% version

function [message,populations]=populations_chain(indices,varargin)
% returns populations coordinates for a chain 
% if there are several models, the mean populations are returned
%
% indices  index vector, length at least 2

global model


if nargin>=2,
    argin=varargin{1};
else
    argin={};
end;

coor=model.structures{indices(1)}(indices(2)).xyz;
n=length(coor);
mindices=[indices(1:2) 1];
[message,populations]=get_chain_model(mindices,'populations',argin);
if n>1,
    for k=2:n,
        mindices(3)=k;
        [message,pop0]=get_chain_model(mindices,'populations',argin);
        populations=populations+pop0;
    end;
end;
populations=populations/n;
