function [message,argout]=get_object(address,property,varargin)
% function [message,argout]=get_object(address,property,varargin)
%
% Gets properties of objects (any type) in MMM
% uses resolve_address for address resolution and get functions of 
% appropriate object types
%
% address   object address string, may address several objects at once or
%           index vector or array of indices
% property  property to be set, given as a string, available properties are 
%           defined in the set functions of the individual object types
% varargin  further input arguments given as a cell array varargin{:}
%
% message   error message structure with fields .error and .text,
%           .error=0 indicates no error
% argout    further output arguments given as structure argout{n}
%
% G. Jeschke, 2009

argout={};

message.error=0;
message.text='';

if isa(address,'char'),
    [indices,message]=resolve_address(address);
else
    indices=address;
end;

[m,n]=size(indices);

if m==0 || n==0 || message.error~=0,
    if message.error==0 || message.error==2,
       message.error=2;
       message.text='Object does not exist';
    end;
    return
end;

argout={};
if nargin>=3,
    argin=varargin{1};
else
    argin='';
end;
for k=1:m, % loop over all objects
    idepth=length(find(indices(k,:)>0)); % determine type of current object
    cindices=indices(k,1:idepth);
    switch idepth
        case 1
            [message,ao]=get_structure(cindices,property,argin);
        case 2
            [message,ao]=get_chain(cindices,property,argin);
        case 3
            [message,ao]=get_chain_model(cindices,property,argin);
        case 4
            [message,ao]=get_residue(cindices,property,argin);
        case 5
            [message,ao]=get_atom(cindices,property,argin);
        case 6
            [message,ao]=get_location(cindices,property,argin);
        otherwise
            message.error=4;
            message.text='Unknown object type';
    end;
    argout{k}=ao;
end;

if m==1,
    argout=ao;
end;
