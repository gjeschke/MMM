function [message,argout]=set_object(address,property,argin)
% function [message,argout]=set_object(address,property,argin)
%
% Sets properties of objects (any type) in MMM
% uses resolve_address for address resolution and set functions of 
% appropriate object types
%
% address   object address string, may address several objects at once
% property  property to be set, given as a string, available properties are 
%           defined in the set functions of the individual object types
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

if strcmpi(property,'show') && strcmpi(argin{1},'cartoon'),
    if (n ~= 4 || m < 2),
        message.error = 3;
        message.text = 'Cartoon display is defined only for ranges of residues.';
        add_msg_board('ERROR: Cartoon display is defined only for ranges of residues.');
        return
    else
        plot_cartoon(indices);
    end;
end;
for k=1:m, % loop over all objects
    idepth=length(find(indices(k,:)>0)); % determine type of current object
    cindices=indices(k,1:idepth);
    switch idepth
        case 1
            set_structure(cindices,property,argin);
        case 2
            set_chain(cindices,property,argin);
        case 3
            set_chain_model(cindices,property,argin);
        case 4
            set_residue(cindices,property,argin);
        case 5
            set_atom(cindices,property,argin);
        case 6
            set_location(cindices,property,argin);
        otherwise
            message.error=4;
            message.text='Unknown object type';
    end;
end;

