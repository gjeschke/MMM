function message=domain(tag,create,address)
% function message=domain(tag,create,address)
%
% Defines a domain in an MMM model
% by either a robust address or a secondary structure tag
%
% tag       tag (name) of the domain
% create    flag, if 1, the domain is newly defined (an old domain with the
%           same name is overwritten), if create=0 the function attempts to
%           add to the existing domain with the same name or creates a new
%           domain if the name does not yet exist
%           use 2 to check whether the domain already exists
%           (message.error=8, message.text='Domain exists';), otherwise
%           message.error<8;
% address   address of the objects to be added to the domain
%
% message   error message structure with fields
%           .error
%           .text
%
% G. Jeschke, 2009

global model

message.error=0;
message.text='';

[indices,message]=resolve_address(address);
if message.error>0,
    return;
end;

if isempty(indices),
    message.error=3;
    message.text='No objects addressed.';
    return
end;

[m,n]=size(indices);
mix=0;
for k=1:m,
    cindices=indices(k,:);
    if k==1,
        level=length(cindices(cindices>0));
    else
        if level~=length(cindices(cindices>0)),
            mix=1; break;
        end;
    end;
end;
if mix,
    message.error=2;
    message.text='Objects on different levels cannot be mixed in a domain.';
    return;
end;
if level~=4,
    message.error=7;
    message.text='Domains should be defined at residue level.';
    return
end;
snum=indices(1,1);
if sum(indices(:,1)-snum),
    message.error=4;
    message.text='Domain must not contain objects from different structures.';
end;

if ~isfield(model,'domains') % this is the first domain to be defined
    model.domains{snum}(1).name=tag;
    model.domains{snum}(1).objects=indices;
    dnum=1;
    subindices=1:m;
elseif length(model.domains)<snum % its also new
    model.domains{snum}(1).name=tag;
    model.domains{snum}(1).objects=indices;
    dnum=1;
    subindices=1:m;
else
    doms=length(model.domains{snum});
    new_define=true;
    dnum=doms+1;
    if doms>0,
        for k=1:doms,
            if strcmpi(tag,model.domains{snum}(k).name) % there is already a domain with that name
                dnum=k;
                if ~create,
                    new_define=false;
                end;
                if create==2,
                    message.error=8;
                    message.text='Domain exists.';
                    return;
                end;
            end;
        end;
    end;
    if create==2, return; end;
    if new_define,
        model.domains{snum}(dnum).name=tag;
        model.domains{snum}(dnum).objects=indices;
        subindices=1:m;
    else
        old_indices=model.domains{snum}(dnum).objects;
        [m0,n0]=size(old_indices);
        if n~=n0,
            message.error=5;
            message.text='Hierarchy level mismatch.';
            return
        else
            indices=[old_indices; indices];
            model.domains{snum}(dnum).objects=indices;
            subindices=m0+1:m0+m;
        end;
    end;
end;

% Process name to see if secondary structure element is to be defined as
% subdomain

secondary_request=false;
t_start=strfind(address,'<');
t_mid=strfind(address,'>');
if ~isempty(t_start) && ~isempty(t_mid) % could be secondary structure tag
    if t_mid<t_start+2,
        message.error=6;
        message.text='Invalid secondary structure tag.';
        return
    else
        subtag=address(t_start+1:t_mid-1);
        secondary_request=true;
    end;
end;

if secondary_request,
    if ~isfield(model.domains{snum}(dnum),'subdomains') % this is the first subdomain to be defined
        model.domains{snum}(dnum).subdomains{1}.indices=subindices;
        model.domains{snum}(dnum).subdomains{1}.name=subtag;
        sdnum=1;
    else
        subdoms=length(model.domains{snum}(dnum).subdomains);
        new_define=true;
        sdnum=subdoms+1;
        if subdoms>0,
            for k=1:subdoms,
                if strcmpi(subtag,model.domains{snum}(dnum).subdomains{k}.name) % there is already a subdomain with that name
                    sdnum=k;
                    if ~create,
                        new_define=false;
                    end;
                end;
            end;
        end;
        if new_define,
            model.domains{snum}(dnum).subdomains{sdnum}.name=subtag;
            model.domains{snum}(dnum).subdomains{sdnum}.indices=subindices;
        else
            old_subindices=model.domains{snum}(dnum).subdomains{sdnum}.indices;
            subindices=[old_subindices subindices];
            model.domains{snum}(dnum).subdomains{sdnum}.indices=subindices;
        end;
    end;
end;

