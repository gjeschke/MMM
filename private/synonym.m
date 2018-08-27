function message=synonym(address,syn_tag)
% function message=synonym(address,syn_tag)
%
% Definition of a synonym for a structure or chain with an MMM model,
% existing synonyms are overwritten without warning
%
% address   address of the form [stag], [stag](ctag), or (ctag), where stag
%           is an alredy defined structure tag and ctag an already defined
%           chain tag, objects on coordinate set, residue, atom, or
%           location level cannot be given synonyms
% syn_tag   synonymous tag for the structure or chain
% message   error message, if any
%
% G. Jeschke, 2009

global model

[indices,message]=resolve_address(address);
if message.error,
    add_msg_board('ERROR: Wrong address format.');
    add_msg_board(message.text);
    add_msg_board('Use [old_structure_name] or (old_chain_name).');
    return
end;

if isempty(indices) % nothing proper was addressed
    message=no_object;
    add_msg_board('ERROR: Object does not exist.');
    add_msg_board('Use [old_structure_name] or (old_chain_name).');
    return;
end;

[m,n]=size(indices); % something was addressed

if n>2,
    message.error=3;
    message.text='No synonyms can be defined below chain level';
    add_msg_board(message.text);
    return
end;

if n==1, % structure synonym
    if m>1, 
        message.error=4;
        message.text='A synonym cannot apply to more than one structure';
        add_msg_board(message.text);
        return
    end;
    if isfield(model,'structure_tags') && isfield(model,'structure_ids') % check whether synonyms are defined
            tags=model.structure_tags;
            ids=model.structure_ids;
    else % generate tags string and code list
        tags=':';
        for k=1:length(model.structures),
            tags=[tags sprintf('%i:',k)];
        end;
        ids=1:length(model.structures);
    end;
    syn_ind=tag2id(upper(syn_tag),upper(tags),ids);
    if isempty(syn_ind) % new synonym
        tags=[tags sprintf('%s:',syn_tag)];
        model.structure_tags=tags;
        ids=[ids indices(1,1)];
    else
        ids(syn_ind)=indices(1,1);
    end;
    model.structure_ids=ids;
end;

if n==2, % chain synonym
    struct_test=zeros(1,length(model.structures));
    for k=1:m,
        snum=indices(k,1);
        if struct_test(snum) % for this structure, a synonym was lready defined and is requested again
            message.error=5;
            message.text='A synonym cannot apply to more than one chain in the same structure';
            add_msg_board(message.text);
            return
        else
            struct_test(snum)=1;
            tags=model.chain_tags{snum};
            ids=model.chain_ids{snum};
            syn_ind=tag2id(upper(syn_tag),upper(tags),ids);
        end;
        if isempty(syn_ind) % new synonym
            tags=[tags sprintf('%s:',syn_tag)];
            model.chain_tags{snum}=tags;
            ids=[ids indices(k,2)];
        else
            ids(syn_ind)=indices(k,2);
        end;
        model.chain_ids{snum}=ids;
    end;
end;

function message=syntax_error
% returns the error message structure for a syntax error

message.error=1;
message.text='Invalid address.';

function message=no_object
% returns the error message structure if no object is addressed

message.error=2;
message.text='Object does not exist.';
add_msg_board(message.text);
