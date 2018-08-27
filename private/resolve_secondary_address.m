function secondary_indices=resolve_secondary_address(address,indices)
% function secondary_indices=resolve_secondary_address(address,indices)
%
% Secondary structure element indices for the input address
% 
% address               object address in MMM format
% indices               object indices as obtained by resolve_address to
%                       identify structure, chain, and chain model
%
% secondary_indices     indices that identify secondary object within the
%                       given structure

global model

secondary_indices=[];

if length(indices)<2, % secondary structure elements can only be addressed if structure and chain are known
    return;
else
    snum=indices(1);
    cnum=indices(2);
end;

secondary_request=0;

t_start=strfind(address,'<');
t_mid=strfind(address,'>');
t_stop=strfind(address,'.');

if isempty(t_start)  % a secondary structure element request does not exist
    return
else
    secondary_request=1;
    if isempty(t_mid) || t_mid<t_start,
        return
    end;
    if t_mid-t_start<4, % no place for type and name
        return
    end;
    if isempty(t_stop), % no dot separator between type and name
        return
    end;
    if t_stop(1)-t_start~=2
        return
    else
        sec_type=address(t_start+1);
        switch sec_type
            case {'L','l'}
                sec=0;
            case {'H','h'}
                sec=1;
            case {'E','e'}
                sec=2;
            otherwise
                return
        end;
    end;
    if t_mid-t_stop(1)<2,
        return
    end;
    sec_tag=address(t_stop(1)+1:t_mid-1);
    if length(t_stop)<2,
        t_stop=length(address)+1;
    else
        t_stop=t_stop(2);
    end;
    if t_stop-t_mid>2
        return % residues or atoms are specified, hence no whole sec. structure element selected
    end;
end;

sec_tags=':';
switch upper(sec_type)
    case 'H'
        if isfield(model.structures{snum}(cnum),'helix_defs'),
            secnum=length(model.structures{snum}(cnum).helix_defs);
            sec_ids=1:secnum;
            for sk=1:secnum,
                sec_tags=[sec_tags model.structures{snum}(cnum).helix_defs{sk}.name ':'];
            end;
            sec_ind=tag2id(sec_tag,sec_tags,sec_ids);
        end;
    case 'E'
        if isfield(model.structures{snum}(cnum),'sheet_defs'),
            secnum=length(model.structures{snum}(cnum).sheet_defs);
            sec_ids=1:secnum;
            for sk=1:secnum,
                sec_tags=[sec_tags model.structures{snum}(cnum).sheet_defs{sk}.name ':'];
            end;
            sec_ind=tag2id(sec_tag,sec_tags,sec_ids);
        end;
    case 'L'
        if isfield(model.structures{snum}(cnum),'loop_defs'),
            secnum=length(model.structures{snum}(cnum).loop_defs);
            sec_ids=1:secnum;
            for sk=1:secnum,
                sec_tags=[sec_tags model.structures{snum}(cnum).loop_defs{sk}.name ':'];
            end;
            sec_ind=tag2id(sec_tag,sec_tags,sec_ids);
        end;
end;

if ~isempty(sec_ind),
    secondary_indices=[sec sec_ind];
end;