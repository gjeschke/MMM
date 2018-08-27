function message=secondary(address,type,tag,range,TM)
% function message=secondary(address,type,tag,range,TM)
%
% Defines a secondary structure element (SSE) in an MMM model
% this can be a helix, sheet or loop
% the function does not test whether these residues indeed exist
%
% address   address of the chain(s) for which the element is to be defined
% type      can be 'helix','sheet',or 'loop'
% tag       tag for the SSE
% range     continuous range of residues that define the SSE, array[1,2],
%           only first and last residue are given
% TM        optional flag indicating a transmembrane helix or sheet,
%           defaults to 0, is ignored for loops
%
% message   error message structure with fields
%           .error
%           .text
%
% G. Jeschke, 2009

global model
global graph_settings

message.error=0;
message.text='';

if nargin<5,
    TM=0;
end;

if length(range)~=2,
    message.error=3;
    message.text='Wrong residue range dimension.';
    return;
end;

if address(end)==']',
    address=sprintf('%s(%c)',address,model.current_chain);
end;

[indices,message]=resolve_address(address);
if message.error>0,
    return;
end;

[m,n]=size(indices);
if n~=2,
    message.error=4;
    message.text='First argument does not address a chain.';
    return
end;

for k=1:m,
    snum=indices(k,1);
    cnum=indices(k,2);
    switch type
        case {'helix','HELIX'}
            if isfield(model.structures{snum}(cnum),'helix_defs')
                poi=length(model.structures{snum}(cnum).helix_defs)+1;
                for kk=1:poi-1, % avoid double definition of same tag
                    if strcmp(model.structures{snum}(cnum).helix_defs{kk}.name,tag)
                        poi=kk;
                    end;
                end;
                model.structures{snum}(cnum).helix_defs{poi}.name=tag;
                model.structures{snum}(cnum).helix_defs{poi}.range=range;
                model.structures{snum}(cnum).helix_defs{poi}.TM=TM;
            else
                model.structures{snum}(cnum).helix_defs{1}.name=tag;
                model.structures{snum}(cnum).helix_defs{1}.range=range;
                model.structures{snum}(cnum).helix_defs{1}.TM=TM;
                poi=1;
            end;
            if TM,
                colind=resolve_address(sprintf('%s{:}%i-%i',mk_address([snum cnum]),range(1),range(2)));
                set_object(colind,'ribboncolor',{graph_settings.TM_helix_color});
            end;
            coor_sets=length(model.structures{snum}(cnum).residues);
            for kk=1:coor_sets,
                for curr_res=range(1):range(2),
                    cr_tag=sprintf('%i',curr_res);
                    id=tag2id(cr_tag,model.structures{snum}(cnum).residues{kk}.residue_tags);
                    if ~isempty(id),
                        model.structures{snum}(cnum).residues{kk}.info(id).secondary=1;
                    end;
                end;
            end;
        case {'sheet','SHEET'}
            if isfield(model.structures{snum}(cnum),'sheet_defs')
                poi=length(model.structures{snum}(cnum).sheet_defs)+1;
                for kk=1:poi-1, % avoid double definition of same tag
                    if strcmp(model.structures{snum}(cnum).sheet_defs{kk}.name,tag)
                        poi=kk;
                    end;
                end;
                model.structures{snum}(cnum).sheet_defs{poi}.name=tag;
                model.structures{snum}(cnum).sheet_defs{poi}.range=range;
                model.structures{snum}(cnum).sheet_defs{poi}.TM=TM;
            else
                model.structures{snum}(cnum).sheet_defs{1}.name=tag;
                model.structures{snum}(cnum).sheet_defs{1}.range=range;
                model.structures{snum}(cnum).sheet_defs{1}.TM=TM;
                poi=1;
            end;
            if TM,
                colind=resolve_address(sprintf('%s{:}%i-%i',mk_address([snum cnum]),range(1),range(2)));
                set_object(colind,'ribboncolor',{graph_settings.TM_sheet_color});
            end;
            coor_sets=length(model.structures{snum}(cnum).residues);
            for kk=1:coor_sets,
                for curr_res=range(1):range(2),
                    cr_tag=sprintf('%i',curr_res);
                    id=tag2id(cr_tag,model.structures{snum}(cnum).residues{kk}.residue_tags);
                    if ~isempty(id),
                        model.structures{snum}(cnum).residues{kk}.info(id).secondary=2;
                    end;
                end;
            end;
        case {'loop','LOOP'}
            if isfield(model.structures{snum}(cnum),'loop_defs')
                poi=length(model.structures{snum}(cnum).loop_defs)+1;
                for kk=1:poi-1, % avoid double definition of same tag
                    if strcmp(model.structures{snum}(cnum).loop_defs{kk}.name,tag)
                        poi=kk;
                    end;
                end;
                model.structures{snum}(cnum).loop_defs{poi}.name=tag;
                model.structures{snum}(cnum).loop_defs{poi}.range=range;
            else
                model.structures{snum}(cnum).loop_defs{1}.name=tag;
                model.structures{snum}(cnum).loop_defs{1}.range=range;
            end;
            coor_sets=length(model.structures{snum}(cnum).residues);
            for kk=1:coor_sets,
                for curr_res=range(1):range(2),
                    cr_tag=sprintf('%i',curr_res);
                    id=tag2id(cr_tag,model.structures{snum}(cnum).residues{kk}.residue_tags);
                    if ~isempty(id),
                        model.structures{snum}(cnum).residues{kk}.info(id).secondary=0;
                    end;
                end;
            end;
        otherwise
            message.error=5;
            message.text=sprintf('Unknown secondary structure type %s.',type);
            return
    end;
end;