function select(indices)
% selects objects addressed by an index array

global model

[m,n]=size(indices);
for k=1:m,
    idepth=length(find(indices(k,:)>0)); % determine type of current object
    cindices=indices(k,1:idepth);
    sec_indices=resolve_secondary_address(address,indices);
    if ~isempty(sec_indices),
        snum=indices(1);
        cnum=indices(2);
        switch sec_indices(1)
            case 0
                model.structures{snum}(cnum).loop_defs{sec_indices(2)}.selected=1;
            case 1
                model.structures{snum}(cnum).helix_defs{sec_indices(2)}.selected=1;
            case 2
                model.structures{snum}(cnum).sheet_defs{sec_indices(2)}.selected=1;
        end;
    end;
    if model.selections==0,
        model.selections=1;
        model.selected{1}=cindices;
    else
        store=1;
        for kk=1:model.selections,
            sindices=model.selected{kk};
            sdepth=length(find(sindices>0)); % determine type of current object
            sindices=sindices(1:sdepth);
            if sdepth==idepth
                if cindices==sindices
                    store=0;
                    add_msg_board('Double selection of same object ignored');
                end;
            end;
        end;
        if store
            model.selections=model.selections+1;
            model.selected{model.selections}=cindices;
        end;
    end;
end;

highlight_selection;
