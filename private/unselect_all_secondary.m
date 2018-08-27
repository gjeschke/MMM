function unselect_all_secondary
% function unselect_all_secondary
%
% unselects all secondary structure elements in the model
%
% G. Jeschke, 2009

global model

structures=length(model.structures);

for snum=1:structures,
    chains=length(model.structures{snum});
    for cnum=1:chains,
        loops=length(model.structures{snum}(cnum).loop_defs);
        if loops>0,
            for k=1:loops,
                model.structures{snum}(cnum).loop_defs{k}.selected=0;
            end;
        end;
        helices=length(model.structures{snum}(cnum).helix_defs);
        if helices>0,
            for k=1:helices,
                model.structures{snum}(cnum).helix_defs{k}.selected=0;
            end;
        end;
        sheets=length(model.structures{snum}(cnum).sheet_defs);
        if sheets>0,
            for k=1:sheets,
                model.structures{snum}(cnum).sheet_defs{k}.selected=0;
            end;
        end;
    end;
end;