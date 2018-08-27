function flag=is_selected(indices)
% function flag=is_selected(indices)
%
% Determines whether an object with given indices is selected,
% selection of an ancestor of this object is considered
%
% indices  indices of the object
% flag     1 object is selected, 0 object is not selected
%
% G. Jeschke, 2009

global model

flag=0;

if isempty(model.selected),
    return;
end;

for k=1:length(model.selected),
    cindices=model.selected{model.selections};
    milen=min([length(cindices),length(indices)]);
    if cindices(1:milen)==indices(1:milen),
        flag=true;
        break;
    end;
end;
