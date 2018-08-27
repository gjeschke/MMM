function indices = sorted_selection
% Provides a sorted list of the indices of all selected objects down to
% atom location level

global model

indices = zeros(100000,6);
poi = 0;

% extend indices to location array

for ko = 1:length(model.selected), % loop over all objects
    oindices = model.selected{ko};
    idepth = length(find(oindices>0)); % determine type of current object
    cindices = oindices(1:idepth);
    allindices = cindices;
    if length(cindices)<4,
        if length(cindices)<2,
            adr  =mk_address(cindices(1));
            adr = sprintf('%s(:){:}:',adr);
            allindices = resolve_address(adr);
        elseif length(cindices)<3,
            adr = mk_address(cindices(1:2));
            adr = sprintf('%s{:}:',adr);
            allindices = resolve_address(adr);
        else
            adr = mk_address(cindices(1:3));
            adr = sprintf('%s:',adr);
            allindices = resolve_address(adr);
        end;
    end;
    idepth = length(find(allindices(1,:)>0));
    [ma,~] = size(allindices);
    for kr = 1:ma,
        newindices = allindices(kr,:);
        if idepth == 4,
            [~,newindices] = get_residue(allindices(kr,:),'descendants');
        elseif idepth == 5
            [~,newindices] = get_atom(allindices(kr,:),'descendants');
        end;
        [~,nn] = size(newindices); 
        if nn == 6
            newindices = newindices(newindices(:,6)~=0,:);
            [mn,~] = size(newindices);
            indices(poi+1:poi+mn,:) = newindices;
            poi = poi + mn;
        end
    end;
end

indices = indices(1:poi,:);

[~,order] = sort(indices(:,1));
indices = indices(order,:);
mins = min(indices(:,1));
maxs = max(indices(:,1));
for k1 = mins:maxs
    sel1 = indices(indices(:,1) == k1,:);
    [~,ord2] = sort(sel1(:,2));
    sel1 = sel1(ord2,:);
    mins2 = min(sel1(:,2));
    maxs2 = max(sel1(:,2));
    for k2 = mins2:maxs2
        sel2 = sel1(sel1(:,2) == k2,:);
        [~,ord3] = sort(sel2(:,3));
        sel2 = sel2(ord3,:);
        mins3 = min(sel2(:,3));
        maxs3 = max(sel2(:,3));
        for k3 = mins3:maxs3
            sel3 = sel2(sel2(:,3) == k3,:);
            [~,ord4] = sort(sel3(:,4));
            sel3 = sel3(ord4,:);
            mins4 = min(sel3(:,4));
            maxs4 = max(sel3(:,4));
            for k4 = mins4:maxs4
                sel4 = sel3(sel3(:,4) == k4,:);
                [~,ord5] = sort(sel4(:,5));
                sel4 = sel4(ord5,:);
                mins5 = min(sel4(:,5));
                maxs5 = max(sel4(:,5));
                for k5 = mins5:maxs5
                    sel5 = sel4(sel4(:,4) == k5,:);
                    [~,ord6] = sort(sel5(:,6));
                    sel4(sel4(:,4) == k5,:) = sel5(ord6,:);
                end;
                sel3(sel3(:,4) == k4,:) = sel4;
            end;
            sel2(sel2(:,3) == k3,:) = sel3;
        end;
        sel1(sel1(:,2) == k2,:) = sel2;
    end;
    indices(indices(:,1) == k1,:) = sel1;
end;

