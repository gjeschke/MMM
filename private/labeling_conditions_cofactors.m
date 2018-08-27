function labeling_conditions_cofactors(indices)

global hMain
global label_defs
global rotamer_libraries

% make list of cofactor label definitions
to_check = zeros(1,length(label_defs.residues));
poi = 0;
for k = 1:length(label_defs.residues),
    attach = label_defs.residues(k).attachment;
    if attach(1)==':' && attach(end)==':',
        poi = poi+1;
        to_check(poi) = k;
    end;
end;
to_check = to_check(1:poi);

[m,n] = size(indices);

hMain.temperature = [];
hMain.library = '';

if n ~= 4,
    add_msg_board('ERROR: Cofactor labeling is allowed only for selection of residues');
    return
end;

for k = 1:m,
    libs{k} = [];
    [msg,name] = get_object(indices(k,:),'name');
    if msg.error,
        add_msg_board(msg.text);
        return;
    end;
    found = false;
    for kd = 1:length(to_check),
        attach = label_defs.residues(to_check(kd)).attachment;
        id = tag2id(name,attach);
        if ~isempty(id),
            resname = label_defs.residues(to_check(kd)).resname;
            for kr = 1:length(rotamer_libraries),
                if strcmpi(resname,rotamer_libraries(kr).tc),
                    found = true;
                    libs{k} = [libs{k} kr];
                end;
            end;
        end;
    end;
    if ~found,
        adr = mk_address(indices(k,:));
        add_msg_board(sprintf('ERROR: No label definition for residue type %s of selected residue %s.\n',name,adr));
        return
    end;
end;
common_libs = libs{1};
for k = 2:m,
    new_cl = zeros(1,length(common_libs));
    poi = 0;
    for kl = 1:length(common_libs),
        for kc = 1:length(libs{k}),
            if common_libs(kl) == libs{k}(kc),
                poi = poi+1;
                new_cl(poi)= common_libs(kl);
            end;
        end;
        if poi == 0,
            add_msg_board('ERROR: No common rotamer library for all selected residues.');
            return
        else
            common_libs = new_cl(1:poi);
        end;
    end;
end;

poi = 0;
for k = 1:length(common_libs),
    for kl = 1:length(rotamer_libraries(common_libs(k)).T),
        poi = poi + 1;
        T = rotamer_libraries(common_libs(k)).T;
        library = id2tag(kl,rotamer_libraries(common_libs(k)).files);
        lib_selection{poi} = sprintf('%s at %iK',library,T);
    end;
end;

if poi == 1,
    hMain.temperature = T;
    hMain.library = library;
end;