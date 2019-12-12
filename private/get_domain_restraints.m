function [restraints,restrain,aux] = get_domain_restraints(fname)

global model
global residue_defs

snum = [];
if isfield(model,'current_structure')
    snum=model.current_structure;
end

restraints=rd_restraints(fname);

aux.ensemble=restraints.ensemble;

if isfield(restraints,'PDB')
    if ~strcmpi(model.info{snum}.idCode,restraints.PDB)
        id=tag2id(restraints.PDB,model.structure_tags);
        if ~isempty(id)
            snum=id;
            model.current_structure=snum;
        else
            button = questdlg(sprintf('Restraint file specifies structure %s, while current structure is %s. Do you want to load specified template?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
            if strcmp(button,'Yes')
                fname=get_pdb_file(restraints.PDB);
                [message,snum]=add_pdb(fname);
                if message.error
                    add_msg_board(sprintf('ERROR: Specified structure PDB file %s could not be retrieved from server',restraints.PDB));
                    add_msg_board(message.text);
                    set(hfig,'Pointer','arrow');
                    return
                else
                    model.current_structure=snum;
                end
            end
        end
    end
end

restraints.template = snum;
if ~isempty(snum)
    template_seq = model.structures{snum}(1).sequence;
else
    template_seq = '';
end

[restrain,monitor,restraints,number,bnumber,number_monitor,cancelled] = process_domain_restraints(restraints);

if cancelled
    add_msg_board('Processing of restraints cancelled.');
    set(hfig,'Pointer','arrow');
    return
end
aux.template_seq = template_seq;
aux.monitor = monitor;
aux.n_restraints = number;
aux.b_restraints = bnumber;
aux.n_monitor = number_monitor;
add_msg_board(sprintf('%i distribution restraints have been read.',number));

if isfield(restraints,'ensemble')
    aux.ensemble = restraints.ensemble;
else
    aux.ensemble = 0;
end
if isfield(restraints,'uncertainty')
    aux.p_model = restraints.uncertainty/restraints.rescale;
end
Na_indices = [];
if isfield(restraints,'Nanchor')
    [indices,message]=resolve_address(restraints.Nanchor);
    Na_indices = indices;
    if message.error
        add_msg_board(sprintf('ERROR: N-terminal anchor residue %s does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    end
    Nname = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
    id = tag2id(Nname,upper(residue_defs.restags));
    N_slc = residue_defs.single_letter_code(id);
    indices_p = indices;
    indices_p(4) = indices_p(4) - 1; % index of previous residue before N terminal anchor
    Npname = model.structures{indices_p(1)}(indices_p(2)).residues{indices_p(3)}.info(indices_p(4)).name;
    id = tag2id(Npname,upper(residue_defs.restags));
    Np_slc = residue_defs.single_letter_code(id);
    restraints.Nanchor_p = mk_address(indices_p);
    anchorNp = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor N atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorNp(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor CA atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorNp(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor C atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorNp(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Nanchor_p),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor O atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorNp(4,:) = coor;
    end
    anchorN = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Nanchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s N atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorN(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Nanchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s CA atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorN(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Nanchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s C atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorN(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Nanchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s O atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorN(4,:) = coor;
    end
    restraints.anchorN = anchorN;
    restraints.anchorNp = anchorNp;
    restraints.Nseq = [Np_slc N_slc];
else
    restraints.anchorN = [];
    restraints.anchorNp = [];
    restraints.Nseq = '';
end
Ca_indices = [];
if isfield(restraints,'Canchor'),
    [indices,message]=resolve_address(restraints.Canchor);
    Ca_indices = indices;
    if message.error,
        add_msg_board(sprintf('ERROR: C-terminal anchor residue %s does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    end
    Cname = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
    id = tag2id(Cname,upper(residue_defs.restags));
    C_slc = residue_defs.single_letter_code(id);
    indices_n = indices;
    indices_n(4) = indices_n(4) + 1; % index of next residue after C-terminal anchor
    Cnname = model.structures{indices_n(1)}(indices_n(2)).residues{indices_n(3)}.info(indices_n(4)).name;
    id = tag2id(Cnname,upper(residue_defs.restags));
    Cn_slc = residue_defs.single_letter_code(id);
    restraints.Canchor_n = mk_address(indices_n);
    anchorCn = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Canchor_n),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor N atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorCn(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Canchor_n),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor CA atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorCn(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Canchor_n),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor C atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorCn(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Canchor_n),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor O atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorCn(4,:) = coor;
    end
    anchorC = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Canchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s N atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorC(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Canchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s CA atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorC(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Canchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s C atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorC(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Canchor),'coor');
    if msg.error,
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s O atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    else
        anchorC(4,:) = coor;
    end
    restraints.anchorC = anchorC;
    restraints.anchorCn = anchorCn;
    restraints.Cseq = [C_slc Cn_slc];
else
    restraints.anchorC = [];
    restraints.anchorCn = [];
    restraints.Cseq = '';
end
restraints.Na_indices = Na_indices;
restraints.Ca_indices = Ca_indices;


function [restrain,monitor,restraints,number,bnumber,number_monitor,cancelled]=process_domain_restraints(restraints)

global rotamer_libraries
global ligand_libraries

cancelled=false;

restrain = [];
number = 0;
bnumber = 0;
number_monitor = 0;

if ~isfield(restraints,'sequence')
    add_msg_board('ERROR: Sequence specification is missing.'); 
    restrain=[];
    cancelled = true;
    return;
end

for k = 1:length(restraints.sequence)
    restrain(k).secondary = 0;
    restrain(k).label = [];
    restrain(k).cis = 0;
    restrain(k).r_beacon = [];
    restrain(k).r_intern = [];
    restrain(k).oligomer = [];
    restrain(k).depth = [];
end

monitor = restrain;

if ~isfield(restraints,'domain_start') || ~isfield(restraints,'domain_end'),
    add_msg_board('ERROR: Domain is not specified in restraint file.'); 
    restrain=[];
    cancelled = true;
    return;
else
    poi1 = strfind(restraints.domain_start,'}');
    if ~isempty(poi1);
        poi2 = strfind(restraints.domain_end,'}');
        if ~strcmpi(restraints.domain_start(1:poi1),restraints.domain_end(1:poi2)),
            add_msg_board('ERROR: Addresses of first and last residue of the domain are inconsistent.'); 
            add_msg_board(sprintf('%s does not match %s',restraints.domain_start(1:poi1),restraints.domain_end(1:poi2))); 
            restrain=[];
            cancelled = true;
            return;
        else
            chain_model = restraints.domain_start(1:poi1);
        end
    else
        poi1 = strfind(restraints.domain_start,')');
        if ~isempty(poi1);
            poi2 = strfind(restraints.domain_end,')');
            if ~strcmpi(restraints.domain_start(1:poi1),restraints.domain_end(1:poi2)),
                add_msg_board('ERROR: Addresses of first and last residue of the domain are inconsistent.'); 
                add_msg_board(sprintf('%s does not match %s',restraints.domain_start(1:poi1),restraints.domain_end(1:poi2))); 
                restrain=[];
                cancelled = true;
                return;
            else
                chain_model = restraints.domain_start(1:poi1);
            end
        else
            chain_model = '';
            poi1 = 0;
        end
    end
    res1 = str2double(restraints.domain_start(poi1+1:end));
    rese = str2double(restraints.domain_end(poi1+1:end));
    if isnan(res1) || isnan(rese),
        add_msg_board('ERROR: Residue numbers of domain could not be recognized.'); 
        restrain=[];
        cancelled = true;
        return;
    end
end

restraints.chain_model = chain_model;
restraints.res1 = res1;
restraints.rese = rese;

for k = res1:rese,
    adr = sprintf('%s%i',chain_model,k);
    [indices,message]=resolve_address(adr);
    if message.error ~=2 && message.error ~=13,
        add_msg_board(sprintf('ERROR: Residue %s does exist in the structure. Modelling impossible as this would create a clash.',adr)); 
        address = mk_address(indices,1);
        add_msg_board(address);
        restrain=[];
        cancelled = true;
        return;
    end
end

if isfield(restraints,'DEER'),
    poi = 0;
    llist = cell(0);
    label = cell(0);
    T = zeros(1,length(restraints.DEER));
    for k = 1:length(restraints.DEER),
        sep = strfind(restraints.DEER(k).label,'|');
        if isempty(sep)
            label1 = restraints.DEER(k).label;
            label2 = label1;
        else
            label1 = restraints.DEER(k).label(1:sep-1);
            label2 = restraints.DEER(k).label(sep+1:end);
        end
        [indices,message]=resolve_address(restraints.DEER(k).adr1);
        if message.error ~= 2 && message.error ~=13, % this site exists and is thus a beacon
            if message.error,
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr1)); 
                address = mk_address(indices,1);
                add_msg_board(address);
                restrain=[];
                cancelled = true;
                return;
            end
            poi = poi + 1;
            llist{poi} = restraints.DEER(k).adr1;
            T(poi) = restraints.DEER(k).T;
            label{poi} = label1;
            restraints.DEER(k).type1 = 1;
            restraints.DEER(k).indices1 = indices;
        else % this site is in the loop to be modelled
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,label1) || strcmpi(rotamer_libraries(kr).tc,label1),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel1 = NO;
                    restraints.DEER(k).type1 = 0;
                end
            end
            for kr = 1:length(ligand_libraries),
                if strcmpi(ligand_libraries(kr).label,label1) || strcmpi(ligand_libraries(kr).tc,label1),
                    Tvec = ligand_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,ligand_libraries(kr).files);
                    spin = get_relative_bilabel(libname,ligand_libraries(kr).V_ksi);
                    restraints.DEER(k).NO_rel1 = spin;
                    restraints.DEER(k).type1 = 0;
                end
            end
        end
        [indices,message]=resolve_address(restraints.DEER(k).adr2);
        if message.error ~= 2 && message.error ~= 13, % this site exists and is thus a beacon
            if message.error,
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr2)); 
                address = mk_address(indices,1);
                add_msg_board(address);
                restrain=[];
                cancelled = true;
                return;
            end
            poi = poi + 1;
            llist{poi} = restraints.DEER(k).adr2;
            T(poi) = restraints.DEER(k).T;
            label{poi} = label2;
            restraints.DEER(k).type2 = 1;
            restraints.DEER(k).indices2 = indices;
        else % this site is in the loop to be modelled
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,label2) || strcmpi(rotamer_libraries(kr).tc,label2)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel2 = NO;
                    restraints.DEER(k).type2 = 0;
                end
            end
            for kr = 1:length(ligand_libraries),
                if strcmpi(ligand_libraries(kr).label,label2) || strcmpi(ligand_libraries(kr).tc,label2),
                    Tvec = ligand_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.DEER(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,ligand_libraries(kr).files);
                    spin = get_relative_bilabel(libname,ligand_libraries(kr).V_ksi);
                    restraints.DEER(k).NO_rel2 = spin;
                    restraints.DEER(k).type2 = 0;
                end
            end
        end
    end
    labels = get_labels(llist,label,T);
    for k = 1:length(restraints.DEER),
        clabel = restraints.DEER(k).label;
        cT = restraints.DEER(k).T;
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 0, % internal restraint
            resa = separate_address(restraints.DEER(k).adr1);
            resb = separate_address(restraints.DEER(k).adr2);
            if isempty(resa) || isempty(resb),
                add_msg_board(sprintf('ERROR: Residue address %s or %s is invalid.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO1 = restraints.DEER(k).NO_rel1;
            NO2 = restraints.DEER(k).NO_rel2;
            if resa < resb,
                exch = resa; resa = resb; resb = exch;
                exch = NO1; NO1 = NO2; NO2 = exch;
            end
            if restraints.DEER(k).r ~= 0,
                [restrain,number,bnumber] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,bnumber,clabel,cT);
            else
                [monitor,number_monitor] = mk_internal_restraint(monitor,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,bnmon,clabel,cT);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 ==0, % beacon restraint, first residue beacon
            indices = restraints.DEER(k).indices1;
            resb = separate_address(restraints.DEER(k).adr1);
            for kl = 1:length(labels),
                if sum(abs(indices-labels(kl).indices)) == 0,
                    xyz_beacon = labels(kl).xyz;
                end
            end
            res_loop = separate_address(restraints.DEER(k).adr2);
            if isempty(res_loop),
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr2)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO = restraints.DEER(k).NO_rel2;
            if restraints.DEER(k).r ~= 0,
                [restrain,number,bnumber] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,bnumber,clabel,cT,indices,resb);
            else
                [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,bnmon,clabel,cT,indices,resb);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 1, % beacon restraint, second residue beacon
            indices = restraints.DEER(k).indices2;
            resb = separate_address(restraints.DEER(k).adr2);
            for kl = 1:length(labels),
                if sum(abs(indices-labels(kl).indices)) == 0,
                    xyz_beacon = labels(kl).xyz;
                end
            end
            res_loop = separate_address(restraints.DEER(k).adr1);
            if isempty(res_loop),
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr1)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO = restraints.DEER(k).NO_rel1;
            if restraints.DEER(k).r ~= 0,
               [restrain,number,bnumber] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,bnumber,clabel,cT,indices,resb);
            else
               [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,bnmon,clabel,cT,indices,resb);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 == 1, % nonsensical restraint inside defined structure
            add_msg_board(sprintf('Warning: Restraint between residues %s and %s inside defined structure will be ignored.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
        end
    end
end

if isfield(restraints,'depth'),
    for k = 1:length(restraints.depth),
        clabel = restraints.depth(k).label;
        cT = restraints.depth(k).T;
        if strcmpi(restraints.depth(k).label,'CA'),
            NO = [];
        else
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.depth(k).label),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.depth(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end
            end
        end
        res = separate_address(restraints.depth(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.depth(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        if restraints.depth(k).r ~= 0,
            [restrain,number,bnumber] = mk_depth_restraint(restrain,NO,res,restraints.depth(k).r,restraints.depth(k).sigr,res1,number,bnumber,clabel,cT);
        else
            [monitor,number_monitor] = mk_depth_restraint(monitor,NO,res,restraints.depth(k).r,restraints.depth(k).sigr,res1,number_monitor,bnmon,clabel,cT);
            number_monitor = number_monitor + 1;
        end
    end
end

if isfield(restraints,'oligomer'),
    for k = 1:length(restraints.oligomer),
        clabel = restraints.oligomer(k).label;
        cT = restraints.oligomer(k).T;
        if strcmpi(restraints.oligomer(k).label,'CA'),
            NO = [];
        else
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.oligomer(k).label),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-restraints.oligomer(k).T));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end
            end
            for kr = 1:length(ligand_libraries),
                if strcmpi(ligand_libraries(kr).label,restraints.oligomer(k).label) || strcmpi(ligand_libraries(kr).tc,restraints.oligomer(k).label),
                    libname = id2tag(1,ligand_libraries(kr).files);
                    NO = get_relative_bilabel(libname,ligand_libraries(kr).V_ksi);
                end
            end
        end
        res = separate_address(restraints.oligomer(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.oligomer(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        if restraints.oligomer(k).r ~= 0,
            [restrain,number,bnumber] = mk_oligomer_restraint(restrain,NO,res,restraints.oligomer(k).multi,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number,bnumber,clabel,cT);
        else
            [monitor,number_monitor] = mk_oligomer_restraint(monitor,NO,res,restraints.oligomer(k).multi,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number_monitor,bnmon,clabel,cT);
            number_monitor = number_monitor + 1;
        end
    end
end

if isfield(restraints,'cispeptides'),
    for k = 1:length(restraints.cispeptides),
        kr = restraints.cispeptides(k) - res1 + 1;
        restrain(kr).cis = 1;
    end
end

if isfield(restraints,'aprop'),
    for k = 1:length(restraints.aprop),
        res = separate_address(restraints.aprop(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.aprop(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).aprop = restraints.aprop(k).prop;
    end
end

if isfield(restraints,'bprop'),
    for k = 1:length(restraints.bprop),
        res = separate_address(restraints.bprop(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.bprop(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).bprop = restraints.bprop(k).prop;
    end
end

if isfield(restraints,'pprop'),
    for k = 1:length(restraints.pprop),
        res = separate_address(restraints.pprop(k).adr);
        if isempty(res),
            add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.pprop(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).pprop = restraints.pprop(k).prop;
    end
end

if isfield(restraints,'helices'),
    for k = 1:length(restraints.helices),
        hpoi = strfind(restraints.helices(k).adr,'-');
        ha = str2double(restraints.helices(k).adr(1:hpoi-1));
        he = str2double(restraints.helices(k).adr(hpoi+1:end));
        if isnan(ha) || isnan(he),
            add_msg_board(sprintf('ERROR: Wrong helix specification %s.',restraints.helices(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        for ks = 1:length(restraints.sequence),
            off1 = ks - 1 + res1 - ha;
            off2 = he - (ks - 1 + res1);
            if off1 >=0 && off2 >=0,
                restrain(ks).secondary = 1;
            end
            if off1 >=2 && off2 >=2,
                restrain(ks).secondary = 3;
            end
        end
    end
end

if isfield(restraints,'strands'),
    for k = 1:length(restraints.strands),
        hpoi = strfind(restraints.strands(k).adr,'-');
        ha = str2double(restraints.strands(k).adr(1:hpoi-1));
        he = str2double(restraints.strands(k).adr(hpoi+1:end));
        if isnan(ha) || isnan(he),
            add_msg_board(sprintf('ERROR: Wrong strand specification %s.',restraints.strands(k).adr)); 
            restrain=[];
            cancelled = true;
            return;
        end
        for ks = 1:length(restraints.sequence),
            off1 = ks - 1 + res1 - ha;
            off2 = he - (ks - 1 + res1);
            if off1 >=0 && off2 >=0,
                restrain(ks).secondary = 2;
            end
        end
    end
end

function labels = get_labels(llist,label,T)

global model
global hMain

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end

% check whether sites are already labelled and whether all restraint sites
% do exist
lindices=zeros(length(labels),4);
for k=1:length(labels),
    cindices=labels(k).indices;
    if ~isempty(cindices),
        lindices(k,:)=cindices;
    end
end
poi=0;
to_do_list{1}=' ';
for k=1:length(llist),
    adr1=llist{k};
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
        add_msg_board(sprintf('This site does not exist in current structure %s',mk_address(1)));
    end
    found=false;
    for l=1:length(labels),
        diff=ind1-lindices(l,:);
        if sum(abs(diff))==0,
            found=true;
        end
    end
    if ~found,
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}),
                found=true;
            end
        end
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end
    end
end

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label{k},T(k));
        hMain.store_undo=false;
        hMain.dynamic_rotamers=false;
        cmd(hMain,command);
    end
end

if isfield(model,'sites'),
    labels = label_information(model.sites);
else
    labels = cell(0);
end

function NO = get_relative_label(libname)

load(libname);
pops = rot_lib.calibration.pop;
pops = pops/sum(pops);
if isfield(rot_lib.usefull_atoms,'midNO')
    midNO = rot_lib.usefull_atoms.midNO;
    NO = zeros(1,3);
    for k = 1:length(rot_lib.library),
        coor = rot_lib.library(k).ecoor;
        NO = NO + pops(k)*(coor(midNO(1),2:4) + coor(midNO(2),2:4))/2;
    end
elseif isfield(rot_lib,'spin_density')
    NO = zeros(1,3);
    for k = 1:length(rot_lib.library),
        coor = rot_lib.library(k).ecoor;
        coor = coor(rot_lib.spin_density(:,1),2:4);
        coor = rot_lib.spin_density(:,2)'*coor;
        NO = NO + pops(k)*coor;
    end
end
NO = NO/(sum(pops));

function spin_coor = get_relative_bilabel(libname,V_ksi)

gas_un = 8.314472;    % universal gas constant in CI (J/(mol*K)       
T = 298;

load(libname);
pops = zeros(1,length(ligand_lib.library));
all_coor = zeros(length(ligand_lib.library),3);
if isfield(ligand_lib,'spin_density')
    spin_density = ligand_lib.spin_density;
else
    spin_density = [1 1];
end
sdp = spin_density(:,2)';
for k = 1:length(ligand_lib.library)
    ksi = pi*ligand_lib.library(k).dihedrals/pi;
    int_en = V_ksi*(1-cos(ksi))/2;
    pops(k) = exp(-int_en/(gas_un*T));
    all_coor(k,:) = sdp*ligand_lib.library(k).ecoor(spin_density(:,1),2:4);
end
spin_coor = pops*all_coor;
spin_coor = spin_coor/(sum(pops));

function labels=label_information(sites)

global model
global label_defs

poi=0;
for k0=1:length(sites),
    for k1=1:length(sites{k0}),
        for k=1:length(sites{k0}(k1).residue),
            poi=poi+1;
            labels(poi).indices=sites{k0}(k1).residue(k).indices;
            id=tag2id(sites{k0}(k1).residue(k).label,label_defs.restags);
            labels(poi).name=label_defs.residues(id).short_name;
            labels(poi).T=sites{k0}(k1).residue(k).T;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            x=sum(NOpos(:,1).*NOpos(:,4));
            y=sum(NOpos(:,2).*NOpos(:,4));
            z=sum(NOpos(:,3).*NOpos(:,4));
            labels(poi).xyz=[x y z];
            labels(poi).rmsd=NOpos_rmsd(NOpos);
        end
    end
end

function [rmsd,xyz]=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
xyz = [xmean,ymean,zmean];
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function [res,chain_model] = separate_address(adr)

spoi = strfind(adr,'|');
if isempty(spoi)
    cadr{1} = adr;
    res = 0;
else
    cadr{1} = adr(1:spoi-1);
    cadr{2} = adr(spoi+1:end);
    res = [0,0];
end
for k = 1:length(cadr)
    adr = cadr{k};
    poi = 0;
    poi1 = strfind(adr,']');
    if ~isempty(poi1),
        poi = poi1;
    end
    poi1 = strfind(adr,')');
    if ~isempty(poi1),
        poi = poi1;
    end
    poi1 = strfind(adr,'}');
    if ~isempty(poi1),
        poi = poi1;
    end
    if poi > 0,
        chain_model = adr(1:poi);
    else
        chain_model = '';
    end
    resstr = adr(poi+1:end);
    cres = str2double(resstr);
    if isnan(cres),
        res = [];
        return
    end
    if cres - floor(cres) > eps,
        res = [];
        return
    end
    res(k) = cres;
end

function [restrain,number,bnumber] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,rmean,sigr,res1,number,bnumber,label,T,bindices,resb)

scale_units = 10;

grace = 0.5; % 5 Å uncertainty of label position

k = res_loop - res1 + 1;
kr = length(restrain(k).r_beacon)+1;
restrain(k).label = NO;
restrain(k).r_beacon(kr).xyz = xyz_beacon;
restrain(k).r_beacon(kr).label_type = label;
restrain(k).r_beacon(kr).label_T = T;
restrain(k).r_beacon(kr).bindices = bindices;
restrain(k).r_beacon(kr).resb = resb;
if rmean > 0 && sigr > 0
    restrain(k).r_beacon(kr).type = 'Gaussian';
    restrain(k).r_beacon(kr).par1 = rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = sqrt(sigr^2 + grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_beacon(kr).type = 'bounds';
    restrain(k).r_beacon(kr).par1 = -rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end

function [restrain,number,bnumber] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,rmean,sigr,res1,number,bnumber,label,T)

scale_units = 10;

grace = 0.5; % 5 Å uncertainty of label position

if resa < resb, % restraint must be stored at the later site
    exch = resa;
    resa = resb;
    resb = exch;
    exch = NO1;
    NO1 = NO2;
    NO2 = exch;
end
k = resa - res1 + 1;
k2 = resb - res1 + 1;
kr = length(restrain(k).r_intern)+1;
restrain(k).label = NO1;
restrain(k2).label = NO2;
restrain(k).r_intern(kr).site = k2;
restrain(k).r_intern(kr).label_type = label;
restrain(k).r_intern(kr).label_T = T;
restrain(k).r_intern(kr).resb = resb;
if rmean > 0 && sigr > 0,
    restrain(k).r_intern(kr).type = 'Gaussian';
    restrain(k).r_intern(kr).par1 = rmean*scale_units;
    restrain(k).r_intern(kr).par2 = sqrt(sigr^2 + 2*grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_intern(kr).type = 'bounds';
    restrain(k).r_intern(kr).par1 = -rmean*scale_units;
    restrain(k).r_intern(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end

function [restrain,number,bnumber] = mk_depth_restraint(restrain,NO,res,rmean,sigr,res1,number,bnumber,label,T)

scale_units = 10;

k = res - res1 + 1;
kr = length(restrain(k).depth)+1;
restrain(k).depth(kr).label_type = label;
restrain(k).depth(kr).label_T = T;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).depth(kr).site = 'label';
else
    restrain(k).depth(kr).site = 'CA';
end
if rmean > 0 && sigr > 0,
    restrain(k).depth(kr).type = 'Gaussian';
    restrain(k).depth(kr).par1 = rmean*scale_units;
    restrain(k).depth(kr).par2 = sigr*scale_units;
    number = number + 1;
else
    restrain(k).depth(kr).type = 'bounds';
    restrain(k).depth(kr).par1 = -rmean*scale_units;
    restrain(k).depth(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end

function [restrain,number,bnumber] = mk_oligomer_restraint(restrain,NO,res,n,rmean,sigr,res1,number,bnumber,label,T)

scale_units= 10;

k = res(end) - res1 + 1;
kr = length(restrain(k).oligomer)+1;
restrain(k).oligomer(kr).res0 = res(1)  - res1 + 1;
restrain(k).oligomer(kr).label_type = label;
restrain(k).oligomer(kr).label_T = T;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).oligomer(kr).site = 'label';
else
    restrain(k).oligomer(kr).site = 'CA';
end
if rmean > 0 && sigr > 0,
    restrain(k).oligomer(kr).type = 'Gaussian';
    restrain(k).oligomer(kr).par1 = rmean*scale_units;
    restrain(k).oligomer(kr).par2 = sigr*scale_units;
    number = number + 1;
else
    restrain(k).oligomer(kr).type = 'bounds';
    restrain(k).oligomer(kr).par1 = -rmean*scale_units;
    restrain(k).oligomer(kr).par2 = -sigr*scale_units;
    bnumber = bnumber + 1;
end
restrain(k).oligomer(kr).n = n;

