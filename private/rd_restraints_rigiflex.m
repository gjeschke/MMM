function [restraints,failed] = rd_restraints_rigiflex(fname,unprocessed)
% function [restraints,failed] = rd_restraints_rigiflex(fname)
%
% reads mixed experimental restraints for a RigiFlex model from an ASCII file
% restraints are returned in a structure whose fields correspond to
% different restraint types
%
% the restraints are pre-processed with respect to the current structure in
% MMM, which must conform to the rigid body definitions
%
% flag unprocessed suppresses spin labelling and other processing, defaults
% to false
%
% currently implemented types are:
%
% MAXTRIALS maximum number of trials for exhaustive scan
% RIGID     rigid-body and reference point definitions
% DEER      distance distribution restraints between spin labels
% ODEER     distance distribution restraints between two different spin labels
% RLINK     RNA link (<= 7 Å per nt)
% PLINK     protein link (<= 3.8 Å per aa)
% CROSSLINK crosslinks, currently the linker type is stored, but not
%           further evaluated
% SANS      SANS data file and resolution parameter file
% SAXS      SAXS data file and resolution parameter file
% PFLEX     flexible peptide segment with residue numbers of first and last
%           residue, N-terminal and C-terminal anchor address (must refer
%           to a rigid body), use an asterisk (*) for terminal segments
% RFLEX     flexible RNA segment with residue numbers of first and last
%           residue, C5'-terminal and C3'-terminal anchor address (must
%           refer to a rigid body), use asterisk (*) for terminal segments
% 
%
% additionally, the following keywords are recognized
%
% ENSEMBLE      number of models in the final ensemble (first argument) and
%               probability threshold (second argument)
% UNITS         can be NM (nanometers, default) or A (Angstroem)
% MAXTIME       maximum time for a computation step, corresponds to the
%               full computation for rigid-body arrangements and model
%               assembly, but to a single linker in a single arrangement
%               for the flexible part, use FLEXTIME to specify the time for
%               flexible linkers differently
% FLEXTIME      maximum computation time for an ensemble of flexible
%               linkers/termini in a single rigid-body arrangement
% MODELS        number of rigid-body arrangements to be generated
% COLOR         color specifications for chains in the script file for
%               display of the rigid-body arrangement
% PDB           PDB file name of the rigid-body construct
% SUPERIMPOSE   which chain (of the rigid bodies) is superimposed in
%               ensembles
%
% G. Jeschke, 24.3.2017

global model

if ~exist('unprocessed','var') || isempty(unprocessed)
    unprocessed = false;
end

lpr_nt = 7; % maximum length per RNA nucleotide
lpr_aa = 3.8; % maximum length per amino acid residue

failed = true;

fid=fopen(fname);
if fid==-1
    add_msg_board('ERROR: Restraint file does not exist');
    return;
end

clear restraints

restraints.search = false;
restraints.maxtrials = -1;
restraints.ensemble = [];
restraints.p_model = [];
restraints.rb(1).chains = [];
restraints.rb(1).ref = [];
restraints.rb(1).indices = [];
restraints.rb(1).points = [];
restraints.rb(1).label = {''};
restraints.DEER(1).r = [];
restraints.DEER(1).sigr = [];
restraints.DEER(1).label = '';
restraints.DEER(1).ecoor1 = [];
restraints.DEER(1).indices1 = [];
restraints.DEER(1).ecoor2 = [];
restraints.DEER(1).indices2 = [];
restraints.links(1).maxr = [];
restraints.links(1).coor1 = [];
restraints.links(1).indices1 = [];
restraints.links(1).coor2 = [];
restraints.links(1).indices2 = [];
restraints.links(1).nta = [];
restraints.links(1).nte = [];
restraints.links(1).ref_indices = [];
restraints.xlinks(1).type = '';
restraints.xlinks(1).maxr = [];
restraints.xlinks(1).coor1 = [];
restraints.xlinks(1).indices1 = [];
restraints.xlinks(1).coor2 = [];
restraints.xlinks(1).indices2 = [];

restraints.stemloop_links = [];
restraints.SL_DEER = [];
restraints.SAXS = [];
restraints.SANS = [];
restraints.pflex = [];
restraints.rflex = [];
restraints.xlinks = [];
restraints.RNA = [];
restraints.max_time = [];
restraints.solutions = '';
restraints.solution_mode = '';
restraints.substitute = [];

DEER_poi=0;
rb_poi=0;
link_poi=0;
xlink_poi = 0;
SANS_poi = 0;
SAXS_poi = 0;
pflex_poi = 0;
rflex_poi = 0;
RNA_poi = 0;
max_chain = 0;
subspoi = 0;

label_adr = ':';
label_types = ':';
label_coor = zeros(1000,5);
label_indices = zeros(1000,4);
mode=0;
submode = 0;
scale_units=1;
stemlibs = {};
slpoi = 0;

nl = 0;
while 1
    tline = fgetl(fid);
    nl = nl + 1;
    if ~ischar(tline) || mode<0, break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            if k(1)>1
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end
        end
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#')
            switch upper(char(args(2)))
                case 'SEARCH'
                    mode=0;
                    restraints.search = true;
                case 'MAXTRIALS'
                    mode=0;
                    restraints.maxtrials=str2double(char(args(3)));
                case 'ENSEMBLE'
                    mode=0;
                    restraints.ensemble=str2double(char(args(3)));
                    if length(args) > 3
                        restraints.p_model = str2double(char(args(4)));
                    end
                case 'SOLUTIONS'
                    mode=0;
                    restraints.solutions = char(args(3));
                    if length(args) > 3
                        restraints.solution_mode = char(args(4));
                    else
                        restraints.solution_mode = 'single';
                    end
                case 'STEMLIB'
                    slpoi = slpoi + 1;
                    load(char(args(3)));
                    library.nta = str2double(char(args(4)));
                    library.nte = str2double(char(args(5)));
                    library.rrm = char(args(6));
                    motif = char(args(7));
                    if ~strcmpi(motif,library.chaintag)
                        add_msg_board(sprintf('ERROR: Stemloop library %s does not match binding motif in line %i. Aborting.',char(args(3)),nl));
                        fclose(fid);
                        return
                    end
                    ind1 = resolve_address(library.chaintag);
                    ind2 = resolve_address(library.rrm);
                    rba = 0;
                    for rb = 1:rb_poi
                        chains = restraints.rb(rb).chains;
                        if min(abs(chains-ind1(2))) == 0 && min(abs(chains-ind2(2))) == 0
                            rba = rb;
                        end
                    end
                    if rba == 0
                        add_msg_board(sprintf('ERROR: Rigid body for stemloop library %s not (yet) consistently defined in line %i. Aborting.',char(args(3)),nl));
                        fclose(fid);
                        return
                    end
                    library.rba = rba;
                    library.cind = ind1;
                    stemlibs{slpoi} = library;
                case 'DEER'
                    mode=1;
                    label = char(args(3));
                    label = check_label(label);
                    if isempty(label)
                        add_msg_board(sprintf('Warning: Label %s is unknown in line %i. Reverting to MTSL.',char(args(3))),nl);
                        label = 'R1A';
                    end
                case 'ODEER'
                    mode=10;
                    label1 = char(args(3));
                    label1 = check_label(label1);
                    if isempty(label1)
                        add_msg_board(sprintf('Warning: Label(1) %s is unknown in line %i. Reverting to MTSL.',char(args(3))),nl);
                        label = 'R1A';
                    end
                    label2 = char(args(4));
                    label2 = check_label(label2);
                    if isempty(label2)
                        add_msg_board(sprintf('Warning: Label(2) %s is unknown in line %i. Reverting to MTSL.',char(args(4))),nl);
                        label = 'R1A';
                    end
                case 'RIGID'
                    mode=2;
                    if length(args)<= 2
                        add_msg_board(sprintf('ERROR: Rigid-body definition in line %i misses chain identifiers. Aborting.',nl));
                        fclose(fid);
                        return
                    end
                    rb_poi = rb_poi + 1;
                    chains = zeros(1,length(args)-2);
                    for kc = 1:length(args)-2
                        indices = resolve_address(char(args(kc+2)));
                        if length(indices) ~= 2
                            add_msg_board(sprintf('ERROR: Rigid-body definition in line %i contains a wrong chain address. Aborting.',nl));
                            fclose(fid);
                            return
                        end
                        chains(kc) = indices(2);
                    end
                    restraints.rb(rb_poi).chains  = chains; 
                    max_chain = max([max_chain chains]);
                    restraints.rb(rb_poi).ref = zeros(3,5);
                    restraints.rb(rb_poi).indices = zeros(3,4);
                    ref_poi = 0;
                case 'RLINK'
                    mode = 3;
                case 'PLINK'
                    mode = 4;
                case 'CROSSLINK'
                    mode = 7;
                    linker = char(args(3));
                case 'UNITS'
                    mode=0;
                    unitstr=char(args(3));
                    switch upper(unitstr)
                        case {'A','ANGSTROEM'}
                            scale_units=0.1;
                        case {'NM','NANOMETER','NANOMETERS'}
                            scale_units=1;
                        otherwise
                            add_msg_board('Warning: Unknown unit identifier. Reverting to nanometers.');
                            scale_units=1;
                    end
                case 'BUILDTIME'
                    mode=0;
                    restraints.build_time = str2double(char(args(3))); 
                case 'FLEXTIME'
                    mode=0;
                    restraints.flex_time = str2double(char(args(3))); 
                case 'MAXTIME'
                    mode=0;
                    restraints.max_time = str2double(char(args(3))); 
                case 'MODELS'
                    mode=0;
                    restraints.models = str2double(char(args(3))); 
                case 'PROTEIN'
                    mode=0;
                    restraints.protein = char(args(3)); 
                case 'SANS'
                    mode=5;
                    SANS_poi = SANS_poi + 1;
                    restraints.SANS(SANS_poi).data = char(args(3)); 
                    if length(args) > 3
                        restraints.SANS(SANS_poi).illres = char(args(4)); 
                    else
                        restraints.SANS(SANS_poi).illres = ''; 
                    end
                    if length(args) > 4
                        restraints.SANS(SANS_poi).D2O = str2double(char(args(5))); 
                    else
                        restraints.SANS(SANS_poi).D2O = 0; 
                    end
                case 'SAXS'
                    mode=6;
                    SAXS_poi = SAXS_poi + 1;
                    restraints.SAXS(SAXS_poi).data = char(args(3));
                    if length(args) > 3
                        restraints.SAXS(SAXS_poi).sm = str2double(char(args(4))); 
                    else
                        restraints.SAXS(SAXS_poi).sm = []; 
                    end
                case 'PFLEX'
                    mode=8;
                    submode = 0;
                    pflex_poi = pflex_poi + 1;
                    restraints.pflex(pflex_poi).start = str2double(char(args(3)));
                    restraints.pflex(pflex_poi).end = str2double(char(args(4))); 
                    restraints.pflex(pflex_poi).Nanchor = char(args(5));
                    restraints.pflex(pflex_poi).Canchor = char(args(6));
                    restraints.pflex(pflex_poi).models = 20; % default setting
                    restraints.pflex(pflex_poi).time = 1; % default setting
                    restraints.pflex(pflex_poi).prob = 0.5; % default setting
                    pflex_DEER_poi = 0;
                    pflex_oligomer_poi = 0;
                    pflex_depth_poi = 0;
                    pflex_helix_poi = 0;
                    pflex_strand_poi = 0;
                    pflex_aprop_poi = 0;
                    pflex_bprop_poi = 0;
                    pflex_pprop_poi = 0;
                    if ~strcmp(restraints.pflex(pflex_poi).Nanchor,'*')
                        indices = resolve_address(restraints.pflex(pflex_poi).Nanchor);
                        if length(indices) ~= 4
                            add_msg_board(sprintf('Warning: N-terminal anchor address %s in line %i cannot be resolved.',restraints.pflex(pflex_poi).Nanchor,nl));
                            restraints.pflex(pflex_poi).Na_indices = [];
                        else
                            restraints.pflex(pflex_poi).Na_indices = indices;
                        end
                    else
                        restraints.pflex(pflex_poi).Na_indices = [];
                    end
                    if ~strcmp(restraints.pflex(pflex_poi).Canchor,'*')
                        indices = resolve_address(restraints.pflex(pflex_poi).Canchor);
                        if length(indices) ~= 4
                            add_msg_board(sprintf('Warning: C-terminal anchor address %s in line %i cannot be resolved.',restraints.pflex(pflex_poi).Canchor,nl));
                            restraints.pflex(pflex_poi).Ca_indices = [];
                        else
                            restraints.pflex(pflex_poi).Ca_indices = indices;
                        end
                    else
                        restraints.pflex(pflex_poi).Ca_indices = [];
                    end
                case 'RFLEX'
                    mode=9;
                    submode = 0;
                    rflex_poi = rflex_poi + 1;
                    restraints.rflex(rflex_poi).start = str2double(char(args(3)));
                    restraints.rflex(rflex_poi).end = str2double(char(args(4))); 
                    restraints.rflex(rflex_poi).anchor5p = char(args(5));
                    restraints.rflex(rflex_poi).anchor3p = char(args(6));
                    restraints.rflex(rflex_poi).models = 20; % default setting
                    restraints.rflex(rflex_poi).prob = 0.5; % default setting
                    rflex_DEER_poi = 0;
                    rflex_oligomer_poi = 0;
                    if ~strcmp(restraints.rflex(rflex_poi).anchor5p,'*')
                        indices = resolve_address(restraints.rflex(rflex_poi).anchor5p);
                        if length(indices) ~= 4
                            add_msg_board(sprintf('Warning: 5''-terminal anchor address %s in line %i cannot be resolved.',restraints.rflex(rflex_poi).anchor5p,nl));
                            restraints.rflex(rflex_poi).a5p_indices = [];
                        else
                            restraints.rflex(rflex_poi).a5p_indices = indices;
                        end
                    else
                        restraints.rflex(rflex_poi).a5p_indices = [];
                    end
                    if ~strcmp(restraints.rflex(rflex_poi).anchor3p,'*')
                        indices = resolve_address(restraints.rflex(rflex_poi).anchor3p);
                        if length(indices) ~= 4
                            add_msg_board(sprintf('Warning: 3''-terminal anchor address %s in line %i cannot be resolved.',restraints.rflex(rflex_poi).anchor3p,nl));
                            restraints.rflex(rflex_poi).a3p_indices = [];
                        else
                            restraints.rflex(rflex_poi).a3p_indices = indices;
                        end
                    else
                        restraints.rflex(rflex_poi).a3p_indices = [];
                    end
                case 'COLOR'
                    mode=0;
                    indices = resolve_address(char(args(3)));
                    if length(indices) ~= 2
                        add_msg_board(sprintf('Warning: Color definition in line %i contains a wrong chain address. Ignored.',nl));
                    else
                        restraints.color{indices(2)} = char(args(4));
                    end
                case 'SUPERIMPOSE'
                    mode=0;
                    indices = resolve_address(char(args(3)));
                    if length(indices) ~= 2
                        add_msg_board(sprintf('Warning: Superposition key has a wrong chain address %s in line %i. Stored as address.',char(args(3)),nl));
                        restraints.superimpose = char(args(3));
                    else
                        restraints.superimpose = indices(2);
                    end
                 case 'PDB'
                    mode=0;
                    pdbid = strtok(char(args(3)),':');
                    restraints.PDB = pdbid; % to conform to old format
                    if length(args) > 3
                        restraints.newID = char(args(4));
                    else
                        restraints.newID = 'RIFL';
                    end
                    if exist('model','var') && isfield(model,'structure_tags')
                        id = tag2id(restraints.PDB,model.structure_tags);
                    else
                        id = [];
                    end
                    if ~isempty(id)
                        snum=id;
                        model.current_structure=snum;
                    else
                        fname=get_pdb_file(restraints.PDB,true);
                        if isempty(fname)
                            fname = [restraints.PDB '.pdb'];
                        end
                        [message,snum]=add_pdb(fname);
                        if message.error
                            add_msg_board(sprintf('ERROR: Specified structure PDB file %s could not be retrieved from server',restraints.PDB));
                            add_msg_board(message.text);
                            return
                        else
                            model.current_structure=snum;
                        end
                    end
                case 'RNA'
                    mode=11;
                    submode = 0;
                    RNA_poi = RNA_poi + 1;
                    restraints.RNA(RNA_poi).models = 20; % default setting
                    restraints.RNA(RNA_poi).prob = 0.5; % default setting
                    restraints.RNA(RNA_poi).chain_id = char(args(3));
                    if length(args) > 3
                        restraints.RNA(RNA_poi).maxtime = str2double(char(args(4)));
                    else
                        restraints.RNA(RNA_poi).maxtime = 1;
                    end
                    RNA_DEER_poi = 0;
                    RNA_oligomer_poi = 0;
                    RNA_bind_poi = 0;
                    RNA_stem_poi = 0;
                case 'SUBSTITUTE'
                    mode=12;
                case 'RIGIFLEX'
                    mode = 0;
                case 'END'
                    mode=-1;
                otherwise
                    mode=0;
                    add_msg_board(sprintf('Warning: Unknown restraint mode %s',upper(char(args(2)))));
            end
        elseif mode>0 && ~strncmp(strtrim(char(args(1))),'%',1)
% %             keyboard
            switch mode
                case 1 % DEER
                    DEER_poi = DEER_poi + 1;
                    restraints.DEER(DEER_poi).r = 10*scale_units*str2double(char(args(3)));
                    restraints.DEER(DEER_poi).sigr = 10*scale_units*str2double(char(args(4)));
                    restraints.DEER(DEER_poi).label = label;
                    restraints.DEER(DEER_poi).adr1 = char(args(1));
                    restraints.DEER(DEER_poi).adr2 = char(args(2));
                    if length(args) > 4
                        restraints.DEER(DEER_poi).file = char(args(5));
                    end
                    if ~unprocessed
                        [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(char(args(1)),label,label_adr,label_types,label_coor,label_indices,stemlibs,restraints.rb);
                        if isempty(ecoor)
                            add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(1))));
                            fclose(fid);
                            return
                        end
                        restraints.DEER(DEER_poi).ecoor1 = ecoor;
                        restraints.DEER(DEER_poi).indices1 = indices;
                        [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(char(args(2)),label,label_adr,label_types,label_coor,label_indices,stemlibs,restraints.rb);
                        if isempty(ecoor)
                            add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(2))));
                            fclose(fid);
                            return
                        end
                        restraints.DEER(DEER_poi).ecoor2 = ecoor;
                        restraints.DEER(DEER_poi).indices2 = indices;
                    end
                case 10 % ODEER
                    DEER_poi = DEER_poi + 1;
                    if ~unprocessed
                        [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(char(args(1)),label1,label_adr,label_types,label_coor,label_indices,stemlibs,restraints.rb);
                        if isempty(ecoor)
                            add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(1))));
                            fclose(fid);
                            return
                        end
                        restraints.DEER(DEER_poi).ecoor1 = ecoor;
                        restraints.DEER(DEER_poi).indices1 = indices;
                    end
                    restraints.DEER(DEER_poi).r = 10*scale_units*str2double(char(args(3)));
                    restraints.DEER(DEER_poi).sigr = 10*scale_units*str2double(char(args(4)));
                    restraints.DEER(DEER_poi).label = [label1 '|' label2];
                    if length(args) > 4
                        restraints.DEER(DEER_poi).file = char(args(5));
                    end
                    if ~unprocessed
                        [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(char(args(2)),label2,label_adr,label_types,label_coor,label_indices,stemlibs,restraints.rb);
                        if isempty(ecoor)
                            add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(2))));
                            fclose(fid);
                            return
                        end
                        restraints.DEER(DEER_poi).ecoor2 = ecoor;
                        restraints.DEER(DEER_poi).indices2 = indices;
                    end
                    restraints.DEER(DEER_poi).adr1 = char(args(1));
                    restraints.DEER(DEER_poi).adr2 = char(args(2));
                case 2
                    adr = char(args(1));
                    label = char(args(2));
                    label = check_label(label);
                    if isempty(label)
                        add_msg_board(sprintf('Warning: Label %s is unknown in line %i. Reverting to MTSL.',char(args(3)),nl));
                        label = 'R1A';
                    end
                    ref_poi = ref_poi + 1;
                    if ~unprocessed
                        [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(adr,label,label_adr,label_types,label_coor,label_indices,[],restraints.rb);
                        if isempty(ecoor)
                            add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(2))));
                            fclose(fid);
                            return
                        end
                        if min(abs(indices(2)-restraints.rb(rb_poi).chains)) ~= 0
                            add_msg_board(sprintf('ERROR: Reference point %s is not in a chain of this rigid body in line %i. Aborting.',adr,nl));
                            fclose(fid);
                            return
                        end
                        restraints.rb(rb_poi).ref(ref_poi,:) = ecoor;
                        restraints.rb(rb_poi).indices(ref_poi,:) = indices;
                    end
                    restraints.rb(rb_poi).label{ref_poi} = label;
                    restraints.rb(rb_poi).points = ref_poi;
                case 3
                    adr1 = char(args(1));
                    adr2 = char(args(2));
                    resdiff = str2double(char(args(3)));
                    if length(args) > 3
                        maxlength = str2double(char(args(4)));
                    else
                        maxlength = lpr_nt*resdiff;
                    end
                    ind1 = [];
                    ind2 = [];
                    lib1 = 0;
                    lib2 = 0;
                    [msg,coor1] = get_object(sprintf('%s.C5''',adr1),'coor');
                    if msg.error
                        [coor1,ind1,lib1] = get_RNA_anchor(adr1,stemlibs);
                        if isempty(coor1)
                            add_msg_board(sprintf('ERROR: RNA anchor %s.C5'' coordinate could not be obtained in line %i. Aborting.\n',adr1,nl));
                            fclose(fid);
                            return
                        end
                    end
                    [msg,coor2] = get_object(sprintf('%s.C5''',adr2),'coor');
                    if msg.error
                        [coor2,ind2,lib2] = get_RNA_anchor(adr2,stemlibs);
                        if isempty(coor2)
                            add_msg_board(sprintf('ERROR: RNA anchor %s.C5'' coordinate could not be obtained in line %i. Aborting.\n',adr2,nl));
                            fclose(fid);
                            return
                        end
                    end
                    if isempty(ind1)
                        [ind1,message] = resolve_address(adr1);
                        if message.error
                            add_msg_board(sprintf('ERROR: RNA anchor %s indices not found in line %i. Aborting.\n',adr1,nl));
                            fclose(fid);
                            return
                        end
                    end
                    if isempty(ind2)
                        [ind2,message] = resolve_address(adr2);
                        if message.error
                            add_msg_board(sprintf('ERROR: RNA anchor %s indices not found in line %i. Aborting.\n',adr2,nl));
                            fclose(fid);
                            return
                        end
                    end
                    link_poi = link_poi + 1;
                    restraints.links(link_poi).coor1 = coor1;
                    restraints.links(link_poi).coor2 = coor2;
                    restraints.links(link_poi).maxr = maxlength;
                    restraints.links(link_poi).indices1 = ind1;
                    restraints.links(link_poi).indices2 = ind2;
                    restraints.links(link_poi).lib1 = lib1;
                    restraints.links(link_poi).lib2 = lib2;
                    restraints.links(link_poi).type = 'R';
                    npoi = strfind(adr1,')');
                    if isempty(npoi)
                        npoi = 1;
                    else
                        npoi = npoi+1;
                    end
                    restraints.links(link_poi).nta = str2double(adr1(npoi:end));
                    npoi = strfind(adr2,')');
                    if isempty(npoi)
                        npoi = 1;
                    else
                        npoi = npoi+1;
                    end
                    restraints.links(link_poi).nte = str2double(adr2(npoi:end));
                case 4
                    adr1 = char(args(1));
                    adr2 = char(args(2));
                    resdiff = str2double(char(args(3)));
                    [msg,coor1] = get_object(sprintf('%s.CA',adr1),'coor');
                    if msg.error
                        add_msg_board(sprintf('ERROR: Protein anchor %s.CA coordinate could not be obtained in line %i. Aborting.\n',adr1,nl));
                        fclose(fid);
                        return
                    end
                    [msg,coor2] = get_object(sprintf('%s.CA',adr2),'coor');
                    if msg.error
                        add_msg_board(sprintf('ERROR: Protein anchor %s.CA coordinate could not be obtained in line %i. Aborting.\n',adr2,nl));
                        fclose(fid);
                        return
                    end
                    [ind1,message] = resolve_address(adr1);
                    if message.error
                        add_msg_board(sprintf('ERROR: Protein anchor indices not found in line %i. Aborting.\n',adr1,nl));
                        fclose(fid);
                        return
                    end
                    [ind2,message] = resolve_address(adr2);
                    if message.error
                        add_msg_board(sprintf('ERROR: Protein anchor indices not found in line %i. Aborting.\n',adr2,nl));
                        fclose(fid);
                        return
                    end
                    link_poi = link_poi + 1;
                    restraints.links(link_poi).coor1 = coor1;
                    restraints.links(link_poi).coor2 = coor2;
                    restraints.links(link_poi).maxr = resdiff*lpr_aa;
                    restraints.links(link_poi).indices1 = ind1;
                    restraints.links(link_poi).indices2 = ind2;
                    restraints.links(link_poi).type = 'C5''';
                case 5
                    for kc = 1:length(args)
                        adr = char(args(kc));
                        if adr(1) == '!'
                            chains(kc) = -str2double(adr(2:end));
                        else
                            indices = resolve_address(adr);
                            if length(indices) ~= 2
                                add_msg_board(sprintf('ERROR: SANS definition in line %i contains a wrong chain address. Aborting.',nl));
                                fclose(fid);
                                return
                            end
                            chains(kc) = indices(2);
                        end
                    end
                    restraints.SANS(SANS_poi).chains  = chains;
                case 6
                    for kc = 1:length(args)
                        adr = char(args(kc));
                        if adr(1) == '!'
                            chains(kc) = -str2double(adr(2:end));
                        else
                            indices = resolve_address(adr);
                            if length(indices) ~= 2
                                add_msg_board(sprintf('ERROR: SAXS definition in line %i contains a wrong chain address. Aborting.',nl));
                                fclose(fid);
                                return
                            end
                            chains(kc) = indices(2);
                        end
                    end
                    restraints.SAXS(SAXS_poi).chains  = chains;
                case 7
                    adr1 = char(args(1));
                    adr2 = char(args(2));
                    [msg,coor1] = get_object(sprintf('%s.CA',adr1),'coor');
                    if msg.error
                        add_msg_board(sprintf('ERROR: Protein anchor %s.CA coordinate could not be obtained in line %i. Aborting.\n',adr1,nl));
                        fclose(fid);
                        return
                    end
                    [msg,coor2] = get_object(sprintf('%s.CA',adr2),'coor');
                    if msg.error
                        add_msg_board(sprintf('ERROR: Protein anchor %s.CA coordinate could not be obtained in line %i. Aborting.\n',adr2,nl));
                        fclose(fid);
                        return
                    end
                    [ind1,message] = resolve_address(adr1);
                    if message.error
                        add_msg_board(sprintf('ERROR: Protein anchor indices not found in line %i. Aborting.\n',adr1,nl));
                        fclose(fid);
                        return
                    end
                    [ind2,message] = resolve_address(adr2);
                    if message.error
                        add_msg_board(sprintf('ERROR: Protein anchor indices not found in line %i. Aborting.\n',adr2,nl));
                        fclose(fid);
                        return
                    end
                    xlink_poi = xlink_poi + 1;
                    restraints.xlinks(xlink_poi).adr1 = adr1;
                    restraints.xlinks(xlink_poi).adr2 = adr2;
                    restraints.xlinks(xlink_poi).coor1 = coor1;
                    restraints.xlinks(xlink_poi).coor2 = coor2;
                    restraints.xlinks(xlink_poi).maxr = 30;
                    restraints.xlinks(xlink_poi).indices1 = ind1;
                    restraints.xlinks(xlink_poi).indices2 = ind2;
                    restraints.xlinks(xlink_poi).type = linker;           
                case 8 % PFLEX
                    key = strtrim(char(args(1)));
                    if key(1) == ':' % new subkey
                        switch key
                            case ':SEQ'
                                restraints.pflex(pflex_poi).sequence = strtrim(char(args(2)));
                                submode = 0;
                            case ':MODELS'
                                restraints.pflex(pflex_poi).models = str2double((char(args(2))));
                                if length(args) > 2
                                    restraints.pflex(pflex_poi).prob = str2double(char(args(3)));
                                end
                                submode = 0;
                            case ':TIME'
                                restraints.pflex(pflex_poi).time = str2double((char(args(2))));
                                submode = 0;
                            case ':DEER'
                                pflex_label1 = char(args(2));
                                if length(args) > 2
                                    pflex_label2 = char(args(3));
                                else
                                    pflex_label2 = char(args(2));
                                end
                                submode = 1;
                            case ':OLIGOMER'
                                oligo_mult = str2double(char(args(2)));
                                oligo_label = char(args(3));
                                submode = 2;
                            case ':DEPTH'
                                depth_label = char(args(2));
                                submode = 3;
                            case ':HELICES'
                                submode = 4;
                            case ':STRANDS'
                                submode = 5;
                            case ':CISPEPTIDE'
                                cpp  = 0;
                                cispeptides = zeros(1,1000);
                                for ka = 2:length(args)
                                    cpp = cpp + 1;
                                    cispeptides(cpp)=str2double(char(args(ka)));
                                end
                                restraints.pflex(pflex_poi).cispeptides = cispeptides(1:cpp);
                                submode = 0;
                            case ':APROP'
                                submode = 6;
                            case ':BPROP'
                                submode = 7;
                            case ':PPROP'
                                submode = 8;
                            otherwise
                                add_msg_board(sprintf('Warning: Unknown subkey %s encountered. Ignored.',key));
                                submode = 0;
                        end
                    elseif submode>0 && ~strncmp(strtrim(char(args(1))),'%',1)
                        switch submode
                            case 1 % DEER
                                pflex_DEER_poi = pflex_DEER_poi + 1;
                                restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).label1 = pflex_label1;
                                restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).label2 = pflex_label2;
                                restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).adr1 = char(args(1));
                                restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).adr2 = char(args(2));
                                restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).r = 10*scale_units*str2double(char(args(3)));
                                restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).sigr = 10*scale_units*str2double(char(args(4)));
                                if length(args) > 4
                                    restraints.pflex(pflex_poi).DEER(pflex_DEER_poi).file = char(args(5));
                                end
                            case 2 % OLIGOMER
                                pflex_oligomer_poi = pflex_oligomer_poi + 1;
                                restraints.pflex(pflex_poi).oligomer(pflex_oligomer_poi).label = oligo_label;
                                restraints.pflex(pflex_poi).oligomer(pflex_oligomer_poi).mult = oligo_mult;
                                restraints.pflex(pflex_poi).oligomer(pflex_oligomer_poi).num = str2double(char(args(1)));
                                restraints.pflex(pflex_poi).oligomer(pflex_oligomer_poi).r = 10*scale_units*str2double(char(args(2)));
                                restraints.pflex(pflex_poi).oligomer(pflex_oligomer_poi).sigr = 10*scale_units*str2double(char(args(3)));
                            case 3 % DEPTH
                                pflex_depth_poi = pflex_depth_poi + 1;
                                restraints.pflex(pflex_poi).depth(pflex_depth_poi).label = depth_label;
                                restraints.pflex(pflex_poi).depth(pflex_depth_poi).num = str2double(char(args(1)));
                                restraints.pflex(pflex_poi).depth(pflex_depth_poi).z = 10*scale_units*str2double(char(args(2)));
                                restraints.pflex(pflex_poi).depth(pflex_depth_poi).sigz = 10*scale_units*str2double(char(args(3)));
                            case 4 % HELICES
                                adr = char(args(1));
                                hp = strfind(adr,'-');
                                if isempty(hp)
                                    add_msg_board(sprintf('Warning: Wrong address format (%s) for HELICES restraint in line %i. Ignored.\n',adr,nl));
                                else
                                    pflex_helix_poi = pflex_helix_poi + 1;
                                    restraints.pflex(pflex_poi).helices(pflex_helix_poi).Nterm = str2double(adr(1:hp-1));
                                    restraints.pflex(pflex_poi).helices(pflex_helix_poi).Cterm = str2double(adr(hp+1:end));
                                end
                            case 5 % STRANDS
                                adr = char(args(1));
                                hp = strfind(adr,'-');
                                if isempty(hp)
                                    add_msg_board(sprintf('Warning: Wrong address format (%s) for STRANDS restraint in line %i. Ignored.\n',adr,nl));
                                else
                                    pflex_strand_poi = pflex_strand_poi + 1;
                                    restraints.pflex(pflex_poi).strands(pflex_strand_poi).Nterm = str2double(adr(1:hp-1));
                                    restraints.pflex(pflex_poi).strands(pflex_strand_poi).Cterm = str2double(adr(hp+1:end));
                                end
                            case 6 % APROP
                                pflex_aprop_poi = pflex_aprop_poi+1;
                                restraints.pflex(pflex_poi).aprop(pflex_aprop_poi).num = str2double(char(args(1)));
                                restraints.pflex(pflex_poi).aprop(pflex_aprop_poi).prop = str2double(char(args(2)));
                            case 7 % BPROP
                                pflex_bprop_poi = pflex_bprop_poi+1;
                                restraints.pflex(pflex_poi).bprop(pflex_bprop_poi).num = str2double(char(args(1)));
                                restraints.pflex(pflex_poi).bprop(pflex_bprop_poi).prop = str2double(char(args(2)));
                            case 8 % PPROP
                                pflex_pprop_poi = pflex_pprop_poi+1;
                                restraints.pflex(pflex_poi).pprop(pflex_pprop_poi).num = str2double(char(args(1)));
                                restraints.pflex(pflex_poi).pprop(pflex_pprop_poi).prop = str2double(char(args(2)));
                        end
                    end
                case 9 % RFLEX
                    key = strtrim(char(args(1)));
                    if key(1) == ':' % new subkey
                        switch key
                            case ':SEQ'
                                restraints.rflex(rflex_poi).sequence = strtrim(char(args(2)));
                                submode = 0;
                            case ':MODELS'
                                restraints.rflex(rflex_poi).models = str2double((char(args(2))));
                                if length(args) > 2
                                    restraints.rflex(rflex_poi).prob = str2double(char(args(3)));
                                end
                                submode = 0;
                            case ':DEER'
                                rflex_label1 = char(args(2));
                                if length(args) > 2
                                    rflex_label2 = char(args(3));
                                else
                                    rflex_label2 = char(args(2));
                                end
                                submode = 1;
                            case ':OLIGOMER'
                                oligo_mult = str2double(char(args(2)));
                                oligo_label = char(args(3));
                                submode = 2;
                        end
                    elseif submode>0 && ~strncmp(strtrim(char(args(1))),'%',1)
                        switch submode
                            case 1 % DEER
                                rflex_DEER_poi = rflex_DEER_poi + 1;
                                restraints.rflex(rflex_poi).DEER(rflex_DEER_poi).label1 = rflex_label1;
                                restraints.rflex(rflex_poi).DEER(rflex_DEER_poi).label2 = rflex_label2;
                                restraints.rflex(rflex_poi).DEER(rflex_DEER_poi).adr1 = char(args(1));
                                restraints.rflex(rflex_poi).DEER(rflex_DEER_poi).adr2 = char(args(2));
                                restraints.rflex(rflex_poi).DEER(rflex_DEER_poi).r = 10*scale_units*str2double(char(args(3)));
                                restraints.rflex(rflex_poi).DEER(rflex_DEER_poi).sigr = 10*scale_units*str2double(char(args(4)));
                            case 2 % OLIGOMER
                                rflex_oligomer_poi = rflex_oligomer_poi + 1;
                                restraints.rflex(rflex_poi).oligomer(rflex_oligomer_poi).label = oligo_label;
                                restraints.rflex(rflex_poi).oligomer(rflex_oligomer_poi).mult = oligo_mult;
                                restraints.rflex(rflex_poi).oligomer(rflex_oligomer_poi).adr = char(args(1));
                                restraints.rflex(rflex_poi).oligomer(rflex_oligomer_poi).r = 10*scale_units*str2double(char(args(2)));
                                restraints.rflex(rflex_poi).oligomer(rflex_oligomer_poi).sigr = 10*scale_units*str2double(char(args(3)));
                        end
                    end
                case 11 % RNA
                    key = strtrim(char(args(1)));
                    if key(1) == ':' % new subkey
                        switch key
                            case ':SEQ'
                                restraints.RNA(RNA_poi).nta = str2double((char(args(2))));
                                restraints.RNA(RNA_poi).nte = str2double((char(args(3))));
                                restraints.RNA(RNA_poi).sequence = strtrim(char(args(4)));
                                submode = 0;
                            case ':MODELS'
                                restraints.RNA(RNA_poi).models = str2double((char(args(2))));
                                if length(args) > 2
                                    restraints.RNA(RNA_poi).prob = str2double(char(args(3)));
                                end
                                submode = 0;
                            case {':DEER',':ODEER'}
                                RNA_label1 = char(args(2));
                                if length(args) > 2
                                    RNA_label2 = char(args(3));
                                else
                                    RNA_label2 = char(args(2));
                                end
                                submode = 1;
                            case ':OLIGOMER'
                                oligo_mult = str2double(char(args(2)));
                                oligo_label = char(args(3));
                                submode = 2;
                            case ':STEM'
                                submode = 4;
                            case ':BIND'
                                submode = 3;
                            otherwise
                                add_msg_board(sprintf('Warning: Unknown subkey %s encountered in #RNA. Ignored.',key));
                                submode = 0;
                        end
                    elseif submode>0 && ~strncmp(strtrim(char(args(1))),'%',1)
                        switch submode
                            case 1 % DEER
                                RNA_DEER_poi = RNA_DEER_poi + 1;
                                restraints.RNA(RNA_poi).DEER(RNA_DEER_poi).label1 = RNA_label1;
                                restraints.RNA(RNA_poi).DEER(RNA_DEER_poi).label2 = RNA_label2;
                                restraints.RNA(RNA_poi).DEER(RNA_DEER_poi).adr1 = char(args(1));
                                restraints.RNA(RNA_poi).DEER(RNA_DEER_poi).adr2 = char(args(2));
                                restraints.RNA(RNA_poi).DEER(RNA_DEER_poi).r = 10*scale_units*str2double(char(args(3)));
                                restraints.RNA(RNA_poi).DEER(RNA_DEER_poi).sigr = 10*scale_units*str2double(char(args(4)));
                            case 2 % OLIGOMER
                                RNA_oligomer_poi = RNA_oligomer_poi + 1;
                                restraints.RNA(RNA_poi).oligomer(RNA_oligomer_poi).label = oligo_label;
                                restraints.RNA(RNA_poi).oligomer(RNA_oligomer_poi).mult = oligo_mult;
                                adr = char(args(1));
                                if adr(1) == 'R'
                                    adr = adr(2:end);
                                end
                                restraints.RNA(RNA_poi).oligomer(RNA_oligomer_poi).num = str2double(adr);
                                restraints.RNA(RNA_poi).oligomer(RNA_oligomer_poi).r = 10*scale_units*str2double(char(args(2)));
                                restraints.RNA(RNA_poi).oligomer(RNA_oligomer_poi).sigr = 10*scale_units*str2double(char(args(3)));
                            case 3 % BIND
                                RNA_bind_poi = RNA_bind_poi + 1;
                                restraints.RNA(RNA_poi).bind(RNA_bind_poi).nta = str2double(char(args(1)));
                                restraints.RNA(RNA_poi).bind(RNA_bind_poi).anchora = char(args(2));
                                restraints.RNA(RNA_poi).bind(RNA_bind_poi).nte = str2double(char(args(3)));
                                restraints.RNA(RNA_poi).bind(RNA_bind_poi).anchore = char(args(4));
                            case 4 % STEM
                                RNA_stem_poi = RNA_stem_poi + 1;
                                restraints.RNA(RNA_poi).stems(RNA_stem_poi).C5pf = str2double(char(args(1)));
                                restraints.RNA(RNA_poi).stems(RNA_stem_poi).C3pf = str2double(char(args(2)));
                                restraints.RNA(RNA_poi).stems(RNA_stem_poi).C5pb = str2double(char(args(3)));
                                restraints.RNA(RNA_poi).stems(RNA_stem_poi).C3pb = str2double(char(args(4)));
                                restraints.RNA(RNA_poi).stems(RNA_stem_poi).lib = char(args(5));
                        end
                    end
                        case 12 % SUBSTITUTE
                            subspoi = subspoi + 1;
                            for karg = 1:length(args)
                                restraints.substitute(subspoi).chains{karg} = char(args(karg));
                            end
            end
        end
    end
end

fclose(fid);

restraints.stemlibs = stemlibs;

if unprocessed
    failed = false;
    return
end

% Assemble the points defined in each rigid body and make separated
% restraint lists

restraints.points = cell(1,rb_poi);
restraints.pindices = cell(1,rb_poi);
restraints.plabels = cell(1,rb_poi);
restraints.dmat0 = zeros(3*rb_poi);
restraints.auxiliary = zeros(1000,6);
restraints.stemloop_rlinks(1).rba = [];
restraints.stemloop_rlinks(1).coor1 = [];
restraints.stemloop_rlinks(1).coor2 = [];
restraints.stemloop_rlinks(1).lib = [];

aux_poi = 0;
core_poi = 0;
sll = 0;
% Assemble reference points and make intra-rigid-body part of the distance
% matrix
rb_pointers = 3*ones(1,rb_poi);
for kr = 1:rb_poi % 
    coor = zeros(500,3);
    indices = zeros(500,4);
    if restraints.rb(kr).points ~= 3
        add_msg_board(sprintf('ERROR: Rigid body %i has %i reference points, but should have 3. Aborting.\n',kr,restraints.rb(kr).points));
        fclose(fid);
        return
    end;
    for kp = 1:3
        coor(kp,:) = restraints.rb(kr).ref(kp,2:4);
        indices(kp,:) = restraints.rb(kr).indices(kp,:);
    end
    restraints.points{kr} = coor;
    restraints.pindices{kr} = indices;
    restraints.plabels{kr} = restraints.rb(kr).label;
    bas = (kr-1)*3;
    for k1 = 1:2
        for k2 = k1:3
            r = norm(coor(k1,:)-coor(k2,:));
            restraints.dmat0(bas+k1,bas+k2) = r;
            restraints.dmat0(bas+k2,bas+k1) = r;
        end
    end
end
restraints.lb = restraints.dmat0; % initialize lower bounds, reference point distances are strictly preserved
restraints.ub = restraints.dmat0; % initialize upper bounds, reference point distances are strictly preserved

% Separate distance restraints into those between reference points,
% auxiliary ones, and delayed restraints that require stemloop attachment,
% add the former ones to the distance matrix

sl_deer = 0;
for kr = 1:DEER_poi
    label1 = restraints.DEER(kr).label;
    label2 = label1;
    sep = strfind(restraints.DEER(kr).label,'|');
    if ~isempty(sep)
        label1 = restraints.DEER(kr).label(1:sep-1);
        label2 = restraints.DEER(kr).label(sep+1:end);
    end
    if isempty(restraints.DEER(kr).indices1) || isempty(restraints.DEER(kr).indices2) % delayed DEER restraint
        sl_deer = sl_deer+1;
        known = [0,0];
        if ~isempty(restraints.DEER(kr).indices1)
            [krb1,kp1] = identify(restraints.DEER(kr).indices1,label1,restraints);
            if krb1 == 0
                wadr = mk_address(restraints.DEER(kr).indices1);
                add_msg_board(sprintf('Warning: Label %s is not in a rigid body. Restraint will be ignored in this step.\n',wadr));
                continue
            end
            known(1) = 1;
        end
        if ~isempty(restraints.DEER(kr).indices2)
            [krb2,kp2] = identify(restraints.DEER(kr).indices2,label2,restraints);
            if krb2 == 0
                wadr = mk_address(restraints.DEER(kr).indices2);
                add_msg_board(sprintf('Warning: Label %s is not in a rigid body. Restraint will be ignored in this step.\n',wadr));
                continue
            end
            known(2) = 1;
        end
        restraints.SL_DEER(sl_deer).r = restraints.DEER(kr).r;
        restraints.SL_DEER(sl_deer).sigr = restraints.DEER(kr).sigr;
        restraints.SL_DEER(sl_deer).coor1 = restraints.DEER(kr).ecoor1(:,2:4);
        restraints.SL_DEER(sl_deer).sig1 = restraints.DEER(kr).ecoor1(:,5);
        restraints.SL_DEER(sl_deer).coor2 = restraints.DEER(kr).ecoor2(:,2:4);
        restraints.SL_DEER(sl_deer).sig2 = restraints.DEER(kr).ecoor2(:,5);
        restraints.SL_DEER(sl_deer).rba = [restraints.DEER(kr).ecoor1(1,1) restraints.DEER(kr).ecoor2(1,1)];
        restraints.SL_DEER(sl_deer).stemlib = [0,0];
        if ~known(1)
            for lib = 1:length(restraints.stemlibs)
                sllib = restraints.stemlibs{lib};
                if sllib.rba == restraints.SL_DEER(sl_deer).rba(1)
                    restraints.SL_DEER(sl_deer).stemlib(1) = lib;
                end
            end
        end
        if ~known(2)
            for lib = 1:length(restraints.stemlibs)
                sllib = restraints.stemlibs{lib};
                if sllib.rba == restraints.SL_DEER(sl_deer).rba(2)
                    restraints.SL_DEER(sl_deer).stemlib(2) = lib;
                end
            end
        end
        continue
    end
    [krb1,kp1] = identify(restraints.DEER(kr).indices1,label1,restraints);
    if krb1 == 0
        wadr = mk_address(restraints.DEER(kr).indices1);
        add_msg_board(sprintf('Warning: Label %s is not in a rigid body. Restraint will be ignored in this step.\n',wadr));
    else
        k1 = 0;
        if kp1 > 0 && kp1 <= 3
            k1 = 3*(krb1-1)+kp1;
        end
        [krb2,kp2] = identify(restraints.DEER(kr).indices2,label2,restraints);
        if krb2 == 0
            wadr = mk_address(restraints.DEER(kr).indices2);
            add_msg_board(sprintf('Warning: Label %s is not in a rigid body. Restraint will be ignored in this step.\n',wadr));
        else
            k2 = 0;
            if kp2 > 0 && kp2 <= 3
                k2 = 3*(krb2-1)+kp2;
            end
            if k1 ~= 0 && k2 ~= 0 % both residues are reference points
                if krb1 == krb2
                    add_msg_board('Warning: DEER restraint within the same rigid body will be ignored');
                else
                     restraints.dmat0(k1,k2) = restraints.DEER(kr).r;
                     restraints.dmat0(k2,k1) = restraints.DEER(kr).r;
                     restraints.lb(k1,k2) = restraints.DEER(kr).r - restraints.DEER(kr).sigr;
                     restraints.lb(k2,k1) = restraints.DEER(kr).r - restraints.DEER(kr).sigr;
                     restraints.ub(k1,k2) = restraints.DEER(kr).r + restraints.DEER(kr).sigr;
                     restraints.ub(k2,k1) = restraints.DEER(kr).r + restraints.DEER(kr).sigr;
                     core_poi = core_poi + 1;
                     restraints.core(core_poi,:) = [krb1 kp1 krb2 kp2 restraints.DEER(kr).r restraints.DEER(kr).sigr];
                end;
            else % this is an auxiliary restraint
                if kp1 == 0
                    kp1 = rb_pointers(krb1) + 1;
                    rb_pointers(krb1) = kp1;
                    coor = restraints.points{krb1};
                    coor(kp1,:) = restraints.DEER(kr).ecoor1(2:4);
                    restraints.points{krb1} = coor;
                    indices = restraints.pindices{krb1};
                    indices(kp1,:) = restraints.DEER(kr).indices1;
                    restraints.pindices{krb1} = indices;
                    labels = restraints.plabels{krb1};
                    labels{kp1} = label1;
                    restraints.plabels{krb1} = labels;
                end
                if kp2 == 0
                    kp2 = rb_pointers(krb2) + 1;
                    rb_pointers(krb2) = kp2;
                    coor = restraints.points{krb2};
                    coor(kp2,:) = restraints.DEER(kr).ecoor2(2:4);
                    restraints.points{krb2} = coor;
                    indices = restraints.pindices{krb2};
                    indices(kp2,:) = restraints.DEER(kr).indices2;
                    restraints.pindices{krb2} = indices;
                    labels = restraints.plabels{krb2};
                    labels{kp2} = label2;
                    restraints.plabels{krb2} = labels;
                end
                aux_poi = aux_poi + 1;
                restraints.auxiliary(aux_poi,:) = [krb1 kp1 krb2 kp2 restraints.DEER(kr).r restraints.DEER(kr).sigr];
            end
        end
    end
end
restraints.auxiliary = restraints.auxiliary(1:aux_poi,:);
% Assign link restraints to rigid bodies
for kl = 1:link_poi
    if length(restraints.links(kl).indices1) == 1 || length(restraints.links(kl).indices2) == 1 % delayed link restraint
        sll = sll + 1;
        if length(restraints.links(kl).indices1) > 1
            krb1 = identify(restraints.links(kl).indices1,restraints.links(kl).type,restraints);
        else
            krb1 = -restraints.links(kl).indices1;
        end
        if length(restraints.links(kl).indices2) > 1
            krb2 = identify(restraints.links(kl).indices2,restraints.links(kl).type,restraints);
        else
            krb2 = -restraints.links(kl).indices2;
        end
        restraints.stemloop_links(sll).rba = [krb1,krb2];
        restraints.stemloop_links(sll).coor1 = restraints.links(kl).coor1;
        restraints.stemloop_links(sll).coor2 = restraints.links(kl).coor2;
        restraints.stemloop_links(sll).lib = [restraints.links(kl).lib1,restraints.links(kl).lib2];
        restraints.stemloop_links(sll).maxr = restraints.links(kl).maxr;
        restraints.stemloop_links(sll).nta = restraints.links(kl).nta;
        restraints.stemloop_links(sll).nte = restraints.links(kl).nte;
        continue
    end
    % link restraints cannot yet be in the point list
    krb1 = identify(restraints.links(kl).indices1,restraints.links(kl).type,restraints);
    krb2 = identify(restraints.links(kl).indices2,restraints.links(kl).type,restraints);
    kp1 = rb_pointers(krb1) + 1;
    rb_pointers(krb1) = kp1;
    coor = restraints.points{krb1};
    coor(kp1,:) = restraints.links(kl).coor1;
    restraints.points{krb1} = coor;
    indices = restraints.pindices{krb1};
    indices(kp1,:) = restraints.links(kl).indices1;
    restraints.pindices{krb1} = indices;
    labels = restraints.plabels{krb1};
    labels{kp1} = restraints.links(kl).type;
    restraints.plabels{krb1} = labels;
    kp2 = rb_pointers(krb2) + 1;
    rb_pointers(krb2) = kp2;
    coor = restraints.points{krb2};
    coor(kp2,:) = restraints.links(kl).coor2;
    restraints.points{krb2} = coor;
    indices = restraints.pindices{krb2};
    indices(kp2,:) = restraints.links(kl).indices2;
    restraints.pindices{krb2} = indices;
    labels = restraints.plabels{krb2};
    labels{kp2} = restraints.links(kl).type;
    restraints.plabels{krb2} = labels;
    restraints.links(kl).ref_indices = [krb1 kp1 krb2 kp2];
end

% Assign crosslink restraints to rigid bodies
for kl = 1:xlink_poi
    [krb1,kp1] = identify(restraints.xlinks(kl).indices1,restraints.xlinks(kl).type,restraints);
    [krb2,kp2] = identify(restraints.xlinks(kl).indices2,restraints.xlinks(kl).type,restraints);
    if kp1 == 0
        kp1 = rb_pointers(krb1) + 1;
        rb_pointers(krb1) = kp1;
        coor = restraints.points{krb1};
        coor(kp1,:) = restraints.xlinks(kl).coor1;
        restraints.points{krb1} = coor;
        indices = restraints.pindices{krb1};
        indices(kp1,:) = restraints.xlinks(kl).indices1;
        restraints.pindices{krb1} = indices;
        labels = restraints.plabels{krb1};
        labels{kp1} = restraints.xlinks(kl).type;
        restraints.plabels{krb1} = labels;
    end
    if kp2 == 0
        kp2 = rb_pointers(krb2) + 1;
        rb_pointers(krb2) = kp2;
        coor = restraints.points{krb2};
        coor(kp2,:) = restraints.xlinks(kl).coor2;
        restraints.points{krb2} = coor;
        indices = restraints.pindices{krb2};
        indices(kp2,:) = restraints.xlinks(kl).indices2;
        restraints.pindices{krb2} = indices;
        labels = restraints.plabels{krb2};
        labels{kp2} = restraints.xlinks(kl).type;
        restraints.plabels{krb2} = labels;
    end
   restraints.xlinks(kl).ref_indices = [krb1 kp1 krb2 kp2];
end
for krb = 1:length(restraints.points)
    points = restraints.points{krb};
    points = points(1:rb_pointers(krb),:);
    restraints.points{krb} = points;
    pindices = restraints.pindices{krb};
    pindices = pindices(1:rb_pointers(krb),:);
    restraints.pindices{krb} = pindices;
end;

% Make topology of chain connections by flexible linkers
% peptide part

peptide_tags = cell(1,20);
pep_chains = 0;
for kl = 1:length(restraints.pflex)
    % identify the chains where this linker is anchored
    ctag1 = get_ctag(restraints.pflex(kl).Nanchor); 
    ctag2 = get_ctag(restraints.pflex(kl).Canchor);
    % identify the contiguous chain to which this linker belongs
    pc = [];
    for kc = 1:pep_chains
        id = tag2id(ctag1,peptide_tags{kc});
        if ~isempty(id)
            if ~isempty(pc)
                add_msg_board(sprintf('ERROR: Peptide linker %i belongs to two contiguous chains by linking to chain %s. Aborting.\n',kl,ctag1));
                fclose(fid);
                return
            else
                pc = kc;
            end
        end
        id = tag2id(ctag2,peptide_tags{kc});
        if ~isempty(id)
            if ~isempty(pc)
                add_msg_board(sprintf('ERROR: Peptide linker %i belongs to two contiguous chains by linking to chains %s and %s. Aborting.\n',kl,ctag1,ctag2));
                fclose(fid);
                return
            else
                pc = kc;
            end
        end
    end
    if isempty(pc) % this is a new contiguous chain
        pep_chains = pep_chains + 1;
        nctags = sprintf(':%i:',kl);
        if ~isempty(ctag1)
            nctags = sprintf(':%s%s',ctag1,nctags);
        end
        if ~isempty(ctag2)
            nctags = sprintf('%s%s:',nctags,ctag2);
        end
        peptide_tags{pep_chains} = nctags;
    else % the linker must be inserted at the right place
        ctags = peptide_tags{pc};
        Napoi = strfind(ctags,ctag1);
        Capoi = strfind(ctags,ctag2);
        if ~isempty(Napoi) && ~isempty(Capoi)
            nctags = sprintf('%s:%i:%s',ctags(1:Napoi),kl,ctags(Capoi:end));
        elseif ~isempty(Napoi)
            nctags = sprintf('%s:%i:',ctags(1:Napoi),kl);
            if ~isempty(ctag2)
                nctags = sprintf('%s%s:',nctags,ctag2);
            end
        else
            nctags = sprintf(':%i:%s',kl,ctags(Capoi:end));
            if ~isempty(ctag1)
                nctags = sprintf(':%s%s',ctag1,nctags);
            end
        end
        peptide_tags{pc} = nctags;
    end
end
restraints.peptide_tags = peptide_tags(1:pep_chains);

% Make topology of chain connections by flexible linkers
% RNA part

% old version with # RFLEX
% RNA_tags = cell(1,20);
% RNA_chains = 0;
% for kl = 1:length(restraints.rflex)
%     klr = kl + length(restraints.pflex);
%     % identify the chains where this linker is anchored
%     ctag1 = get_ctag(restraints.rflex(kl).anchor5p); 
%     ctag2 = get_ctag(restraints.rflex(kl).anchor3p);
%     % identify the contiguous chain to which this linker belongs
%     pc = [];
%     for kc = 1:RNA_chains
%         id = tag2id(ctag1,RNA_tags{kc});
%         if ~isempty(id)
%             if ~isempty(pc)
%                 add_msg_board(sprintf('ERROR: RNA linker %i belongs to two contiguous chains by linking to chain %s. Aborting.\n',kl,ctag1));
%                 fclose(fid);
%                 return
%             else
%                 pc = kc;
%             end
%         end
%         id = tag2id(ctag2,RNA_tags{kc});
%         if ~isempty(id)
%             if ~isempty(pc)
%                 add_msg_board(sprintf('ERROR: RNA linker %i belongs to two contiguous chains by linking to chains %s and %s. Aborting.\n',kl,ctag1,ctag2));
%                 fclose(fid);
%                 return
%             else
%                 pc = kc;
%             end
%         end
%     end
%     if isempty(pc) % this is a new contiguous chain
%         RNA_chains = RNA_chains + 1;
%         nctags = sprintf(':%i:',klr);
%         if ~isempty(ctag1)
%             nctags = sprintf(':%s%s',ctag1,nctags);
%         end
%         if ~isempty(ctag2)
%             nctags = sprintf('%s%s:',nctags,ctag2);
%         end
%         RNA_tags{RNA_chains} = nctags;
%     else % the linker must be inserted at the right place
%         ctags = RNA_tags{pc};
%         p5poi = strfind(ctags,ctag1);
%         p3poi = strfind(ctags,ctag2);
%         if ~isempty(p5poi) && ~isempty(p3poi)
%             nctags = sprintf('%s:%i:%s',ctags(1:p5poi),klr,ctags(p3poi:end));
%         elseif ~isempty(p5poi)
%             nctags = sprintf('%s:%i:',ctags(1:p5poi),klr);
%             if ~isempty(ctag2)
%                 nctags = sprintf('%s%s:',nctags,ctag2);
%             end
%         else
%             nctags = sprintf(':%i:%s',klr,ctags(p3poi:end));
%             if ~isempty(ctag1)
%                 nctags = sprintf(':%s%s',ctag1,nctags);
%             end
%         end
%         RNA_tags{pc} = nctags;
%     end
% end
if isfield(restraints,'RNA')
    if ~isfield(restraints.RNA,'bind')
        restraints.RNA.bind = [];
    end
    restraints.RNA_tags = cell(1,length(restraints.RNA.bind));
    for krna = 1:length(restraints.RNA.bind)
        poia = strfind(restraints.RNA.bind(krna).anchora,'(');
        poie = strfind(restraints.RNA.bind(krna).anchora,')');
        restraints.RNA_tags{krna} = restraints.RNA.bind(krna).anchora(poia:poie);
    end
end
if ~isfield(restraints,'models') && isfield(restraints,'ensemble')
    restraints.models = restraints.ensemble;
end
failed = false;

function label_std = check_label(label)
% checks if a rotamer library exists for a requested spin label and returns
% the MMM-internal three-letter code for this label
% an empty string is returned, if the label does not exist

global rotamer_libraries

label_std = '';

for k = 1:length(rotamer_libraries)
    if strcmpi(label,rotamer_libraries(k).label) || strcmpi(label,rotamer_libraries(k).tc)
        label_std = rotamer_libraries(k).tc;
    end;
end

function [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(adr,label,label_adr,label_types,label_coor,label_indices,stemlibs,rbdef)
% returns chain number, label coordinate, and label position uncertainty as
% well as labeled residue indices
% (ecoor empty, if not defined) 

global hMain
global model

ml = sum(label_coor(:,1) ~= 0);

labnum = [];
labnums = tag2ids(adr,label_adr);
if ~isempty(labnum)
    for k = 1:length(labnums)
        type = id2tag(labnums(k),label_types);
        if strcmpi(label,type)
            labnum = labnums(k);
        end
    end
end

if ~isempty(labnum)
    ecoor = label_coor(labnum,:);
    indices = label_indices(labnum,:);
else
    [indices,message] = resolve_address(adr);
    if message.error % this could be a label site in a stemloop library
        indices = [];
        ecoor = [];
        for lib = 1:length(stemlibs)
            library = stemlibs{lib};
            for lab = 1:length(library.labelsites)
                if strcmp(library.labelsites(lab).adr,adr)
                    [ml,~] = size(library.labelsites(lab).coor);
                    ecoor = library.rba*ones(ml,5);
                    ecoor(:,2:4) = library.labelsites(lab).coor;
                    ecoor(:,5) = library.labelsites(lab).rmsd';
                end
            end
        end
        return
    else
        for rb = 1:length(rbdef)
            chains = rbdef(rb).chains;
            if min(abs(chains-indices(2))) == 0
                rba = rb;
            end
        end
        ecoor(1) = rba;
        command=sprintf('rotamers %s %s %i',adr,label,298);
        hMain.store_undo=false;
        hMain.dynamic_rotamers=false;
        cmd(hMain,command);
        labels = label_information(model.sites);
        found = false;
        for k = 1:length(labels)
            if abs(sum(indices - labels(k).indices)) == 0
                found = true;
                ecoor(2:4) = labels(k).xyz;
                ecoor(5) = labels(k).rmsd;
            end
        end
        if found
            label_coor(ml+1,:) = ecoor;
            label_indices(ml+1,:) = indices;
            label_adr = [label_adr adr ':'];
            label_types = [label_types label ':'];
        else
            ecoor = [];
        end
    end;
end
    

function labels = label_information(sites)

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

function [krb,kp] = identify(indices,label,restraints)

krb = 0;
kp = 0;
for kr = 1:length(restraints.pindices)
    cindices = restraints.pindices{kr};
    if min(abs(indices(2)-restraints.rb(kr).chains)) == 0 % chain is in this rigid body
        krb = kr;
        [mc,~] = size(cindices);
        for kpc = 1:mc
            if abs(sum(cindices(kpc,:)-indices)) == 0 % same residue indices
                clabels = restraints.plabels{kr};
                if strcmpi(label,clabels{kpc}) % same spin label
                    kp = kpc;
                end
            end
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

function  ctag = get_ctag(adr)

poi1 = strfind(adr,'(');
poi2 = strfind(adr,')');
if ~isempty(poi1) && ~isempty(poi2) && poi2 > poi1
    ctag = adr(poi1+1:poi2-1);
else
    ctag = '';
end


function [coor,ind,sllib] = get_RNA_anchor(adr,stemlibs)

coor = [];
ind = [];
sllib = [];

for lib = 1:length(stemlibs)
    library = stemlibs{lib};
    for link = 1:length(library.linksites)
        if strcmpi(adr,library.linksites(link).adr)
            coor = library.linksites(link).coor;
            ind = -library.rba;
            sllib = lib;
        end
    end
end
