function restraints=rd_restraints(fname)
% function restraints=rd_restraints(fname)
%
% reads mixed experimental restraints from an ASCI file
% restraints are returned in a structure whose fields correspond to
% different restraint types
%
%
% currently implemented types are:
%
% DEER      distance between two labelled sites
% DIRECT    distance between C_alpha atoms

default_helix_mindist=0.65; % default value for the closest approach of two helix axes

restraints=[];

fid=fopen(fname);
if fid==-1
    add_msg_board('ERROR: Constraint file does not exist');
    return;
end

clear restraints
restraints.ensemble=1;
restraints.uncertainty=0;
restraints.exclude=true;
restraints.realign=false;
restraints.target='';
restraints.helix_mindist=default_helix_mindist;
restraints.ref_chain='';
restraints.SAXS = [];
restraints.max_time = 1; % default


DEER_poi=0;
direct_poi=0;
atom_poi = 0;
align_poi=0;
helix_poi=0;
strand_poi=0;
sheet_poi=0;
disp_poi=0;
locate_poi=0;
network_poi=0;
ref_poi=0;
oligo_poi=0;
depth_poi=0;
aprop_poi = 0;
bprop_poi = 0;
pprop_poi = 0;
SAXS_poi = 0;
mode=0;
scale_units=1;

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
        % fprintf(1,'%s\n',tline);
        args=myline{1};
        if strcmp(char(args(1)),'#')
            switch upper(char(args(2)))
                case 'FLEX'
                    mode = 0;
                case 'GENERIC'
                    mode = 0;
                case 'PATH'
                    mode = 0;
                    datapath = char(args(3));
                    genpath(datapath);
                case 'TRANSITION'
                    mode=0;
                    restraints.initial=char(args(3));
                    restraints.final=char(args(4));
                case 'PDB'
                    mode=0;
                    pdbid=strtok(char(args(3)),':');
                    restraints.PDB=pdbid; % to conform to old format
                    for k=3:length(args)
                        [pdbid,chains]=strtok(char(args(k)),':');
                        restraints.template(k-2).pdbid=pdbid;
                        if length(chains)>1
                            nonsense=textscan(chains(2:end),'%s','Delimiter',',');
                            chains=nonsense{1};
                            chaintags=':';
                            for kk=1:length(chains)
                                chaintags=[chaintags char(chains(kk)) ':'];
                            end
                            restraints.template(k-2).chaintags=chaintags;
                        else
                            restraints.template(k-2).chaintags='';
                        end
                    end
                case {'ALIGN','ALIGNMENT'}
                    mode=0;
                    restraints.alignment=char(args(3));
                case 'BUNDLE'
                    mode=0;
                    bundle=zeros(1,length(args)-2);
                    for k=3:length(args)
                        bundle(k-2)=str2double(char(args(k)));
                    end
                    restraints.bundle=bundle;
                case 'CHAINS'
                    mode=0;
                    chains=':';
                    for k=3:length(args)
                        nonsense=textscan(char(args(k)),'%s','Delimiter',',');
                        chainlist=nonsense{1};
                        for kk=1:length(chainlist)
                            chains=[chains char(chainlist(kk)) ':'];
                        end
                    end
                    restraints.chains=chains;
                case 'SAXS'
                    SAXS_poi = SAXS_poi + 1;
                    restraints.SAXS(SAXS_poi).data = char(args(3));
                    if length(args) > 3
                        restraints.SAXS(SAXS_poi).sm = str2double(char(args(4)));
                    else
                        restraints.SAXS(SAXS_poi).sm = [];
                    end

                case 'ENSEMBLE'
                    mode=0;
                    restraints.ensemble=str2double(char(args(3)));
                    if length(args)>3
                        restraints.uncertainty=scale_units*str2double(char(args(4)));
                        restraints.rescale = scale_units;
                    else
                        restraints.uncertainty=0;
                    end
                    if length(args)>4
                        if strcmpi(char(args(5)),'all')
                            restraints.exclude=false;
                        else
                            restraints.exclude=true;
                        end
                    else
                        restraints.exclude=true;
                    end
                case 'CANCHOR'
                    mode=0;
                    restraints.Canchor=char(args(3));
                case 'NANCHOR'
                    mode=0;
                    restraints.Nanchor=char(args(3));
                case 'DEER'
                    mode=1;
                    label=char(args(3));
                    if length(args) > 3
                        T=char(args(4));
                    else
                        T= '298';
                    end
                case 'ODEER'
                    mode=16;
                    label1 = char(args(3));
                    [label1,type1] = check_label(label1);
                    if isempty(label1)
                        add_msg_board(sprintf('Warning: Label(1) %s is unknown in line %i. Reverting to MTSL.',char(args(3))),nl);
                        label1 = 'R1A';
                        type1 = 1;
                    end
                    label2 = char(args(4));
                    [label2,type2] = check_label(label2);
                    if isempty(label2)
                        add_msg_board(sprintf('Warning: Label(2) %s is unknown in line %i. Reverting to MTSL.',char(args(4))),nl);
                        label2 = 'R1A';
                        type2 = 1;
                    end
                case 'OLIGOMER'
                    mode=11;
                    multiplicity = str2double(char(args(3)));
                    label=char(args(4));
                    [label_std,label_type] = check_label(label);
                    if length(args)>4
                        T=char(args(5));
                    else
                        T = '298';
                    end
                case 'DEPTH'
                    mode=12;
                    label=char(args(3));
                    if length(args)>4
                        T=char(args(5));
                    else
                        T = '298';
                    end
                case 'CISPEPTIDE'
                    cpp  = 0;
                    cispeptides = zeros(1,1000);
                    for k=3:length(args)
                        cpp = cpp + 1;
                        cispeptides(cpp)=str2double(char(args(k)));
                    end
                    restraints.cispeptides = cispeptides(1:cpp);
                case 'OUTPUT'
                    restraints.output='';
                    if length(args)>2
                        restraints.output=char(args(3));
                    end
                case 'DIRECT'
                    mode=2;
                case 'DRAG' % non-peptide residues to be dragged along with an elastic network
                    mode=17;
                case 'LOCATE'
                    mode=8;
                    if length(args)>2
                        tag=char(args(3));
                    else
                        tag='loc';
                    end
                    if length(args)>3
                        probability=str2double(char(args(4)));
                    else
                        probability=0.5;
                    end
                    if length(args)>4
                        label=char(args(5));
                    else
                        label='MTSL';
                    end
                    if length(args) > 5
                        T=char(args(6));
                    else
                        T='298';
                    end
                    if length(args)>6
                        display_mode=char(args(7));
                    else
                        display_mode='none';
                    end
               case 'NETWORK'
                    mode=9;
                    if length(args)>2
                        pid=char(args(3));
                    else
                        pid='DGM1';
                    end
                    if length(args)>3
                        probability=str2double(char(args(4)));
                    else
                        probability=0.5;
                    end
                    if length(args)>4
                        label=char(args(5));
                    else
                        label='MTSL';
                    end
                    if length(args) > 5
                        T=char(args(6));
                    else
                        T='298';
                    end
                    if length(args)>6
                        display_mode=char(args(7));
                    else
                        display_mode='constraints';
                    end
                case 'REALIGN'
                    mode=3;
                    if length(args)>2
                        restraints.align_threshold=char(args(3));
                    else
                        restraints.align_threshold=0;
                    end
                    restraints.realign=true;
                case 'REFERENCE'
                    mode=10;
                case 'HELICES'
                    mode=4;
                    if length(args)>2 % minimum helix axes distance is provided
                        restraints.helix_mindist=scale_units*str2double(char(args(3)));
                    end
                case 'SEED'
                    mode=0;
                    restraints.seed = round(str2double(char(args(3))));
                case 'STRANDS'
                    mode=5;
                case 'SHEETS'
                    mode=6;
                case 'APROP'
                    mode=13;
                case 'BPROP'
                    mode=14;
                case 'PPROP'
                    mode=15;
                case 'DISPLACEMENTS'
                    mode=7;
                    weights=char(args(3));
                case 'TARGET'
                    mode=0;
                    pdbid=strtok(char(args(3)),':');
                    restraints.target=pdbid; 
                % paramaters specific to elastic network model fitting
                case 'BASIS' % initial active space dimension
                    mode=0;
                    restraints.basis=str2double(char(args(3)));
                case 'TIF' % threshold for active space extension
                    mode=0;
                    restraints.tif=str2double(char(args(3)));
                case 'MMAX' % iteration number at which active space
                    mode=0; % achieves maximum dimension
                    restraints.mmax=str2double(char(args(3)));
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
                case 'LOWER'
                    mode=0;
                    restraints.lower = scale_units*str2double(char(args(3)));
                case 'UPPER'
                    mode=0;
                    restraints.upper = scale_units*str2double(char(args(3)));
                case 'SYMMETRY'
                    mode=0;
                    restraints.Cn=str2double(char(args(3)));
                    restraints.ref_chain=char(args(4));
                    for k=5:length(args)
                        restraints.chain_copies{k-4}=char(args(k));
                    end
                case 'PARALLEL'
                    mode=0;
                    for k=3:length(args)
                        restraints.parallel(k-2)=str2double(char(args(k)));
                    end
                case 'PERPENDICULAR'
                    mode=0;
                    for k=3:length(args)
                        restraints.perpendicular(k-2)=str2double(char(args(k)));
                    end
                case 'EXTEND'
                    mode=0;
                    restraints.initial_model=char(args(3)); 
                case 'MAXTIME'
                    mode=0;
                    restraints.max_time=str2double(char(args(3))); 
                case 'TOLERANCE'
                    mode=0;
                    restraints.tolerance=str2double(char(args(3))); 
                case 'DOMAIN'
                    mode=0;
                    restraints.domain_start = char(args(3)); 
                    restraints.domain_end = char(args(4)); 
                case 'SEQUENCE'
                    mode=0;
                    restraints.sequence = char(args(3)); 
                case 'END'
                    mode=-1;
                otherwise
                    mode=0;
                    add_msg_board('Warning: Unknown restraint mode');
            end
        elseif mode>0 && ~strncmp(strtrim(char(args(1))),'%',1)
            switch mode
                case 1
                    DEER_poi=DEER_poi+1;
                    restraints.DEER(DEER_poi).label=label;
                    restraints.DEER(DEER_poi).T=str2double(T);
                    restraints.DEER(DEER_poi).adr1=char(args(1));
                    restraints.DEER(DEER_poi).adr2=char(args(2));
                    restraints.DEER(DEER_poi).r=scale_units*str2double(char(args(3)));
                    restraints.DEER(DEER_poi).sigr=scale_units*str2double(char(args(4)));
                    restraints.DEER(DEER_poi).type= 1;
                    if length(args) > 4
                        restraints.DEER(DEER_poi).file = char(args(5));
                    end
               case 16
                    DEER_poi = DEER_poi + 1;
                    restraints.DEER(DEER_poi).label = [label1 '|' label2];
                    restraints.DEER(DEER_poi).T = 298;
                    restraints.DEER(DEER_poi).type = [type1 type2];
                    argpoi = 0;
                    if type1 == 1 || type1 == 3
                        restraints.DEER(DEER_poi).adr1=char(args(1));
                        argpoi = argpoi + 1;
                    elseif type1 == 2
                        restraints.DEER(DEER_poi).adr1 = [char(args(1)) '|' char(args(2))];
                        argpoi = argpoi + 2;
                    end
                    if type2 == 1 || type2 == 3
                        restraints.DEER(DEER_poi).adr2=char(args(argpoi+1));
                        argpoi = argpoi + 1;
                    elseif type1 == 2
                        restraints.DEER(DEER_poi).adr2 = [char(args(argpoi+1)) '|' char(args(argpoi+2))];
                        argpoi = argpoi + 2;
                    end
                    restraints.DEER(DEER_poi).r = scale_units*str2double(char(args(argpoi+1)));
                    restraints.DEER(DEER_poi).sigr = scale_units*str2double(char(args(argpoi+2)));
                    if length(args) > argpoi+2
                        restraints.DEER(DEER_poi).file = char(argpoi+3);
                    end
                case 2
                    direct_poi=direct_poi+1;
                    restraints.direct(direct_poi).adr1=char(args(1));
                    restraints.direct(direct_poi).adr2=char(args(2));
                    restraints.direct(direct_poi).r=scale_units*str2double(char(args(3)));
                    restraints.direct(direct_poi).sigr=scale_units*str2double(char(args(4)));
                case 3
                    align_poi=align_poi+1;
                    restraints.aligned(align_poi).adr1=char(args(1));
                    restraints.aligned(align_poi).adr2=char(args(2));
                case 4
                    helix_poi=helix_poi+1;
                    restraints.helices(helix_poi).adr=char(args(1));
                case 5
                    strand_poi=strand_poi+1;
                    restraints.strands(strand_poi).adr=char(args(1));
                case 6
                    sheet_poi=sheet_poi+1;
                    restraints.sheets(sheet_poi).adr1=char(args(1));
                    restraints.sheets(sheet_poi).adr2=char(args(2));
                    restraints.sheets(sheet_poi).length=char(args(3));
                case 7
                    disp_poi=disp_poi+1;                    
                    restraints.displace(disp_poi).adr=char(args(1));
                    cvec=zeros(1,3);
                    cvec(1)=str2double(char(args(2)));
                    cvec(2)=str2double(char(args(3)));
                    cvec(3)=str2double(char(args(4)));
                    restraints.displace(disp_poi).vec=cvec;
                    if ~strcmpi(weights,'uniform') && length(args)>=5
                        restraints.displace(disp_poi).weight=str2double(char(args(5)));
                    else
                        restraints.displace(disp_poi).weight=1;
                    end
                case 8 % localization
                    locate_poi=locate_poi+1;
                    restraints.locate(locate_poi).display_mode=display_mode;
                    restraints.locate(locate_poi).tag=tag;
                    restraints.locate(locate_poi).probability=probability;
                    restraints.locate(locate_poi).adr=char(args(1));
                    restraints.locate(locate_poi).r=scale_units*str2double(char(args(2)));
                    restraints.locate(locate_poi).sigr=scale_units*str2double(char(args(3)));
                    restraints.locate(locate_poi).label = label;
                    restraints.locate(locate_poi).T=str2double(T);
                    if length(args)>3
                        restraints.locate(locate_poi).file=char(args(4));
                        restraints.locate(locate_poi).fulldistr=true;
                    else
                        restraints.locate(locate_poi).fulldistr=false;
                    end
                case 9 % network
                    network_poi=network_poi+1;
                    restraints.network(network_poi).pid=pid;
                    restraints.network(network_poi).label=label;
                    restraints.network(network_poi).T=str2double(T);
                    restraints.network(network_poi).probability=probability;
                    restraints.network(network_poi).adr1=char(args(1));
                    restraints.network(network_poi).adr2=char(args(2));
                    restraints.network(network_poi).r=scale_units*str2double(char(args(3)));
                    restraints.network(network_poi).sigr=scale_units*str2double(char(args(4)));
                    restraints.network(network_poi).display_mode=display_mode;
                    if length(args)>4
                        restraints.network(network_poi).file=char(args(5));
                        restraints.network(network_poi).fulldistr=true;
                    else
                        restraints.network(network_poi).fulldistr=false;
                    end
                case 10 % reference
                    ref_poi=ref_poi+1;
                    restraints.reference(ref_poi).tag=char(args(1));
                    x = str2double(char(args(2)));
                    y = str2double(char(args(3)));
                    z = str2double(char(args(4)));
                    restraints.reference(ref_poi).xyz = [x,y,z];
                    if length(args) > 4
                        restraints.reference(ref_poi).rmsd = str2double(char(args(5)));
                    else
                        restraints.reference(ref_poi).rmsd = 0;
                    end
                case 11 % oligomer
                    oligo_poi=oligo_poi+1;
                    restraints.oligomer(oligo_poi).label=label_std;
                    restraints.oligomer(oligo_poi).T=str2double(T);
                    restraints.oligomer(oligo_poi).multi=multiplicity;
                    restraints.oligomer(oligo_poi).adr=char(args(1));
                    switch label_type
                        case 1
                            restraints.oligomer(oligo_poi).adr=char(args(1));
                            restraints.oligomer(oligo_poi).r=scale_units*str2double(char(args(2)));
                            restraints.oligomer(oligo_poi).sigr=scale_units*str2double(char(args(3)));
                        case 2
                            restraints.oligomer(oligo_poi).adr=[char(args(1)) '|' char(args(2))];
                            restraints.oligomer(oligo_poi).r=scale_units*str2double(char(args(3)));
                            restraints.oligomer(oligo_poi).sigr=scale_units*str2double(char(args(4)));
                    end
                case 12 % depth
                    depth_poi=depth_poi+1;
                    restraints.depth(depth_poi).label=label;
                    restraints.depth(depth_poi).T=str2double(T);
                    restraints.depth(depth_poi).adr=char(args(1));
                    restraints.depth(depth_poi).r=scale_units*str2double(char(args(2)));
                    restraints.depth(depth_poi).sigr=scale_units*str2double(char(args(3)));
                case 13
                    aprop_poi = aprop_poi+1;
                    restraints.aprop(aprop_poi).adr=char(args(1));
                    restraints.aprop(aprop_poi).prop=str2double(char(args(2)));
                case 14
                    bprop_poi = bprop_poi+1;
                    restraints.bprop(bprop_poi).adr=char(args(1));
                    restraints.bprop(bprop_poi).prop=str2double(char(args(2)));
                case 15
                    pprop_poi = pprop_poi+1;
                    restraints.pprop(pprop_poi).adr=char(args(1));
                    restraints.pprop(pprop_poi).prop=str2double(char(args(2)));
                case 17 % residues to be dragged along with elastic network, atom addresses are stored
                    adr = char(args(1));
                    [~,indices] = get_object(adr,'children');
                    [mi,~] = size(indices);
                    for ki = 1:mi
                        [~,ctag,~,resnum,~,atag] = mk_address_parts(indices(ki,:));
                        adr = sprintf('(%s)%i.%s',ctag,resnum,atag);
                        atom_poi = atom_poi + 1;
                        restraints.atoms(atom_poi).adr = adr;
                    end
            end
        end
    end
end

fclose(fid);

function [label_std,type] = check_label(label)
% checks if a rotamer library exists for a requested spin label and returns
% the MMM-internal three-letter code for this label
% an empty string is returned, if the label does not exist
% type 1  spin label
% type 2  ligand
% type 3  atom

global rotamer_libraries
global ligand_libraries

label_std = '';
type = 0;

poi = strfind(lower(label),'atom.');
if ~isempty(poi) % this is an atom request
    type = 3;
    label_std = label(poi+5:end);
end

for k = 1:length(rotamer_libraries)
    if strcmpi(label,rotamer_libraries(k).label) || strcmpi(label,rotamer_libraries(k).tc)
        label_std = rotamer_libraries(k).tc;
        type = 1;
    end
end

if isempty(label_std)
    for k = 1:length(ligand_libraries)
        if strcmpi(label,ligand_libraries(k).label) || strcmpi(label,ligand_libraries(k).tc)
            label_std = ligand_libraries(k).tc;
            type = 2;
        end
    end
end