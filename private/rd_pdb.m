function [structure,info]=rd_pdb(fname,remarks)
% function [structure,info]=rd_pdb(fname,remarks)
%
% Reads structure information from a PDB file into a Matlab structure
% only the 
%   header, (HEADER)
%   title, (TITLE)
%   molecule name, (COMPND, token MOLECULE)
%   organism, (SOURCE, token ORGANISM_SCIENTIFIC)
%   keywords, (KEYWDS)
%   authors, (AUTHOR)
%   primary literature reference (JRNL)
%   data base reference for the sequence, (DBREF)
%   sequence information, (SEQRES)
%   sequence modifications, (MODRES)
%   information on nonstandard residues without their names (HET)
%   secondary structure information (HELIX, SHEET),
%   sites information (SITE)
%   unit cell information for crystal structures (CRYST1, ORIGX1, ORIGX2,
%   ORIGX3, SCALE1, SCALE2, SCALE3) 
%   atom coordinates, (ATOM, HETATM)
%   connection information (CONECT)
%   model numbers (MODEL, ENDMDL)
% are read
% in particular, REMARK information is not fully read, as it is not 
% (always) well formatted, the program checks only for the presence of 
% remarks and provides a list of remark numbers, 
% thus, the client at least knows whether there is additional information,
% unprocessed remarks are stored in a separate file, if rd_pdb is called
% with two arguments
% The following remarks are processed:
%      1    related references
%      2    resolution
%      6    Modeller free text annotation
%    470    missing atoms
%    620    metal coordination
%    800    site declarations
%
% PDB file format v 3.2 is implemented
% the non-conforming MolProbity output files (END line too short, wrong
% tags for hydrogen atoms) are supported
% PDB-like formats with at least 54 columns in ATOM and HETATM records are
% supported, occupancy is then set to 1.0, B-factor to 0.0, element symbols
% are derived from columns 13:14, which usually works, but is not
% guaranteed to work, a warning flag 1 is set
%
% several sets of atom coordinates (PDB jargon: "models") per file and 
% several chains per structure are supported,
% alternate location is supported
% anisotropic temperature factors are not supported,
% charges on atoms are not supported, as they are NOT experimentally known
% 
% fname         input file name, if without extension, it is appended by 
%               '.pdb'
% remarks       optional output file name for unprocessed remarks, this
%               allows to post-process these remarks for annotations after
%               the whole structure is defined
%
% structure     output structure, array of structures for separate chains
%               structure(i).name           chain identifier
%               structure(i).resnum         number of residues
%               structure(i).sequence       sequence (single-letter code)
%               structure(i).seqtype        sequence type (0: unknown, 1: amino acid, 2: nucleic acid)
%               structure(i).header         sequence header (if available)
%               structure(i).restags        sequence (three-letter code and
%                                           MMM tag format), full sequence
%               structure(i).modified       total number of modified residues
%               structure(i).mutations(k)   sequence modifications, array of
%                                           structures .number, .original, 
%                                           .modified
%               structure(i).nonstandard    total number of nonstandard residues
%               structure(i).het(k)         array of nonstandard residues with
%                                           .number number in sequence
%                                           .id     Het identifier without
%                                                   leading blanks
%               structure(i).helices        number of helices
%               structure(i).helix(k)       array of helix records with
%                                           .start first residue
%                                           .end   last residue
%                                           .name  helix name
%                                           .class helix class
%               structure(i).strands        number of strands in sheets
%               structure(i).strand(k)      array of strand records with
%                                           .start first residue
%                                           .end   last residue
%                                           .s_c   chain of first residue
%                                           .e_c   chain of last residue
%                                           .name  sheet name
%                                           the first residue decides under
%                                           which chain the strand is
%                                           stored
%               structure(i).conn           connection table for chain i,
%                                           array(n,maxconn)
%               structure(i).maxconn        maximum number of connections 
%                                           (bonds) for an atom of chain i
%               structure(i).dbref          data bank reference for chain i
%               structure(i).isotopes       atomic number and mass of all
%                                           atoms of chain i, array(n,2)
%                                           with the first index matching
%                                           the one in .xyz
%               structure(i).xyz{j}         atom coordinates for the j-th
%                                           "model" of chain i, array(n,3)
%               structure(i).atoms{j}       number of atoms in the j-th
%                                           model of chain i
%               structure(i).Bfactor{j}     B-factors of atoms in the j-th
%                                           model of chain i
%                                           8*pi^2 coordinate variance
%               structure(i).Btensor{j}     anisotropic temperature factors 
%                                           of atoms in the j-th
%                                           model of chain i, array[1,6]
%                                           10^4 coordinate variance
%               structure(i).residues{j}    residue definitions with fields
%                                           .residue_tags string of defined
%                                           residue numbers, contains only
%                                           residues for which at least
%                                           partial coordinates are known
%                                           .info array(1,m), where m is
%                                           the number of defined residues
%                                           with the following fields:
%                                           .name residue name
%                                           .hetflag flag, hetflag==1
%                                           indicates non-standard residues
%                                           read from HETATM records,
%                                           hetflag==0 standard residues
%                                           read from ATOM records
%                                           .secondary secondary structure
%                                           identifier 0: none, 1: helix,
%                                           2: strand
%                                           .connected flag that tells
%                                           whether internal connections
%                                           for a standard residue were made
%                                           .number residue number
%                                           .type  1 amino acid
%                                                  2 nucleic acid 
%                                           .atom_tags string with atom tags
%                                           .atom_numbers{l} cell array of
%                                           atom numbers in the same
%                                           sequence as the atom tags, for
%                                           models without alternate atom
%                                           positions, each cell contains a
%                                           single number, otherwise it
%                                           contains an array(m,3) of the m
%                                           alternate locations with 1st
%                                           column atom number, 2nd column
%                                           occupancy, 3rd column location
%                                           id
%                                           .elements array with element
%                                           numbers, required and dependable
%                                           only if hetflag==1, if no
%                                           element information is present
%                                           in an atom record, it is given
%                                           as '?'
%                                           .location_tags, string with
%                                           used location tags, ':A:' if no
%                                           alternate locations
%                                           .terminal exists only, if this
%                                           is the terminal residue of a
%                                           chain (TER record), then it has
%                                           the value 1
%                                           .absresnum absolute residue
%                                           number (for PRONOX)
% info          Matlab structure with the following fields (if information
%               is present in file)
%               .idCode        PDB identifier
%               .class         structure class (from PDB header)
%               .depDate       deposition date
%               .title         title line(s) of the PDB structure
%               .remarks       list of remarks (only numbers) present in the PDB file, see
%                               PDB file format definition for details
%               .chain_tags    tags of the chains defined in the structure (string)
%               .chain_ids     list of numeric chain identifiers
%               .warning       warning flag, if element information is missing
%               .center        center coordinates (mean coordinate of all
%                              atoms)
%               .atoms         total number of atoms
%               .residues      total number of residues
%               .B_range       B factor range
%               .missing       missing atoms (cellstring of partial addresses)
%               .site_tags     tag list for sites (SITE records)
%               .sites         site information
%               .SSbonds       list of sulfur-sulfur bonds (SSBOND record)
%               .authors       authors (AUTHOR record)
%               .keywords      keywords (KEYWDS record)
%               .references    references (JRNL record and REMARK 1 records)
%               .molecule      molecule name (from, COMPND/MOLECULE: record)
%               .organism      organism name (from, SOUREC/ORGANISM_SCIENTIFIC: record)
%               .Modeller_obj  Modeller objective function
%               .Modeller_sid  Modeller sequence identity information
%
% Known bugs:
% - if there are several models with different numbers of alternate
%   positions for the same atom, the connection tables will be wrong; this
%   is inherent in the PDB definition v 3.2 that includes only ATOM and HETATM
%   records in the model definitions, but not CONECT records. Hence, in PDB
%   several models (actually atom coordinate sets) do not allow for
%   alternate positions within individual models
%
% G. Jeschke, 2009-2011


add_msg_board('Loading and interpreting PDB file. Please be patient.');
drawnow;

maxat = scan_pdb_fast(fname);

min_resnum=999999;
max_resnum=0;

seqtype=0;

    if exist('onCleanup','class')==8 % this is done only for Matlab version 
                                      % 2008a and later, where onCleanup exists
        c = onCleanup(@myCleanup);
    end
    
    occupancy_error=0.001;  % if occupancy deviates by less than this from 1.00,
                            % no alternate locations are recognized (is implied
                            % by PDB format v 3.2)
    max_atoms = maxat; % for preallocation of coordinate arrays, corresponds to single chain

    atom_corr=zeros(5*max_atoms,3); % atom number correspondence table for proper storage of connections

    global residue_defs
    global chemistry
    global web_adr
    global queries
    global general

    if isempty(strfind(fname,'.')), fname=strcat(fname,'.pdb'); end

    if nargin>1
        rem_out=fopen(remarks,'wt');
    else
        rem_out=[];
    end

    sr=0;
    at=0;
    het=0;

    structure=[]; % initialize all output arguments to avoid error due to missing records in PDB file
    info={};
    pdb_class='';
    idCode='';
    depDate='';
    title=''; 
    db_tags=':';
    db_access={};
    db_cutout=zeros(50,2);
    chain_tags=':';
    chain_id=0;
    chain_ids=[];
    site_tags=':';
    sites=[];
    sites_declared=':';
    site_declarations={};
    site_evidence={};
    sdr_poi=1;
    seq='';
    restags=':';
    B_range=[1e6,0];
    SSbonds=[];

    remarks=zeros(1,1000);
    Btensor=zeros(1,6);
    rem_num=0;
    curr_model=1; % default model number
    structures=1;
    warning=0;
    atoms=0;
    in_remark=0;
    absresnum=0;
    seq_found=zeros(1,50); % flags that tell whether sequence information is present for the chains
    total_residues=0;
    center=zeros(1,3);
    missing={};
    misspoi=0;
    metal=[];
    metalpoi=0;
    authors='';
    keywords='';
    references=[];
    ins_poi=0; % index for inserted residues
    insertions=zeros(50,3); % small preallocation, but will be sufficient in almost all cases
    alternate=false;
    resolution=[];
    element_warn=false;
    rotamers=false;
    molecule='undefined';
    organism='undefined';
    Modeller_obj=[];
    Modeller_sid=[];
    
    chain_type=0;
    terminated_chains=':';
    
    cyana_flag=false;

    fid=fopen(fname);
    if fid==-1
        info.no_file=1;
        return;
    end

    nl=0;
    while 1
        tline = fgetl(fid);
        nl=nl+1;
        if ~ischar(tline), break, end
%         fprintf(1,'%s\n',tline); % echo for debugging
        if length(tline)>=6, % catches too short end line of MolProbity files
            record=tline(1:6);
        else
            record='END   ';
        end;
        switch record
            case 'HEADER'
                if length(tline)>=50,
                    pdb_class=strtrim(tline(11:50));
                elseif length(tline)>=11,
                    pdb_class=strtrim(tline(11:end));
                else
                    pdb_class='';
                end;
                if length(tline)>=59,
                    depDate=tline(51:59);
                else
                    depDate='';
                end;
                if length(tline)>=66,
                    idCode=tline(63:66);
                    if ~strcmpi(strtrim(idCode),idCode), idCode = ''; end;
                else
                    idCode='';
                end;
            case 'TITLE '
                cont=sscanf(tline(9:10),'%i');
                if isnan(cont),
                    title=tline(11:end);
                else
                    title=strvcat(title,tline(11:end));
                end;
            case 'COMPND'
                if length(tline)>=12,
                    [token,rem]=strtok(tline(12:end),':');
                    if strcmpi(strtrim(token),'MOLECULE') && length(rem)>1,
                        molecule=strtrim(rem(2:end));
                        if strcmp(molecule(end),';')
                            molecule=molecule(1:end-1);
                        end;
                        molecule=strtrim(molecule);
                    end;
                end;
            case 'SOURCE'
                if length(tline)>=12,
                    [token,rem]=strtok(tline(12:end),':');
                    if strcmpi(strtrim(token),'ORGANISM_SCIENTIFIC') && length(rem)>1,
                        organism=strtrim(rem(2:end));
                        if strcmp(organism(end),';')
                            organism=organism(1:end-1);
                        end;
                    end;
                end;
            case 'KEYWDS'
                cont=sscanf(tline(9:10),'%i');
                if isnan(cont),
                    keywords=strtrim(tline(11:end));
                else
                    keywords=strtrim([keywords,' ',tline(11:end)]);
                end;
            case 'AUTHOR'
                cont=sscanf(tline(9:10),'%i');
                if isnan(cont),
                    authors=strtrim(tline(11:end));
                else
                    authors=strtrim([authors,' ',tline(11:end)]);
                end;
            case 'JRNL  '
                if isempty(references),
                    references(1).format=6;
                    references(1).type=1;
                    references(1).title='';
                    references(1).book_title='';
                    references(1).authors='';
                    references(1).editors='';
                    references(1).journal='';
                    references(1).year='';
                    references(1).volume='';
                    references(1).issue='';
                    references(1).chapter='';
                    references(1).pages='';
                    references(1).DOI='';
                    references(1).publisher='';
                    references(1).city='';
                    references(1).short=sprintf('PDB_%s_primary',idCode);
                    references(1).URL='';
                    references(1).PMID=0;
                end;
                token='';
                if length(tline)>=16,
                    token=tline(13:16);
                end;
                endl=length(tline);
                if endl>70, endl=70; end;
                switch token
                    case 'AUTH'
                        cont=sscanf(tline(17:18),'%i');
                        if isnan(cont),
                            references(1).authors=strtrim(tline(20:end));
                        else
                            references(1).authors=strtrim([references(1).authors,' ',tline(20:endl)]);
                        end;            
                    case 'TITL'
                        cont=sscanf(tline(17:18),'%i');
                        if isnan(cont),
                            references(1).title=strtrim(tline(20:end));
                        else
                            references(1).title=strtrim([references(1).title,' ',tline(20:endl)]);
                        end;            
                    case 'REF '
                        cont=sscanf(tline(17:18),'%i');
                        if isnan(cont),
                            if length(tline)>=49,
                                source=tline(20:49);
                            else
                                source=tline(20:end);
                            end;
                            references(1).journal=strtrim(source);
                            if length(tline)>=55 && strcmp(tline(50:51),'V.'),
                                references(1).volume=strtrim(tline(52:55));
                            end;
                            if length(tline)>=61,
                                references(1).pages=strcat(strtrim(tline(57:61)),'-?');
                            end;
                            if length(tline)>=66,
                                references(1).year=strtrim(tline(63:66));
                            end;
                        else
                            references(1).journal=strtrim([references(1).journal,' ',tline(20:endl)]);
                        end;
                    case 'PMID'
                        id=str2num(tline(20:end));
                        if ~isnan(id),
                            references(1).PMID=id;
                        end;
                    case 'DOI '
                        references(1).DOI=strtrim(tline(20:end));
                end;
            case 'REMARK'
                if length(tline)>=10,
                    id=sscanf(tline(8:10),'%i');
                    if ~isempty(id) && ~sum(find(remarks==id)) % new remark
                        rem_num=rem_num+1;
                        remarks(rem_num)=id;
                        in_remark=id;
                        start_rd_remark=0;
                    end;
                    switch in_remark,
                        case 1, % literature references
                            if length(tline)>=20,
                                if strcmp(tline(12:20),'REFERENCE'),
                                    rec=length(references)+1;
                                    references(rec).format=6;
                                    references(rec).type=1;
                                    references(rec).title='';
                                    references(rec).book_title='';
                                    references(rec).authors='';
                                    references(rec).editors='';
                                    references(rec).journal='';
                                    references(rec).year='';
                                    references(rec).volume='';
                                    references(rec).issue='';
                                    references(rec).chapter='';
                                    references(rec).pages='';
                                    references(rec).DOI='';
                                    references(rec).publisher='';
                                    references(rec).city='';
                                    references(rec).short=sprintf('PDB_%s_[%i]',idCode,rec);
                                    references(rec).URL='';
                                    references(rec).PMID=0;
                                end;
                            end;
                            token='';
                            if length(tline)>=16,
                                token=tline(13:16);
                            end;
                            endl=length(tline);
                            if endl>70, endl=70; end;
                            switch token
                                case 'AUTH'
                                    cont=sscanf(tline(17:18),'%i');
                                    if isnan(cont),
                                        references(rec).authors=strtrim(tline(20:end));
                                    else
                                        references(rec).authors=strtrim([references(rec).authors ' ' tline(20:endl)]);
                                    end;            
                                case 'TITL'
                                    cont=sscanf(tline(17:18),'%i');
                                    if isnan(cont),
                                        references(rec).title=strtrim(tline(20:end));
                                    else
                                        references(rec).title=strtrim([references(rec).title,' ',tline(20:endl)]);
                                    end;            
                                case 'REF '
                                    cont=sscanf(tline(17:18),'%i');
                                    if isnan(cont),
                                        if length(tline)>=49,
                                            source=tline(20:49);
                                        else
                                            source=tline(20:end);
                                        end;
                                        references(rec).journal=strtrim(source);
                                        if length(tline)>=55 && strcmp(tline(50:51),'V.'),
                                            references(rec).volume=strtrim(tline(52:55));
                                        end;
                                        if length(tline)>=61,
                                            references(rec).pages=strcat(strtrim(tline(57:61)),'-?');
                                        end;
                                        if length(tline)>=66,
                                            references(rec).year=strtrim(tline(63:66));
                                        end;
                                    else
                                        references(rec).journal=strtrim([references(rec).journal,' ',tline(20:endl)]);
                                    end;
                                case 'PMID'
                                    id=sscanf(tline(20:end),'%i');
                                    if ~isnan(id),
                                        references(rec).PMID=id;
                                    end;
                                case 'DOI '
                                    references(rec).DOI=strtrim(tline(20:end));
                            end;
                        case 2, % resolution
                            if length(tline)>=41,
                                tok1=tline(12:22);
                                tok2=tline(32:41);
                                if strcmpi(tok1,'RESOLUTION.') && strcmpi(tok2,'ANGSTROMS.'),
                                    resolution=sscanf(tline(24:30),'%f');
                                    if isnan(resolution),
                                        resolution=[];
                                    end;
                                end;
                            end;
                        case 4, % format
                            poi=strfind(tline,'MMM ROTAMERS');
                            if ~isempty(poi),
                                rotamers=true;
                            end;
                        case 6, % Modeller free text annotation
                            poi=strfind(tline,'MODELLER OBJECTIVE FUNCTION:');
                            if ~isempty(poi),
                                [tok,rem]=strtok(tline,':');
                                Modeller_obj=str2double(rem(2:end));
                            end;
                            poi=strfind(tline,'% SEQ ID');
                            if ~isempty(poi),
                                [tok,rem]=strtok(tline,':');
                                Modeller_sid=str2double(rem(2:end));
                            end;
                        case 7, % format
                            poi=strfind(tline,'CYANA');
                            if ~isempty(poi),
                                cyana_flag=true;
                            end;                        
                        case 470, % missing atoms
                            if ~start_rd_remark && length(tline)>10,
                                pattern=strtrim(tline(11:end));
                                if strcmp(pattern,'M RES CSSEQI  ATOMS'),
                                    start_rd_remark=1;
                                end;
                            else
                                if length(tline)>10,
                                    myargs=textscan(tline(11:end),'%s');
                                    chain_tag=char(myargs{1}(2));
                                    res_num=char(myargs{1}(3));
                                    for ka=4:length(myargs{1}),
                                        misspoi=misspoi+1;
                                        missing{misspoi}=sprintf('(%s){1}%s.%s',chain_tag,res_num,char(myargs{1}(ka)));
                                    end;
                                end;
                            end;
                        case 620, % metal coordination
                            pattern=strtrim(strtok(tline(11:end),':'));
                            if strcmp(pattern,'COORDINATION ANGLES FOR'), % first line of record
                                tline = fgetl(fid);
                                metalpoi=metalpoi+1;
                                if length(tline)>=51,
                                    metal(metalpoi).center=['(' tline(44) ')' '{:}' strtrim(tline(45:48)) '.' strtrim(tline(51:end))];
                                end;
                                tline=fgetl(fid);
                                tline=fgetl(fid);
                                poi=0;
                                while ~strcmp(tline(12),'N') && ~isempty(strtrim(tline(11:end))) && poi<100,
                                    poi=poi+1;
                                    metal(metalpoi).coordinated{poi}=['(' tline(18) ')' '{:}' strtrim(tline(19:22)) '.' strtrim(tline(25:28))];
                                    tline=fgetl(fid);
                                end;
                                if ~isfield(metal(metalpoi),'coordinated'),
                                    metal(metalpoi).coordinated={};
                                end;
                            end;
                        case 800, % site declarations
                            [token,rem]=strtok(tline(11:end),':');
                            switch strtrim(upper(token))
                                case 'SITE_IDENTIFIER'
                                    sdr_poi=length(site_declarations)+1;
                                    if length(rem)>1,
                                        sdr_tag=strtrim(rem(2:end));
                                    else
                                        sdr_tag=sprintf('MMM_%i',sdr_poi);
                                    end;
                                    sites_declared=sprintf('%s%s:',sites_declared,sdr_tag);
                                    sdr_poi=tag2id(sdr_tag,sites_declared);
                                case 'EVIDENCE_CODE'
                                    if length(rem)>1,
                                        site_evidence{sdr_poi}=strtrim(rem(2:end));
                                    else
                                        site_evidence{sdr_poi}='unknown';
                                    end;
                                case 'SITE_DESCRIPTION'
                                    if length(rem)>1,
                                        site_declarations{sdr_poi}=strtrim(rem(2:end));
                                    else
                                        site_declarations{sdr_poi}='nothing declared';
                                    end;
                            end; 
                        otherwise
                            if ~isempty(rem_out),
                                fprintf(rem_out,'%s\n',tline);
                            end;
                    end;
                end;
            case 'DBREF '
                if length(tline)>=67,
                    tag=tline(13);
                    db_id=tag2id(tag,db_tags);                
                    if isempty(db_id),
                        db_id=length(db_access)+1;
                        db_tags=[db_tags tag ':'];
                        db_access{db_id}=[strtrim(tline(27:32)) ':' strtrim(tline(34:41))];
                        sstart=sscanf(tline(56:60),'%i');
                        send=sscanf(tline(63:67),'%i');
                    end;
                    if ~isnan(sstart) && ~isnan(send),
                        db_cutout(db_id,1)=sstart;
                        db_cutout(db_id,2)=send;
                    end;
                end;
            case 'SEQRES' % treats sequence information
                sr=sr+1;
                tag=tline(12);
%                 if double(tag)<double('A') || double(tag)>double('Z') % fix wrong chain tags
%                     tag='Z';
%                 end;
                if tag==' ', tag=get_chain_tag(terminated_chains); end;
                seq_id=tag2id(tag,chain_tags);
                if isempty(seq_id),
                    if chain_id>0,
                        structure(chain_id).sequence=seq;
                        structure(chain_id).restags=restags;
                        seq='';
                        restags=':';
                    end;
                    chain_id=chain_id+1;
                    structure(chain_id).name=tag;
                    chain_tags=[chain_tags tag ':'];
                    terminated_chains=[terminated_chains tag ':'];
                    chain_ids=[chain_ids chain_id];
                    structure(chain_id).seqtype=0;
                    structure(chain_id).modified=0;
                    structure(chain_id).nonstandard=0;
                    structure(chain_id).helices=0;
                    structure(chain_id).strands=0;
                    structure(chain_id).xyz{1}=zeros(max_atoms,3);
                    structure(chain_id).isotopes=zeros(max_atoms,2,'single');
                    structure(chain_id).Bfactor{1}=zeros(1,max_atoms);
                    structure(chain_id).Btensor{1}=zeros(max_atoms,6,'int32');
                    structure(chain_id).atoms{1}=0;
                    structure(chain_id).conn=zeros(max_atoms,30,'int32');
                    structure(chain_id).maxconn=1;
                    structure(chain_id).residues{1}.residue_tags=':';
                    structure(chain_id).residues{1}.info=[];
                end;
                seq_found(chain_id)=1;
                structure(chain_id).resnum=sscanf(tline(14:17),'%i');
                seqinfo=textscan(tline(18:end),'%s');
                num=length(seqinfo{1});
                for k=1:num,
                    tag=char(seqinfo{1}(k));
                    while length(tag)<3,
                        tag=[' ' tag];
                    end;
                    restags=[restags tag ':'];
                    id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
                    if isempty(id), 
                        id=tag2id(tag,upper(residue_defs.nucleotide_tags),residue_defs.nucleotide_slc);
                        if isempty(id),
                            id='?'; 
                        end;
                    end;
                    seq=[seq id];
                end;
            case 'MODRES' % treats chemical or posttranslational modifications
                tag=tline(17);
                c_id=tag2id(tag,chain_tags);
                if ~isempty(c_id)
                    modnum=structure(c_id).modified+1;
                    structure(c_id).mutations(modnum).number=sscanf(tline(19:22),'%i');
                    structure(c_id).mutations(modnum).original=tline(25:27);
                    structure(c_id).mutations(modnum).modified=tline(13:15);
                    if length(tline)>27,
                        structure(c_id).mutations(modnum).comment=strtrim(tline(30:end));
                    end;
                    structure(c_id).modified=modnum;
                end;
            case 'HET   '
                tag=tline(13);
                c_id=tag2id(tag,chain_tags);
                if ~isempty(c_id)
                    modnum=structure(c_id).nonstandard+1;
                    structure(c_id).het(modnum).number=sscanf(tline(14:17),'%i');
                    structure(c_id).het(modnum).id=fliplr(deblank(fliplr(tline(8:10))));
                    structure(c_id).het(modnum).hetrecs=sscanf(tline(21:25),'%i');
                    structure(c_id).nonstandard=modnum;
                    if length(tline)>30,
                        structure(c_id).het(modnum).descriptor=tline(31:end);
                    end;
                end;
            case 'HELIX '
                tag=tline(20);
                c_id=tag2id(tag,chain_tags);
                if ~isempty(c_id)
                    current=structure(c_id).helices+1;
                    structure(c_id).helix(current).start=sscanf(tline(22:25),'%i');
                    structure(c_id).helix(current).end=sscanf(tline(34:37),'%i');
                    structure(c_id).helix(current).name=strtrim(tline(12:14));
                    structure(c_id).helix(current).class=sscanf(tline(39:40),'%i');
                    structure(c_id).helices=current;
                end;
            case 'SHEET '
                tag=tline(22);
                c_id=tag2id(tag,chain_tags);
                if ~isempty(c_id)
                    current=structure(c_id).strands+1;
                    structure(c_id).strand(current).start=sscanf(tline(23:26),'%i');
                    structure(c_id).strand(current).s_c=tag2id(tline(22),chain_tags);
                    structure(c_id).strand(current).end=sscanf(tline(34:37),'%i');
                    structure(c_id).strand(current).e_c=tag2id(tline(33),chain_tags);
                    structure(c_id).strand(current).name=sprintf('%s/%s',strtrim(tline(12:14)),strtrim(tline(8:10))); 
                    structure(c_id).strands=current;
                end;
            case 'SSBOND'
                if length(tline)>=78,
                    sbo=length(SSbonds)+1;
                    SSbonds(sbo).chainID1=tline(16);
                    SSbonds(sbo).res1=sscanf(tline(18:21),'%i');
                    SSbonds(sbo).icode1=tline(22);
                    SSbonds(sbo).chainID2=tline(30);
                    SSbonds(sbo).res2=sscanf(tline(32:35),'%i');
                    SSbonds(sbo).icode2=tline(36);
                    SSbonds(sbo).symop1=tline(60:65);
                    SSbonds(sbo).symop2=tline(67:72);
                    SSbonds(sbo).distance=sscanf(tline(74:78),'%f');
                end;
            case 'SITE  '
                tag=strtrim(tline(12:14));
                site_id=tag2id(tag,site_tags);
                if isempty(site_id),
                    site_id=length(sites)+1;
                    site_tags=sprintf('%s%s:',site_tags,tag);
                    sites(site_id).tag=tag;
                    sites(site_id).residues={};
                end;
                site_str='';
                if length(tline)>19,
                    site_str=tline(19:end);
                end;
                while length(site_str)>9,
                    sc_tag=site_str(5);
                    sr_num=strtrim(site_str(6:9));
                    if ~isempty(sr_num), % avoid nonsense due to space padding
                        sc_poi=length(sites(site_id).residues)+1;
                        sites(site_id).residues{sc_poi}=sprintf('(%s){:}%s',sc_tag,sr_num); % MMM partial address
                    end;
                    if length(site_str)>11,
                        site_str=site_str(12:end);
                    else
                        site_str='';
                    end;
                end;
            case 'CRYST1'
                if length(tline)>=56,
                    cryst.a=sscanf(tline(7:15),'%f');
                    cryst.b=sscanf(tline(16:24),'%f');
                    cryst.c=sscanf(tline(25:33),'%f');
                    cryst.alpha=sscanf(tline(34:40),'%f');
                    cryst.beta=sscanf(tline(41:47),'%f');
                    cryst.gamma=sscanf(tline(48:54),'%f');
                    cryst.sGroup=tline(56:end);
                end;
                if length(tline)>=70,
                    cryst.Z=sscanf(tline(67:70),'%i');
                else
                    cryst.Z=0;
                end;
            case 'ORIGX1'
                if length(tline)>=40,
                    origx=zeros(1,4);
                    origx(1)=sscanf(tline(11:20),'%f');
                    origx(2)=sscanf(tline(21:30),'%f');
                    origx(3)=sscanf(tline(31:40),'%f');
                    if length(tline)>=55,
                        origx(4)=sscanf(tline(46:55),'%f');
                    end;
                    cryst.origx1=origx;
                end;
            case 'ORIGX2'
                if length(tline)>=40,
                    origx=zeros(1,4);
                    origx(1)=sscanf(tline(11:20),'%f');
                    origx(2)=sscanf(tline(21:30),'%f');
                    origx(3)=sscanf(tline(31:40),'%f');
                    if length(tline)>=55,
                        origx(4)=sscanf(tline(46:55),'%f');
                    end;
                    cryst.origx2=origx;
                end;
            case 'ORIGX3'
                if length(tline)>=40,
                    origx=zeros(1,4);
                    origx(1)=sscanf(tline(11:20),'%f');
                    origx(2)=sscanf(tline(21:30),'%f');
                    origx(3)=sscanf(tline(31:40),'%f');
                    if length(tline)>=55,
                        origx(4)=sscanf(tline(46:55),'%f');
                    end;
                    cryst.origx3=origx;
                end;
            case 'SCALE1'
                if length(tline)>=40,
                    origx=zeros(1,4);
                    origx(1)=sscanf(tline(11:20),'%f');
                    origx(2)=sscanf(tline(21:30),'%f');
                    origx(3)=sscanf(tline(31:40),'%f');
                    if length(tline)>=55,
                        origx(4)=sscanf(tline(46:55),'%f');
                    end;
                    cryst.scale1=origx;
                end;
            case 'SCALE2'
                if length(tline)>=40,
                    origx=zeros(1,4);
                    origx(1)=sscanf(tline(11:20),'%f');
                    origx(2)=sscanf(tline(21:30),'%f');
                    origx(3)=sscanf(tline(31:40),'%f');
                    if length(tline)>=55,
                        origx(4)=sscanf(tline(46:55),'%f');
                    end;
                    cryst.scale2=origx;
                end;
            case 'SCALE3'
                if length(tline)>=40,
                    origx=zeros(1,4);
                    origx(1)=sscanf(tline(11:20),'%f');
                    origx(2)=sscanf(tline(21:30),'%f');
                    origx(3)=sscanf(tline(31:40),'%f');
                    if length(tline)>=55,
                        origx(4)=sscanf(tline(46:55),'%f');
                    end;
                    cryst.scale3=origx;
                end;
            case 'MODEL '
                if length(tline)>=11,
                    curr_model=sscanf(tline(11:end),'%i');
                else
                    curr_model=1;
                end;                
                if curr_model>structures, 
                    structures=curr_model; 
                end;
                if curr_model>1,
                    num_chains=length(structure);
                    for k=1:num_chains,
                        structure(k).xyz{curr_model}=zeros(max_atoms,3);
                        structure(k).isotopes=zeros(max_atoms,2,'single');
                        structure(k).Bfactor{curr_model}=zeros(1,max_atoms);
                        structure(k).Btensor{curr_model}=zeros(max_atoms,6,'int32');
                        structure(k).atoms{curr_model}=0;
                        structure(k).residues{curr_model}.residue_tags=':';
                        structure(k).residues{curr_model}.info=[];
                    end;
                    cc_id = 1;
                    chain_type = 0;
                    chain_tag = id2tag(1,chain_tags);
                end;
                model_open=1;
            case 'ENDMDL'
                model_open=0;
                curr_model = [];
            case {'ATOM  ','HETATM'}
                curr_res=sscanf(tline(23:26),'%i');
                if curr_res<min_resnum, min_resnum=curr_res; end;
                if curr_res>max_resnum, max_resnum=curr_res; end;
                % residue information, occupancy is also stored, when needed
                if ~isempty(strtrim(tline(27))), % care about insertion code crap
                    inserted=1;
                else
                    inserted=0;
                end;
                atom_tag=tline(13:16);
                if strcmpi(record,'ATOM  ') && atom_tag(4) == '*',
                    atom_tag(4) = '''';
                end;
                if strcmp(atom_tag(1),'Q') || strcmp(atom_tag(2),'Q'),
                    continue;
                end;
                curr_resname=tline(18:20);
                if strcmpi(curr_resname,'DUM') || strcmpi(curr_resname,'PL ') || strcmpi(curr_resname,'LP '),
                    continue;
                end;
                x=sscanf(tline(31:38),'%f');
                y=sscanf(tline(39:46),'%f');
                z=sscanf(tline(47:54),'%f');
                if x==9999 && y==9999 && z==9999,
                    continue;
                end;
                at=at+1;
                hetflag=strcmp(record,'HETATM');
                if hetflag, het=het+1; end;
                if isstrprop(atom_tag(1),'digit') % wrong hydrogen tag format of MolProbity
                    atom_tag=[deblank(atom_tag(2:4)) atom_tag(1)]; % fix it to proper format
                end;
                atom_tag=strtrim(atom_tag);
                element='';
                if length(tline)>=78, % get element symbol, if present
                    element=strtrim(tline(77:78));
                end;
                if isempty(element) || strcmpi(element,'X'), % the latter is needed because of the wrong Rosie FARFAR format 
                    warning=1;
                    element=tline(13:14);
                    if cyana_flag,
                        element=[' ' tline(13)];
                    end;
                    if isstrprop(element(1),'digit'),
                        element(1)=' ';
                    end;
                    for k=1:2,
                        if isstrprop(element(k),'digit'), element(k)=' '; end;
                    end;
                    if element(1)=='H' && length(atom_tag)==4, % fix Biomer format, actually this is proper format for some H
                        element=' H';
                    end;
                    if element(1)=='C' && length(atom_tag)==4, % fix Biomer format
                        element=' C';
                        element_warn=1;
                    end;
                    if element(1)=='N' && length(atom_tag)==4, % fix Biomer format
                        element=' N';
                        element_warn=1;
                    end;
                    if element(1)=='O' && length(atom_tag)==4, % fix Biomer format
                        element=' O';
                        element_warn=1;
                    end;
                    element=strtrim(element);
                end;
                chain_tag=tline(22);
%                 if double(chain_tag)<double('A') || double(chain_tag)>double('Z') % fix wrong chain tags
%                     chain_tag='Z';
%                 end;
                if chain_tag==' '; 
                    if curr_model == 1
                        chain_tag=get_chain_tag(terminated_chains); 
                    else
                        chain_tag = id2tag(cc_id,chain_tags);
                    end
                end;
                c_id=tag2id(chain_tag,chain_tags);
                if isempty(c_id) 
                    chain_id=chain_id+1;
                    structure(chain_id).name=chain_tag;
                    chain_tags=[chain_tags chain_tag ':'];
                    chain_ids=[chain_ids chain_id];
                    structure(chain_id).modified=0;
                    structure(chain_id).nonstandard=0;
                    structure(chain_id).seqtype=0;
                    structure(chain_id).helices=0;
                    structure(chain_id).strands=0;
                    structure(chain_id).xyz{1}=zeros(max_atoms,3);
                    structure(chain_id).isotopes=zeros(max_atoms,2,'single');
                    structure(chain_id).Bfactor{1}=zeros(1,max_atoms);
                    structure(chain_id).Btensor{1}=zeros(max_atoms,6,'int32');
                    structure(chain_id).atoms{1}=0;
                    structure(chain_id).conn=zeros(max_atoms,30,'int32');
                    structure(chain_id).maxconn=1;
                    structure(chain_id).residues{1}.residue_tags=':';
                    structure(chain_id).residues{1}.info=[];
                    c_id=chain_id;
                end;
                cr_tag=strtrim(tline(23:27)); % this includes the infamous insertion code
                id=tag2id(cr_tag,structure(c_id).residues{curr_model}.residue_tags); % check, if this residue is already defined
                curr_resname=tline(18:20);
                if isempty(id) % if not, initialize residue
                    if ~hetflag,
                        absresnum=absresnum+1;
                    end;
                   amino_id=tag2id(tline(18:20),upper(residue_defs.restags),residue_defs.single_letter_code); % test whether this is an amino acid
                   het_amino_id=tag2id(tline(18:20),upper(residue_defs.hettags)); % test whether this is a known non-standard amino acid
                   nucleotide_id=tag2id(tline(18:20),upper(residue_defs.nucleotide_tags)); % test whether this is a known nucleic acid
                   cyana_id=tag2id(tline(18:21),upper(residue_defs.cyana_nucleotide_tags)); % test whether this is a known nucleic acid
                   nuctag = strtrim(tline(18:20));
                   while length(nuctag)<3,
                       nuctag = [' ' nuctag];
                   end;
                   wrong_nuc_id=tag2id(nuctag,upper(residue_defs.nucleotide_tags)); % test whether this is a known nucleic acid
                   if ~isempty(het_amino_id), 
                       seqtype=1; 
                   end;
                   if ~isempty(amino_id), 
                       seqtype=1;
                   end;
                   if ~isempty(cyana_id),
                       nucleotide_id=cyana_id;
                       curr_resname=id2tag(nucleotide_id,upper(residue_defs.nucleotide_tags));
                   end;
                   if ~isempty(wrong_nuc_id),
                       nucleotide_id=wrong_nuc_id;
                       curr_resname=id2tag(nucleotide_id,upper(residue_defs.nucleotide_tags));
                   end;
                   if ~isempty(nucleotide_id),
                       seqtype=2; 
                   end;
                   if chain_type==0,
                       chain_type=seqtype;
                   else
                       if seqtype~=chain_type 
                           if curr_model == 1,
                               add_msg_board('Warning! Chain type mismatch in this structure.');
                               chain_type=seqtype;
                               terminated_chains=[terminated_chains curr_chaintag ':'];
                                chain_tag=tline(22);
                                if chain_tag==' '; 
                                    chain_tag=get_chain_tag(terminated_chains); 
                                end;
                                add_msg_board(sprintf('New chain %s defined',chain_tag)); 
                                c_id=tag2id(chain_tag,chain_tags);
                                if isempty(c_id), % catch cases without SEQRES records
                                    chain_id=chain_id+1;
                                    structure(chain_id).name=chain_tag;
                                    chain_tags=[chain_tags chain_tag ':'];
                                    chain_ids=[chain_ids chain_id];
                                    structure(chain_id).modified=0;
                                    structure(chain_id).nonstandard=0;
                                    structure(chain_id).seqtype=seqtype;
                                    structure(chain_id).helices=0;
                                    structure(chain_id).strands=0;
                                    structure(chain_id).xyz{1}=zeros(max_atoms,3);
                                    structure(chain_id).isotopes=zeros(max_atoms,2,'single');
                                    structure(chain_id).Bfactor{1}=zeros(1,max_atoms);
                                    structure(chain_id).Btensor{1}=zeros(max_atoms,6,'int32');
                                    structure(chain_id).atoms{1}=0;
                                    structure(chain_id).conn=zeros(max_atoms,30,'int32');
                                    structure(chain_id).maxconn=1;
                                    structure(chain_id).residues{1}.residue_tags=':';
                                    structure(chain_id).residues{1}.info=[];
                                    c_id=chain_id;
                                    atnum=structure(c_id).atoms{curr_model}+1;
                                end;
                               chain_type=seqtype;
                           else
                               cc_id = cc_id + 1;
                               chain_tag = id2tag(cc_id,chain_tags);
                               chain_type=seqtype;
                               c_id = cc_id;
                           end
                       end;
                   end;
                   curr_chaintag=chain_tag;
                   structure(c_id).residues{curr_model}.residue_tags=[structure(c_id).residues{curr_model}.residue_tags cr_tag ':'];
                   id=length(structure(c_id).residues{curr_model}.info)+1;
                   structure(c_id).residues{curr_model}.info(id).name=curr_resname;
                   structure(c_id).residues{curr_model}.info(id).absresnum=absresnum;
                   structure(c_id).residues{curr_model}.info(id).type=0;
                   if ~isempty(amino_id) && ~strcmpi(record,'HETATM'), 
                       structure(c_id).residues{curr_model}.info(id).type=1; 
                       structure(c_id).seqtype=1; 
                   end;
                   if ~isempty(het_amino_id) && ~strcmpi(record,'HETATM'), 
                       structure(c_id).residues{curr_model}.info(id).type=1; 
                       structure(c_id).seqtype=1; 
                   end;
                   if ~isempty(nucleotide_id) && ~strcmpi(record,'HETATM'), 
                       structure(c_id).residues{curr_model}.info(id).type=2; 
                       structure(c_id).seqtype=2; 
                   end;
                   sec=0;
                   for ksec=1:structure(c_id).helices,
                       if curr_res>=structure(c_id).helix(ksec).start && curr_res<=structure(c_id).helix(ksec).end
                           sec=1;
                       end;
                   end;
                   for ksec=1:structure(c_id).strands,
                       if curr_res>=structure(c_id).strand(ksec).start && curr_res<=structure(c_id).strand(ksec).end
                           sec=2;
                       end;
                   end;
                   structure(c_id).residues{curr_model}.info(id).secondary=sec;
                   structure(c_id).residues{curr_model}.info(id).hetflag=hetflag;
                   structure(c_id).residues{curr_model}.info(id).connected=0;
                   structure(c_id).residues{curr_model}.info(id).number=curr_res;
                   structure(c_id).residues{curr_model}.info(id).atom_tags=':';
                   structure(c_id).residues{curr_model}.info(id).elements=[];
                   structure(c_id).residues{curr_model}.info(id).atom_numbers={};
                   structure(c_id).residues{curr_model}.info(id).residue_graphics=[];
                   structure(c_id).residues{curr_model}.info(id).label_graphics=[];
                   structure(c_id).residues{curr_model}.info(id).location_tags=':A:';
                   structure(c_id).residues{curr_model}.info(id).insertion_code=tline(27);
                   if ~hetflag, total_residues=total_residues+1; end;
                   if inserted,
                       ins_poi=ins_poi+1;
                       insertions(ins_poi,:)=[c_id curr_model id];
                   end;
                   if strcmp(tline(13:15),'OXT'),
                       chain_type=999; % guard against unterminated chains
                   end;
                end;
                % store the atom coordinates and B-factor in the appropriate array
                atnum=structure(c_id).atoms{curr_model}+1;
                structure(c_id).atoms{curr_model}=atnum;
                x=sscanf(tline(31:38),'%f');
                y=sscanf(tline(39:46),'%f');
                z=sscanf(tline(47:54),'%f');
                structure(c_id).xyz{curr_model}(atnum,:)=[x,y,z];
                center=center+[x,y,z];
                atoms=atoms+1;
                elnum=tag2id(upper(element),upper(chemistry.element_tags));
                if isempty(elnum), 
                    elnum=0;
                    mass=0;
                else
                    mass=chemistry.pse(elnum).mass;
                end;
                structure(c_id).isotopes(atnum,1)=elnum;
                structure(c_id).isotopes(atnum,2)=mass;
                pdbatnum=sscanf(tline(7:11),'%i');
                atom_corr(pdbatnum,:)=[at,c_id,atnum]; % correspondence of PDB atom number to chain ID/atom number
                if length(tline)>=66,
                    Bfactor=sscanf(tline(61:66),'%f');
                else
                    Bfactor=0;
                end;
                if isempty(Bfactor),
                    Bfactor=0;
                end;
                structure(c_id).Bfactor{curr_model}(atnum)=Bfactor;
                if Bfactor>B_range(2), B_range(2)=Bfactor; end;
                if Bfactor<B_range(1), B_range(1)=Bfactor; end;    
                if inserted
                    seq_found(c_id)=0; % sequence information has to be redone, as it will not match residue numbers
                end
                % store the atom tag and atom number in the
                at_id=tag2id(atom_tag,structure(c_id).residues{curr_model}.info(id).atom_tags); % check, if atom already exists (alternate positions)
                if isempty(at_id) % if not, define new atom, store its tag and element, and initialize the atom number list
                    at_id=length(structure(c_id).residues{curr_model}.info(id).atom_numbers)+1;
                    structure(c_id).residues{curr_model}.info(id).atom_tags=[structure(c_id).residues{curr_model}.info(id).atom_tags atom_tag ':'];
                    elnum=tag2id(upper(element),upper(chemistry.element_tags));
                    if isempty(elnum),
                        elnum=0;
                    end;
                    structure(c_id).residues{curr_model}.info(id).elements(at_id)=elnum;
                    structure(c_id).residues{curr_model}.info(id).atom_numbers{at_id}=[];
                end;
                if length(tline)>=60
                    occupancy=sscanf(tline(55:60),'%f');
                else
                    occupancy=1;
                end,
                if isempty(occupancy), occupancy=1; end;
                if abs(occupancy-1)<occupancy_error || strcmp(tline(17),' '), % defined location of this atom
                    structure(c_id).residues{curr_model}.info(id).atom_numbers{at_id}=atnum; % residue must know where its atom coordinates are
                else % alternate locations of this atom
                    alternate=true;
                    if ~isfield(structure(c_id).residues{curr_model}.info(id),'location_tags') || length(structure(c_id).residues{curr_model}.info(id).location_tags)<2, % no location tags are yet defined
                        structure(c_id).residues{curr_model}.info(id).location_tags=[':' tline(17) ':'];
                        location_id=1;
                    else
                        location_id=tag2id(tline(17),structure(c_id).residues{curr_model}.info(id).location_tags); % check, if location tag already exists
                        if isempty(location_id), % if not, generate it and request id again
                            structure(c_id).residues{curr_model}.info(id).location_tags=[structure(c_id).residues{curr_model}.info(id).location_tags tline(17) ':'];
                            location_id=tag2id(tline(17),structure(c_id).residues{curr_model}.info(id).location_tags);
                        end;
                    end;
                    % check whether another location was already stored
                    if isempty(structure(c_id).residues{curr_model}.info(id).atom_numbers{at_id}) % new atom
                        structure(c_id).residues{curr_model}.info(id).atom_numbers{at_id}=[atnum,occupancy,location_id];
                    else % alternate location for an existing atom
                        structure(c_id).residues{curr_model}.info(id).atom_numbers{at_id}=[structure(c_id).residues{curr_model}.info(id).atom_numbers{at_id};[atnum,occupancy,location_id]];
                    end;
                end;
               if strcmpi(char(tline(13:15)),'OXT'),
                   chain_type=999; % guard against unterminated chains
               end;
            case {'ANISOU'}
                if length(tline)>=70
                  PDB_sernum=sscanf(tline(7:11),'%i');
                  c_id=atom_corr(PDB_sernum,2);
                  at_id=atom_corr(PDB_sernum,3);
                  Btensor(1)=sscanf(tline(29:35),'%i');
                  Btensor(2)=sscanf(tline(36:42),'%i');
                  Btensor(3)=sscanf(tline(43:49),'%i');
                  Btensor(4)=sscanf(tline(50:56),'%i');
                  Btensor(5)=sscanf(tline(57:63),'%i');
                  Btensor(6)=sscanf(tline(64:70),'%i');
                  structure(c_id).Btensor{curr_model}(at_id,:)=Btensor;
                end;
            case 'TER   '
                if ~exist('c_id','var') || isempty(c_id),
                    chain_tag=tline(22);
%                     if double(chain_tag)<double('A') || double(chain_tag)>double('Z') % fix wrong chain tags
%                         chain_tag='A';
%                     end;
                    if chain_tag==' ', chain_tag=get_chain_tag(terminated_chains); end;
                    c_id=tag2id(chain_tag,chain_tags);
                end;
                if ~isempty(c_id), % catch cases without SEQRES records
                    % residue information, occupancy is also stored, when needed
                    curr_res=sscanf(tline(23:26),'%i');
                    cr_tag=sprintf('%i',curr_res);
                    id0=id;
                    id=tag2id(cr_tag,structure(c_id).residues{curr_model}.residue_tags); % check, if this residue is already defined
                    if isempty(id), id=id0; end; % catches wrong TER records, for instance of NMsim
                    structure(c_id).residues{curr_model}.info(id).terminal=1;
                    terminated_chains=[terminated_chains chain_tag ':'];
                    chain_type=0;
                end;
            case 'CONECT'
                curr_atom=sscanf(tline(7:11),'%i');
                c_id=atom_corr(curr_atom,2);
                c_atnum=atom_corr(curr_atom,3);
                if c_id~=0 && c_atnum~=0,
                    bonds=structure(c_id).conn(c_atnum,:);
                    maxbonds=sum(bonds~=0);
                    basis=7;
                    k=1;
                    while length(tline)>=5*k+basis+4, % check if record line is long enough
                        newbond=sscanf(tline(basis+5*k:basis+5*k+4),'%i');
                        k=k+1;
                        if ~isnan(newbond),
                            b_atnum=atom_corr(newbond,3);
                            if b_atnum~=0 && ~isempty(bonds),
                                if isempty(find(bonds==b_atnum,1)) % check if this bond was already stored
                                    maxbonds=maxbonds+1;
                                    if maxbonds > structure(c_id).maxconn % update maximum size of connection table
                                        structure(c_id).maxconn=maxbonds;
                                    end
                                    bonds(maxbonds)=b_atnum;
                                end
                            end
                        end
                    end
                    structure(c_id).conn(c_atnum,:)=bonds;
                end
        end
    end
    % flush sequence stack
    structure(chain_id).sequence=seq;
    structure(chain_id).seqres=seq;
    structure(chain_id).restags=restags;
    
    % process data base reference for sequences
    
    chains=textscan(chain_tags(2:end),'%s','Delimiter',':');
    chains=chains{1};
    for k=1:length(chains),
        tag=char(chains(k));
        seq0=structure(k).sequence;
        id=tag2id(tag,db_tags);
        header='';
        sequence='';
        if isfield(structure(k),'seqtype') && structure(k).seqtype==1 && ~isempty(id) && length(db_access)>=id, % this is done only for amino acid sequences
            access=db_access{id};
            structure(k).dbref=access;
            [dbase,code]=strtok(access,':');
            if strcmpi(dbase,'UNP'), % for UniProt, sequences are read
                query=[web_adr.UniProt code(2:end) queries.UniProt_fasta];
                fname=strcat(general.tmp_files,'current_sequence.fasta');
                [f,status]=urlwrite(query,fname);
                if status,
                    [header,sequence]=read_fasta(fname);
                    if db_cutout(id,1)>0 && db_cutout(id,2)~=0 && db_cutout(id,1)<=length(sequence) && db_cutout(id,2)<=length(sequence),
                        sequence=sequence(db_cutout(id,1):db_cutout(id,2));
                        if isfield(structure(k),'mutations'),
                            for kk=1:length(structure(k).mutations),
                                res=structure(k).mutations(kk).number;
                                tlc=structure(k).mutations(kk).modified;
                                slcid=tag2id(tlc,upper(residue_defs.restags),residue_defs.single_letter_code);
                                if isempty(slcid),
                                    slcid='?';
                                end;
                                sequence(res)=slcid;
                            end;
                        end;
                        if db_cutout(id,1)>1,
                            for kss=2:db_cutout(id,1),
                                sequence=['?' sequence];
                            end;
                        end;
                    else
                        sequence='';
                    end;
                    structure(k).fasta=sequence;
                end;
            end;
        end;
        structure(k).header=header;
%         if ~isempty(sequence),
%             if length(seq0)>length(sequence),
%                 sequence=[sequence seq0(1+length(sequence):end)];
%             end;
%             structure(k).sequence=sequence;
%         end;
    end;

    num_chains=length(structure);
    if ~rotamers, % suppose backbone connection in rotamer pseudo-format
        % make connectivity for standard backbone and intramolecular
        % connectivity for standard residues
        for k=1:num_chains,
            residue_tags=structure(k).residues{1}.residue_tags;
            info=structure(k).residues{1}.info;
            resnum=length(info);
            for kk=1:resnum,
                current=info(kk).number;
                % backbone connections for amino acids, defensive programming
                if info(kk).type==1, % do only for amino acids, CONECT records take care of hetero residues
                    previous_tag=sprintf('%i',current-1);
                    previous=tag2id(previous_tag,residue_tags);
                    if ~isempty(previous) && info(previous).type==1, % if previous residue is found and is an amino acid
                        curr_N_id=tag2id('N',info(kk).atom_tags);
                        if ~isempty(curr_N_id),
                            curr_N=info(kk).atom_numbers{curr_N_id};
                            [m1,n1]=size(curr_N); % necessary because of alternate positions
                            prev_C_id=tag2id('C',info(previous).atom_tags);
                            if ~isempty(prev_C_id),
                                prev_C=info(previous).atom_numbers{prev_C_id};
                                [m2,n2]=size(prev_C); % necessary because of alternate positions
                                for o=1:m1,
                                    for oo=1:m2,
                                        if n1<3 || n2<3 || curr_N(o,3)==prev_C(oo,3) % this assumes that location tags correspond for the two residues
                                            structure(k)=make_bond(structure(k),curr_N(o,1),prev_C(oo,1));
                                            structure(k)=make_bond(structure(k),prev_C(oo,1),curr_N(o,1));
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
                    next_tag=sprintf('%i',current+1);
                    next=tag2id(next_tag,residue_tags);
                    if ~isempty(next) && info(next).type==1, % if next residue is found and is an amino acid
                        curr_C_id=tag2id('C',info(kk).atom_tags);
                        if ~isempty(curr_C_id),
                            curr_C=info(kk).atom_numbers{curr_C_id};
                            [m1,n1]=size(curr_C); % necessary because of alternate positions
                            next_N_id=tag2id('N',info(next).atom_tags);
                            if ~isempty(next_N_id),
                                next_N=info(next).atom_numbers{next_N_id};
                                [m2,n2]=size(next_N); % necessary because of alternate positions
                                for o=1:m1,
                                    for oo=1:m2,
                                        if n1<3 || n2<3 || curr_C(o,3)==next_N(oo,3)  % this assumes that location tags correspond for the two residues
                                            structure(k)=make_bond(structure(k),curr_C(o,1),next_N(oo,1));
                                            structure(k)=make_bond(structure(k),next_N(oo,1),curr_C(o,1));
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
                    % intramolecular connections for amino acids
                    structure(k)=mk_internal_bonds(structure(k),kk,residue_defs);
                end;
                if info(kk).type==2, % do only for nucleic acids, CONECT records take care of hetero residues
                    next_tag=sprintf('%i',current+1);
                    next=tag2id(next_tag,residue_tags);
                    if ~isempty(next) && info(next).type==2, % if next residue is found and is a nucleic acid
                        curr_C_id=tag2id('O3''',info(kk).atom_tags);
                        if isempty(curr_C_id),
                            curr_C_id=tag2id('O3*',info(kk).atom_tags);
                        end;
                        curr_C=info(kk).atom_numbers{curr_C_id};
                        [m1,n1]=size(curr_C); % necessary because of alternate positions
                        next_N_id=tag2id('P',info(next).atom_tags);
                        if ~isempty(next_N_id),
                            next_N=info(next).atom_numbers{next_N_id};
                            [m2,n2]=size(next_N); % necessary because of alternate positions
                            for o=1:m1,
                                for oo=1:m2,
                                    if n1<3 || n2<3 || curr_C(o,3)==next_N(oo,3)  % this assumes that location tags correspond for the two residues
                                        structure(k)=make_bond(structure(k),curr_C(o,1),next_N(oo,1));
                                        structure(k)=make_bond(structure(k),next_N(oo,1),curr_C(o,1));
                                    end;
                                end;
                            end;
                        end;
                    end;
                    % intramolecular connections for nucleic acids)
                    structure(k)=mk_internal_bonds(structure(k),kk,residue_defs,true);
                end;
            end;
        end;
    end;

    % reduce data size to used array sizes
    for k=1:num_chains,
        structure(k).conn=structure(k).conn(1:structure(k).atoms{1},1:structure(k).maxconn);
        structure(k).isotopes=structure(k).isotopes(1:structure(k).atoms{1},:);
        for kk=1:structures,
            structure(k).xyz{kk}=structure(k).xyz{kk}(1:structure(k).atoms{kk},:);
            structure(k).Bfactor{kk}=structure(k).Bfactor{kk}(1:structure(k).atoms{kk});
            structure(k).Btensor{kk}=structure(k).Btensor{kk}(1:structure(k).atoms{kk},:);
        end;
    end;
    remarks=remarks(1:rem_num);
    fclose(fid);

    % initialize sequence information if this was not present in PDB file or if
    % the SEQRES record does not match sequence

    for cid=1:num_chains,
        info=structure(cid).residues{1}.info;
        if seq_found(cid), % check for sequence integrity
            seq=structure(cid).sequence;
            nres=length(info);
            for kres=1:nres,
                icode=strtrim(info(kres).insertion_code);
                if ~isempty(icode), % insertion code? sequence is crap!
                    seq_found(cid)=0;
                    break;
                end;
                tag=info(kres).name;
                num=info(kres).number;
                if num<=0,
                    seq_found(cid)=0;
                    break;
                end;
                if structure(cid).seqtype==1,
                    id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
                    if isempty(id), id='?'; end;
                elseif structure(cid).seqtype==2,
                    id=tag2id(tag,upper(residue_defs.nucleotide_tags),residue_defs.nucleotide_slc);
                    if isempty(id), id='?'; end;
                else
                    id='?';
                end;
                if num>length(seq),
                    seq_found(cid)=false;
                    break;
                end;
                if ~strcmp(id,seq(num)),
                    seq(num)=id;
                end;
                if isfield(info(kres),'terminal'),
                    if length(info)>kres,
                        goon=info(kres+1).type;
                    else
                        goon=0;
                    end;
                    if info(kres).terminal==1,
                        if goon~=1,
                            break;
                        else
                            structure(cid).residues{1}.info(kres).terminal=[];
                        end;
                    end;
                end;
            end;
            structure(cid).sequence=seq;
        end;
        if ~seq_found(cid),
            seq='';
            restags=':';
            info=structure(cid).residues{1}.info; % coordinate set {1}, as all have the same sequence
            nres=length(info);
            act=1;
            for kres=1:nres,
                icode=strtrim(info(kres).insertion_code);
                tag=info(kres).name;
                if isempty(icode), % residues with insertion code are not shown in sequence
                    while act<info(kres).number, % care about gaps
                        seq=[seq '?'];
                        act=act+1;
                    end;
                    if structure(cid).seqtype==1,
                        id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
                        if isempty(id), id='?'; end;
                    elseif structure(cid).seqtype==2,
                        id=tag2id(tag,upper(residue_defs.nucleotide_tags),residue_defs.nucleotide_slc);
                        if isempty(id), id='?'; end;
                    else
                        id='?';
                    end;
                    seq=[seq id];
                    act=act+1;
                    if isfield(info(kres),'terminal'),
                        if info(kres).terminal==1,
                            if kres==nres,
                                restags=[restags tag ':'];
                                break;
                            elseif info(kres+1).type~=1,
                                restags=[restags tag ':'];
                                break;
                            else
                                info(kres).terminal=[];
                            end;
                        end;
                    end;
                end;
                restags=[restags tag ':'];
            end;
            structure(cid).sequence=seq;
            structure(cid).restags=restags;
            structure(cid).residues{1}.info=info;
        end;
    end;

    if element_warn,
        add_msg_board('Element column missing and probably wrong element format');
        add_msg_board('in columns 13,14 (tried to fix it, assuming only H,N,C,O,P,S)');
    end;
    clear info % variable info is used in two ways, sorry

    % repair shoddy PDB references, if there is web access and PubMed IDs are
    % known, 
    % this is done one by one, which is slower, but safer
    if ~isempty(references),
        for k=1:length(references),
            if references(k).PMID~=0,
                query=sprintf('%s%i%s',web_adr.PubMed,references(k).PMID,queries.PMID);
                fname=strcat(general.tmp_files,queries.PubMed_file);
                [f,status]=urlwrite(query,fname);
                if status,
                    newref=rd_MEDLINE(f,1);
                    if ~isempty(newref)
                        newref(1).short=references(k).short;
                        if ~isempty(references(k).DOI),
                            newref(1).DOI=references(k).DOI;
                        end;
                        references(k)=newref(1);
                    end;
                end;
            end;
        end;
    end;

    info.class=pdb_class;
    info.idCode=idCode;
    info.depDate=depDate;
    info.title=strtrim(title);
    info.remarks=remarks;
    info.chain_tags=chain_tags;
    info.chain_ids=chain_ids;
    info.warning=warning;
    info.center=center/atoms;
    info.atoms=atoms;
    info.residues=total_residues;
    info.B_range=B_range;
    info.missing=missing;
    info.references=references;
    info.keywords=keywords;
    info.authors=authors;
    info.metal=metal;
    info.rotamers=rotamers;
    if exist('cryst','var'),
        info.cryst=cryst;
    end;
    info.molecule=molecule;
    info.organism=organism;
    info.Modeller_obj=Modeller_obj;
    info.Modeller_sid=Modeller_sid;

    info.SSbonds=SSbonds;    

    for k=1:length(sites),
        tag=id2tag(k,site_tags);
        id=tag2id(tag,sites_declared);
        if ~isempty(id) && id<=length(site_evidence),
            sites(k).evidence=site_evidence{id};
        else
            sites(k).evidence='unknown';
        end;
        if ~isempty(id) && id<=length(site_declarations),
            sites(k).declaration=site_declarations{id};
        else
            sites(k).declaration='unknown';
        end;
        sites(k).evidence_type=1; % atomistic structure
    end;

    info.site_tags=site_tags;
    info.sites=sites;

    if ins_poi>0,
        insertions=insertions(1:ins_poi,:);
    else
        insertions=[];
    end;
    info.insertions=insertions;
    info.alternate=alternate;
    info.resolution=resolution;

    if nargin>1,
        fclose(rem_out);
    end;


    function myCleanup

        if exist('fid','var') && sum(fid==fopen('all'))>0,
            fclose(fid);
        end;
    end
    % disp(sprintf('%i sequence lines, %i atom lines, thereof %i hetatom lines',sr,at,het));
    
    function chain_tag=get_chain_tag(terminated_chains)
        chain_tag='A';
        exists=true;
        while exists && double(chain_tag)<=double('Z'),
            id=tag2id(chain_tag,terminated_chains);
            if isempty(id),
                exists=false;
            else
                chain_tag=char(double(chain_tag)+1);
            end;
        end;
        if exists,
            chain_tag=[];
        end;
    end
    
end
