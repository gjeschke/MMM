function restraints=rd_restraints_docking(fname)
% function restraints=rd_restraints(fname)
%
% reads mixed experimental restraints from an ASCI file
% restraints are returned in a structure whose fields correspond to
% different restraint types
%
% currently implemented types are:
%
% DEER      distance between two labelled sites
% DIRECT    distance between C_alpha atoms

default_helix_mindist=0.65; % default value for the closest approach of two helix axes

restraints=[];

fid=fopen(fname);
if fid==-1,
    add_msg_board('ERROR: Constraint file does not exist');
    return;
end;

clear restraints
restraints.ensemble=1;
restraints.uncertainty=0;
restraints.exclude=true;
restraints.realign=false;
restraints.target='';
restraints.helix_mindist=default_helix_mindist;
restraints.ref_chain='';

DEER_poi=0;
direct_poi=0;
align_poi=0;
helix_poi=0;
strand_poi=0;
sheet_poi=0;
disp_poi=0;
locate_poi=0;
network_poi=0;
mode=0;
scale_units=1;

while 1
    tline = fgetl(fid);
    if ~ischar(tline) || mode<0, break, end
    if ~isempty(tline),
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k),
            if k(1)>1,
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end;
        end;
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#'),
            switch upper(char(args(2)))
                case 'PATH'
                    mode = 0;
                    datapath = char(args(3));
                    genpath(datapath);
                case 'TRANSITION'
                    mode=0;
                    restraints.initial=char(args(3));
                    restraints.final=char(args(4));
                case 'PDB' % -------- modified by yepo 190713--------------
                    mode=0;
                    for k=3:length(args),
%                         pdbid=strtok(char(args(3)),':');
                        [pdbid,chains]=strtok(char(args(k)),':');
                        restraints.PDB(k-2).pdbid=pdbid;
                    
%                       restraints.PDB=pdbid; % to conform to old format
%                     for k=3:length(args),
%                         [pdbid,chains]=strtok(char(args(k)),':');
                        restraints.template(k-2).pdbid=pdbid;
                        if length(chains)>1,
                            nonsense=textscan(chains(2:end),'%s','Delimiter',',');
                            chains=nonsense{1};
                            chaintags=':';
                            for kk=1:length(chains),
                                chaintags=[chaintags char(chains(kk)) ':'];
                            end;
                            restraints.template(k-2).chaintags=chaintags;
                        else
                            restraints.template(k-2).chaintags='';
                        end;
                    end; % ------------------------------------------------
                    
                    % % % to be extended for reading two PDB's!!!
                    % % % think of homodim cases!!!!!
                    
                case {'ALIGN','ALIGNMENT'}
                    mode=0;
                    restraints.alignment=char(args(3));
                case 'BUNDLE'
                    mode=0;
                    bundle=zeros(1,length(args)-2);
                    for k=3:length(args),
                        bundle(k-2)=str2double(char(args(k)));
                    end;
                    restraints.bundle=bundle;
                case 'CHAINS'
                    mode=0;
                    chains=':';
                    for k=3:length(args),
                        nonsense=textscan(char(args(k)),'%s','Delimiter',',');
                        chainlist=nonsense{1};
                        for kk=1:length(chainlist),
                            chains=[chains char(chainlist(kk)) ':'];
                        end;
                    end;
                    restraints.chains=chains;
                case 'ENSEMBLE'
                    mode=0;
                    restraints.ensemble=str2double(char(args(3)));
                    if length(args)>3,
                        restraints.uncertainty=scale_units*str2double(char(args(4)));
                    else
                        restraints.uncertainty=0;
                    end;
                    if length(args)>4,
                        if strcmpi(char(args(5)),'all')
                            restraints.exclude=false;
                        else
                            restraints.exclude=true;
                        end;
                    else
                        restraints.exclude=true;
                    end;
                case 'DEER'
                    mode=1;
                    label=char(args(3));
                    T=char(args(4));
                case 'DIRECT'
                    mode=2;
                case 'LOCATE'
                    mode=8;
                    if length(args)>2,
                        tag=char(args(3));
                    else
                        tag='loc';
                    end;
                    if length(args)>3,
                        probability=str2double(char(args(4)));
                    else
                        probability=0.5;
                    end;
                    if length(args)>4,
                        label=char(args(5));
                    else
                        label='MTSL';
                    end;
                    if length(args) > 5,
                        T=char(args(6));
                    else
                        T='298';
                    end;
                case 'NETWORK'
                    mode=9;
                    if length(args)>2,
                        label=char(args(3));
                    else
                        label='MTSL';
                    end;
                    if length(args) > 3,
                        T=char(args(4));
                    else
                        T='298';
                    end;
                    if length(args)>4,
                        probability=str2double(char(args(5)));
                    else
                        probability=0.5;
                    end;
                    if length(args)>5,
                        display_mode=char(args(6));
                    else
                        display_mode='full';
                    end;
                case 'REALIGN'
                    mode=3;
                    if length(args)>2,
                        restraints.align_threshold=char(args(3));
                    else
                        restraints.align_threshold=0;
                    end;
                    restraints.realign=true;
                case 'HELICES'
                    mode=4;
                    if length(args)>2, % minimum helix axes distance is provided
                        restraints.helix_mindist=scale_units*str2double(char(args(3)));
                    end;
                case 'STRANDS'
                    mode=5;
                case 'SHEETS'
                    mode=6;
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
                    end;
                case 'SYMMETRY'
                    mode=0;
                    restraints.Cn=str2double(char(args(3)));
                    restraints.ref_chain=char(args(4));
                    for k=5:length(args),
                        restraints.chain_copies{k-4}=char(args(k));
                    end;
                case 'PARALLEL'
                    mode=0;
                    for k=3:length(args),
                        restraints.parallel(k-2)=str2double(char(args(k)));
                    end;
                case 'PERPENDICULAR'
                    mode=0;
                    for k=3:length(args),
                        restraints.perpendicular(k-2)=str2double(char(args(k)));
                    end;
                case 'EXTEND'
                    mode=0;
                    restraints.initial_model=char(args(3)); 
                case 'MAXTIME'
                    mode=0;
                    restraints.max_time=str2double(char(args(3))); 
                case 'TOLERANCE'
                    mode=0;
                    restraints.tolerance=str2double(char(args(3))); 
                case 'END'
                    mode=-1;
                otherwise
                    mode=0;
                    add_msg_board('Warning: Unknown restraint mode');
            end;
        elseif mode>0 && ~strncmp(strtrim(char(args(1))),'%',1),
% %             keyboard
            switch mode
                case 1
                    DEER_poi=DEER_poi+1;
                    switch length(restraints.PDB)
                        case 1 % (homodimer or homoologomer)
                            restraints.DEER(DEER_poi).label=label;
                            restraints.DEER(DEER_poi).T=str2double(T);
                            restraints.DEER(DEER_poi).adr1=char(args(1));
                            restraints.DEER(DEER_poi).r=scale_units*str2double(char(args(2)));
                            restraints.DEER(DEER_poi).sigr=scale_units*str2double(char(args(3)));
                        case 2 % (heterodimer)
                            restraints.DEER(DEER_poi).label=label;
                            restraints.DEER(DEER_poi).T=str2double(T);
                            restraints.DEER(DEER_poi).adr1=char(args(1));
                            restraints.DEER(DEER_poi).adr2=char(args(2));
                            restraints.DEER(DEER_poi).r=scale_units*str2double(char(args(3)));
                            restraints.DEER(DEER_poi).sigr=scale_units*str2double(char(args(4)));
                    end
                case 2
                    direct_poi=direct_poi+1;
                    switch length(restraints.PDB)
                        case 1
                            restraints.direct(direct_poi).adr1=char(args(1));
                            restraints.direct(direct_poi).adr2=char(args(2));
                            restraints.direct(direct_poi).r=scale_units*str2double(char(args(3)));
                        case 2
                            restraints.direct(direct_poi).adr1=char(args(1));
                            restraints.direct(direct_poi).adr2=char(args(2));
                            restraints.direct(direct_poi).r=scale_units*str2double(char(args(3)));
                            restraints.direct(direct_poi).sigr=scale_units*str2double(char(args(4)));
                    end
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
                    if ~strcmpi(weights,'uniform') && length(args)>=5,
                        restraints.displace(disp_poi).weight=str2double(char(args(5)));
                    else
                        restraints.displace(disp_poi).weight=1;
                    end;
                case 8 % localization
                    locate_poi=locate_poi+1;
                    restraints.locate(locate_poi).tag=tag;
                    restraints.locate(locate_poi).probability=probability;
                    restraints.locate(locate_poi).adr=char(args(1));
                    restraints.locate(locate_poi).r=scale_units*str2double(char(args(2)));
                    restraints.locate(locate_poi).sigr=scale_units*str2double(char(args(3)));
                    if length(args)>3,
                        restraints.locate(locate_poi).type=char(args(4));
                    else
                        restraints.locate(locate_poi).type='DEER';
                    end;
                    if length(args)>4,
                        restraints.locate(locate_poi).label=char(args(5));
                    else
                        restraints.locate(locate_poi).label='MTSL';
                    end;
                    if length(args)>5,
                        restraints.locate(locate_poi).T=str2double(char(args(6)));
                    else
                        restraints.locate(locate_poi).T=298;
                    end;
                    if length(args)>6,
                        restraints.locate(locate_poi).file=char(args(7));
                        restraints.locate(locate_poi).fulldistr=true;
                    else
                        restraints.locate(locate_poi).fulldistr=false;
                    end;
                case 9 % network
                    network_poi=network_poi+1;
                    restraints.network(network_poi).label=label;
                    restraints.network(network_poi).T=str2double(T);
                    restraints.network(network_poi).probability=probability;
                    restraints.network(network_poi).adr1=char(args(1));
                    restraints.network(network_poi).adr2=char(args(2));
                    restraints.network(network_poi).r=scale_units*str2double(char(args(3)));
                    restraints.network(network_poi).sigr=scale_units*str2double(char(args(4)));
                    restraints.network(network_poi).display_mode=display_mode;
                    if length(args)>4,
                        restraints.network(network_poi).file=char(args(5));
                        restraints.network(network_poi).fulldistr=true;
                    else
                        restraints.network(network_poi).fulldistr=false;
                    end;
            end;
        end;
    end;
end;

fclose(fid);