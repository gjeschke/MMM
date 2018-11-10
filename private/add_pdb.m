function [message,snum]=add_pdb(fname,idCode)
% function [message,snum]=add_pdb(fname,idCode)
%
% Adds a PDB structure to an MMM model by calling rd_pdb, the PDB
% identifier, if present, is defined as a synonym for the structure
% if no model exists, a new model is created
%
% fname         input file name
% idCode        optional PDB-like identifier code, is read from input file
%               if not given
% message       error message structure with fields
%               .error 0: no error, 1: PDB file does not exist, 2: PDB file
%               corrupted
%               .text
% snum          number of the added structure in the model
%
% G. Jeschke, 2009

global model
global hMain
global general

message.error=0;
message.text='';
snum=[]; % empty structure number, if it fails

if isempty(model) || ~isfield(model,'structures') || isempty(model.structures),
    snum=1; % this is the first structure
    model.structure_tags=':';
    model.structure_ids=[];
    model.unrecord=gobjects(0,0);
    model.annotations=[];
else
    snum=length(model.structures)+1;
end

if ischar(fname)
    [model.structures{snum},info]=rd_pdb(fname,strcat(general.tmp_files,'remarks.txt'));
else
    model.structures{snum} = fname;
    info = idCode;
    idCode = info.idCode;
end
if isempty(info) || isempty(model.structures{snum}) % PDB read was not successful
    if isfield(info,'no_file')
        message.error=1;
        message.text='PDB file does not exist.';
        return
    else
        message.error=2;
        message.text='PDB file corrupted.';
        return
    end    
end
set(hMain.menu_file_save_PDB,'Enable','on');
% redistribute information
model.chain_tags{snum}=info.chain_tags;
model.chain_ids{snum}=info.chain_ids;
model.info{snum}.title=info.title;
model.info{snum}.molecule=lower(info.molecule);
model.info{snum}.organism=lower(info.organism);
model.info{snum}.organism(1)=upper(model.info{snum}.organism(1));
model.info{snum}.remarks=info.remarks;
model.info{snum}.class=info.class;
model.info{snum}.depDate=info.depDate;
if nargin<2,
    model.info{snum}.idCode=info.idCode;
else
    model.info{snum}.idCode=idCode;
end;
model.info{snum}.center=info.center;
model.info{snum}.atoms=info.atoms;
model.info{snum}.residues=info.residues;
model.info{snum}.B_range=info.B_range;
model.info{snum}.SSbonds=info.SSbonds;
model.info{snum}.Modeller_obj=info.Modeller_obj;
model.info{snum}.Modeller_sid=info.Modeller_sid;
if isfield(info,'cryst'),
    model.info{snum}.cryst=info.cryst;
end;
if ~isempty(strtrim(info.idCode)) || nargin>1
    stag=model.info{snum}.idCode;
    id=tag2id(stag,model.structure_tags);
    poi=1;
    while ~isempty(id),
        stag=sprintf('%s_%i',model.info{snum}.idCode,poi);
        poi=poi+1;
        id=tag2id(stag,model.structure_tags);
    end;
    model.structure_tags=sprintf('%s%s:',model.structure_tags,stag);
else
    if info.rotamers,
        resadr=strtok(info.title);
        stag0='';
        for k=1:length(resadr), % replace all reserved characters by underscore
            switch resadr(k),
                case {'[',']','{','}','(',')','.',':'}
                    if ~isempty(stag0) && stag0(end)~='_',
                        stag0=[stag0 '_'];
                    end;
                otherwise
                    stag0=[stag0 resadr(k)];
            end;
        end;
        stag=sprintf('%s_rotamers',stag0);
        model.structure_tags=sprintf('%s%s:',model.structure_tags,stag);
    else
        model.structure_tags=sprintf('%s%i:',model.structure_tags,snum);
        stag=sprintf('%i',snum);
    end;
end;
model.structure_ids=[model.structure_ids snum];
model.current_structure=snum;
model.current_chain=id2tag(1,model.chain_tags{snum});
dssp=[];
dospath=which('dssp-2.0.4-win32.exe');
if isempty(dospath),
    dospath=which('dssp.exe');
end;
if isempty(dospath),
    dospath=which('dsspcmbi.exe');
end;
if ~isempty(dospath) && ~info.rotamers && ischar(fname), % suppress this if DSSP not known or MMM rotamer format or direct input
    add_msg_board('DSSP geometry analysis is performed');
    infile=which(fname);
    poi=strfind(infile,'.pdb');
    if isempty(poi),
        poi=strfind(infile,'.ent');
    end;
    if isempty(poi) || poi(1)<2,
        basname=infile;
    else
        basname=infile(1:poi(1)-1);
    end;
    outfile=[basname '.dssp'];
    cmd=[dospath ' ' infile ' ' outfile];
    [s, w] = dos(cmd);
    save(fullfile(general.tmp_files,'dssp_diagnostics.mat'),'s','w','cmd');
    if s~=0,
        add_msg_board('Warning: DSSP did not run successfully on this PDB file.');    
        % disp(w);
    else
        add_msg_board('Result file:');
        add_msg_board(outfile);
        dssp_file=outfile;
        dssp=rd_dssp(dssp_file);
        post_process_dssp(dssp,snum);
    end;
end;
% initialize MMM secondary structure definitions
chains=length(info.chain_ids);
for cnum=1:chains,
    pdb_secondary(snum,cnum);
end;
incomplete=[];
if ~isempty(info.missing),
    incomplete=zeros(length(info.missing),4);
    poi=0;
    for k=1:length(info.missing),
        info.missing{k}=sprintf('[%s]%s',stag,info.missing{k});
        t=strfind(info.missing{k},'.');
        indices=resolve_address(info.missing{k}(1:t-1));
        if ~isempty(indices),
            found=0;
            if poi>0,
                for kk=1:poi,
                    if sum(abs(indices-incomplete(kk,:)))==0,
                        found=1;
                    end;
                end;
            end;
            if ~found,
                poi=poi+1;
                incomplete(poi,:)=indices;
            end;
        end;
    end;
    incomplete=incomplete(1:poi,:);
    add_msg_board(sprintf('For %i residues some atoms are missing.',poi));
end;
% copy secondary structure assignements per residue from chain model 1 to 
% other chain models, if any
for cnum=1:chains,
    nm = length(model.structures{snum}(cnum).residues);
    for km = 2:nm
        for kr = 1:length(model.structures{snum}(cnum).residues{1}.info)
            model.structures{snum}(cnum).residues{km}.info(kr).secondary = model.structures{snum}(cnum).residues{1}.info(kr).secondary;
        end
    end
end;
model.info{snum}.missing=info.missing;
model.info{snum}.incomplete=incomplete;
model.info{snum}.site_tags=info.site_tags;
model.info{snum}.sites=info.sites;
model.info{snum}.authors=info.authors;
model.info{snum}.keywords=info.keywords;
model.info{snum}.metal=info.metal;
model.info{snum}.resolution=info.resolution;
model.info{snum}.dssp=dssp;
model.info{snum}.rotamers=info.rotamers;
refs=0;
if isfield(model,'references'),
    refs=length(model.references);
end;
if refs==0,
    model.references=info.references;
elseif ~isempty(info.references),
    model.references(refs+1:refs+length(info.references))=info.references;
end;
% make site annotations, if any
if ~isempty(info.sites),
    preamble=['[' stag ']'];
    for k=1:length(info.sites),
        for kk=1:length(info.sites(k).residues),
            adr=sprintf('%s%s',preamble,info.sites(k).residues{kk});
            indices=resolve_address(adr);
            [m,n]=size(indices);
            if n==4, % valid residue address
                for kkk=1:m,
                    text=[info.sites(k).tag ': ' info.sites(k).declaration];
                    text=strvcat(text,sprintf('found in atomistic structure by %s',lower(info.sites(k).evidence)));
                    add_annotation(indices(kkk,:),'Binding',text,{'binding sites'});
                end;
            end;
        end;
    end;
end;
% make SSbond annotations, if any
if ~isempty(info.sites),
    preamble=['[' stag ']'];
    for k=1:length(info.SSbonds),
        adr1=sprintf('%s(%s)%i',preamble,info.SSbonds(k).chainID1,info.SSbonds(k).res1);
        indices1=resolve_address(adr1);
        adr2=sprintf('%s(%s)%i',preamble,info.SSbonds(k).chainID2,info.SSbonds(k).res2);
        indices2=resolve_address(adr2);
        [m1,n1]=size(indices1);
        [m2,n2]=size(indices2);
        if m1==1 && n1==4 && m2==1 && n2==4, % valid residue address
            text1=sprintf('S-S bond to %s with symmetry operator %s and length %5.2 Å',adr2,info.SSbonds(k).symop1,info.SSbonds(k).distance);
            text2=sprintf('S-S bond to %s with symmetry operator %s and length %5.2 Å',adr1,info.SSbonds(k).symop2,info.SSbonds(k).distance);
            add_annotation(indices1,'S-S bond',text1,{'S-S-bonds'});
            add_annotation(indices2,'S-S bond',text2,{'S-S-bonds'});
        end;
    end;
end;
% make missing atom annotations, if any
if ~isempty(model.info{snum}.missing),
    for k=1:length(model.info{snum}.missing),
        adr=model.info{snum}.missing{k};
        [res_adr,atag]=strtok(adr,'.');
        if ~isempty(res_adr),
            indices=resolve_address(res_adr);
            if ~isempty(indices),
                [m,n]=size(indices);
                for kk=1:m,
                    add_annotation(indices(kk,:),'Missing',atag(2:end),{'missing atoms'});
                end;
            end;
        end;
    end;
end;
% make metal coordination annotations, if any
if ~isempty(info.metal),
    for k=1:length(info.metal),
        center=strcat('[',stag,']',info.metal(k).center);
        indices=resolve_address(center);
        if ~isempty(indices) && length(info.metal(k).coordinated)>=1,
            found=1;
            text=strcat(stag,info.metal(k).coordinated{1});
            [m,n]=size(indices);
            if length(info.metal(k).coordinated)>=2,
                for kk=2:length(info.metal(k).coordinated),
                    text=strvcat(text,strcat(stag,info.metal(k).coordinated{kk}));
                end;
                for kk=1:m,
                    add_annotation(indices(kk,:),'Metal',text,{'metal centers'});
                end;
            end;
        else
            found=0;
        end;
        if length(info.metal(k).coordinated)>=1,
            for kk=1:length(info.metal(k).coordinated),
                adr=strcat('[',stag,']',info.metal(k).coordinated{kk});
                cindices=resolve_address(adr);
                if isempty(cindices),
                    add_msg_board(sprintf('Coordinating atom missing for %s',center));
                else
                    text=center;
                    if ~found,
                        text=strvcat(text,'but metal center not found in structure');
                    end;
                    [m,n]=size(cindices);
                    for kkk=1:m,
                        add_annotation(cindices(kkk,:),'Coordinating',text,{'metal centers'});
                    end;
                end;
            end;
        end;
    end;
end;
% make insertion annotations, if any
if ~isempty(info.insertions),
    [mins,nins]=size(info.insertions);
    for k=1:mins,
        indices=[snum info.insertions(k,:)];
        add_annotation(indices,'Inserted','Residue number not unique. Can be addressed only with insertion code.',{'inserted residues'});
    end;
end;
% make alternate location anotations, if any
if info.alternate,
    cnum=length(model.structures{snum});
    for kc=1:cnum,
        modnum=length(model.structures{snum}(kc).residues);
        for km=1:modnum,
            rnum=length(model.structures{snum}(kc).residues{km}.info);
            for kr=1:rnum, % residue loop
                if isfield(model.structures{snum}(kc).residues{km}.info(kr),'location_tags'),
                    if length(model.structures{snum}(kc).residues{km}.info(kr).location_tags)>3,
                        tags=model.structures{snum}(kc).residues{km}.info(kr).location_tags(2:end-1);
                        text='Location tags: ';
                        for ktc=1:length(tags),
                            if ~strcmp(tags(ktc),':'),
                                text=[text tags(ktc)];
                            else
                                text=[text ', '];
                            end;
                        end;
                        text=strtrim(text);
                        if strcmp(text(end),','),
                            text=strtrim(text(1:end-1));
                        end;
                        indices=[snum kc km kr];
                        add_annotation(indices,'Alternate',text,{'alternate locations'});
                    end;
                end;
            end;
        end;
    end;
end;

% make sequence modification (mutation) anotations, if any
cnum=length(model.structures{snum});
for kc=1:cnum,
    if isfield(model.structures{snum}(kc),'mutations') && ~isempty(model.structures{snum}(kc).mutations)
        for km=1:length(model.structures{snum}(kc).mutations),
            mutation=model.structures{snum}(kc).mutations(km);
            text=sprintf('Original residue: %s',mutation.original);
            if ~isempty(mutation.comment),
                text=sprintf('%s (%s)',text,mutation.comment);
            end;
            text=strtrim(text);
            modnum=length(model.structures{snum}(kc).residues);
            for kmodel=1:modnum,        
                indices=[snum kc kmodel mutation.number];
                add_annotation(indices,'Mutation',text,{'mutations'});
            end;
        end;
    end;
end;

set_ensemble_range(snum);

if isfield(model.info{snum},'resolution') && ~isempty(model.info{snum}.resolution),
    resstring=sprintf('%4.2f Å',model.info{snum}.resolution);
else
    resstring='not specified';
end;
    
set(hMain.MMM,'Name',sprintf('MMM - [%s](%s) Resolution %s',stag,model.current_chain,resstring));
set(hMain.menu_file_save_PDB,'Enable','on');