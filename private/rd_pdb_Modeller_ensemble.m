function [message,snum,stag,models]=rd_pdb_Modeller_ensemble(idCode,mod1,modn,pname)
% function [message,snum,models]=rd_pdb_Modeller_ensemble(idCode,mod1,modn,path)
%
% Reads several models in PDB format by Modeller into a single structure by
% calling rd_pdb, the PDB
% identifier, if present, is defined as a synonym for the structure
% if no model exists, a new model is created
%
% idCode        identifier code of Modeller target sequence, coincides with
%               first part of file name
% mod1          number of the first model to be read in
% modn          number of the last model to be read in
% pname         optional path identifier, defaults to current directory
% message       error message structure with fields
%               .error 0: no error, 1: PDB file does not exist, 2: PDB file
%               corrupted
% snum          number of the added structure in the model
% stag          actual structure tag (is idCode, if a structure by that
%               name did not yet exist)
% models        number of models that have been read in
%
% Remarks:
% - DSSP (if installed) is run only on the first model
%
% G. Jeschke, 2011

global model
global hMain
global general

message.error=0;
message.text='';
snum=[]; % empty structure number, if it fails
stag='';
models=0;

if nargin<4,
    pname=pwd;
end;

if modn<mod1,
    message.error=1;
    message.text='Number of last model smaller than number of first model.';
    return
end;
if modn<0 || mod1<0,
    message.error=2;
    message.text='Negative model number.';
    return
end;

if isempty(model) || ~isfield(model,'structures') || isempty(model.structures),
    snum=1; % this is the first structure
    model.structure_tags=':';
    model.structure_ids=[];
    model.unrecord=[];
    model.annotations=[];
else
    snum=length(model.structures)+1;
end;

failed=0;
for k=mod1:modn,
    modnum=sprintf('%i',k);
    while length(modnum)<4,
        modnum=['0' modnum];
    end;
    fname=sprintf('%s.B9999%s.pdb',idCode,modnum);
    cfname=fullfile(pname,fname);
    [cstructure,info]=rd_pdb(cfname,strcat(general.tmp_files,'remarks.txt'));
    if isempty(info) || isempty(cstructure) % PDB read was not successful
        if k==mod1,
            message.error=3;
            message.text='ERROR: PDB file of first model does not exist or is corrupted.';
            return;
        else
            if isfield(info,'no_file');
                message.error=4;
                message.text='Warning: Not all requested models exist.';
                add_msg_board(sprintf('Warning: PDB file of model %i does not exist.',k));
                failed=failed+1;
            else
                message.error=5;
                message.text='Warning: Some requested models are corrupted.';
                add_msg_board(sprintf('Warning: PDB file of model %i is corrupted.',k));
                failed=failed+1;
            end;
            if failed>3,
                message.error=6;
                message.text='Warning: Adding models was aborted after more than three PDB files could not be read.';
                return
            end;
        end;
    else % a new model is added
        set(hMain.menu_file_save_PDB,'Enable','on');
        models=models+1;
        if k==mod1, % the first model, structure is initialized
            info1=info;
            model.structures{snum}=cstructure;
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
            model.info{snum}.idCode=idCode;
            model.info{snum}.center=info.center;
            model.info{snum}.atoms=info.atoms;
            model.info{snum}.residues=info.residues;
            model.info{snum}.B_range=info.B_range;
            model.info{snum}.SSbonds=info.SSbonds;
            stag=model.info{snum}.idCode;
            id=tag2id(stag,model.structure_tags);
            poi=1;
            while ~isempty(id),
                stag=sprintf('%s_%i',model.info{snum}.idCode,poi);
                poi=poi+1;
                id=tag2id(stag,model.structure_tags);
            end;
            model.structure_tags=sprintf('%s%s:',model.structure_tags,stag);
            model.structure_ids=[model.structure_ids snum];
            model.current_structure=snum;
            model.current_chain=id2tag(1,model.chain_tags{snum});
            model.info{snum}.Modeller_obj=info.Modeller_obj;
            model.info{snum}.Modeller_sid=info.Modeller_sid;
            dssp=[];
            dospath=which('dssp.exe');
            if isempty(dospath),
                dospath=which('dsspcmbi.exe');
            end;            
            if ~isempty(dospath), 
                add_msg_board('DSSP geometry analysis is performed');
                infile=cfname;
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
                if s~=0,
                    add_msg_board('Warning: DSSP did not run successfully on this PDB file.');
                else
                    add_msg_board('Result file:');
                    add_msg_board(outfile);
                    dssp_file=outfile;
                    dssp=rd_dssp(dssp_file);
                    post_process_dssp(dssp,snum);
                end;
            end;        
        else % a new model is added to an existing ensemble
            model.info{snum}.Modeller_obj=[model.info{snum}.Modeller_obj info.Modeller_obj];
            for kc=1:length(cstructure),
                xyz0=model.structures{snum}(kc).xyz{1};
                xyz1=cstructure.xyz{1};
                [rmsd,xyz]=rmsd_superimpose(xyz0,xyz1);
                model.structures{snum}(kc).xyz{models}=xyz;
                model.structures{snum}(kc).Bfactor{models}=cstructure.Bfactor{1};
                model.structures{snum}(kc).Btensor{models}=cstructure.Btensor{1};
                model.structures{snum}(kc).atoms{models}=cstructure.atoms{1};
                model.structures{snum}(kc).residues{models}=cstructure.residues{1};
            end;
        end;
    end;    
end;

info=info1; % for consistency, all information records are from the first model


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
else
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