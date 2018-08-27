function analyze_CS_DEER(hObject, eventdata, handles)

targeted=true; % if true, frame with best fit of target (rather than constraints) is selected
converted=false; % if true, trajectory already converted to Calpha .mat file was read

global general

mypath=pwd;
cd(general.pdb_files);
[FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB template file');
if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Loading of PDB file canceled by user.');
    message.error=5;
    message.text='Cancelled.';
    return;
end;
reset_user_paths(PathName);
general.pdb_files=PathName;
initial=fullfile(PathName,FileName);
initial_PDB=FileName(1:4);

cd(general.restraint_files);

[FileName,PathName] = uigetfile({'*.dat'},'Select constraints file');
if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Loading of constraints canceled by user.');
    message.error=5;
    message.text='Cancelled.';
    return;
end;
reset_user_paths(PathName);
general.restraint_files=PathName;
restraints=fullfile(PathName,FileName);

cd(general.pdb_files);
[FileName,PathName] = uigetfile({'*.pdb;*.gz;*.mat'},'Select trajectory or ensemble file');
if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Loading of trajectory file canceled by user.');
    message.error=5;
    message.text='Cancelled.';
    return;
end;
trj=fullfile(PathName,FileName);
[~,~,ext]=fileparts(trj);
if strcmpi(ext,'.gz'),
    set(gcf,'Pointer','watch');
    drawnow;
    filenames = gunzip(trj);
    trj=fullfile(PathName,filenames{1});
    set(gcf,'Pointer','arrow');
end;
if strcmpi(ext,'.mat'),
    converted=true;
    trj2=fullfile(PathName,strcat(FileName(1:end-4),'.pdb'));
else
    trjout=trj(1:end-4);
    trj2=trj;
end;       
cd(mypath);

set(gcf,'Pointer','watch');
hf=gcf;
drawnow;

[message,initial_num]=add_pdb(initial);
template_coor = coarse_residues_and_store(initial_num);

restraint_data=rd_restraints(restraints);

[DEER,cancelled]=process_DEER_restraints(restraint_data);

if cancelled || isempty(DEER),
    add_msg_board('Processing of constraints failed.');
    message.error=1;
    message.text='Constraint file error.';
    set(gcf,'Pointer','arrow');
    return;
else
    fprintf(1,'%i DEER constraints are used.\n',length(DEER));
end;

if ~converted,
    [coor,residues,frames]=rd_pdb_trajectory_coarse(trj);
    save(trjout,'coor','residues','frames');
else
    load(trj);
end;
[m,n]=size(coor);

fprintf(1,'Number of residues: %i\n',residues);
fprintf(1,'Number of frames: %i\n',frames);
fprintf(1,'Size of coordinate array: [%i,%i]\n',m,n);
if m<residues*frames,
    frames=floor(m/residues);
    fprintf(2,'Warning. Incomplete trajectory. Truncating to %i frames.\n',frames);
end;
    
md=length(DEER);
crmsd=zeros(1,frames);
rmsd2template=zeros(1,frames);
for k=1:frames,
    frame=coor((k-1)*residues+1:k*residues,:);
    rmsd2template(k)=rmsd_superimpose(template_coor,frame);
    if md>0, % DEER constraints exist
        fom0=0;
        for kk=1:md,
            vec=frame(DEER(kk).net2,:)-frame(DEER(kk).net1,:);
            rp=norm(vec);
            fom0=fom0+(DEER(kk).r_CA-rp)^2;    
        end;
    end;
    crmsd(k)=sqrt(fom0/md);
end;
[mcrmsd,sel]=min(crmsd);

figure(1); clf;
plot(1:frames,rmsd2template,'b');
hold on;
plot(1:frames,crmsd,'r');
set(gca,'FontSize',14);
xlabel('Frame');
ylabel('rmsd (Å)');

fprintf(1,'Minimum rmsd to template     : %4.2f Å\n',min(rmsd2template));
fprintf(1,'Maximum rmsd to template     : %4.2f Å\n',max(rmsd2template));
fprintf(1,'Select model rmsd to template: %4.2f Å at frame %i\n',rmsd2template(sel),sel);
fprintf(1,'Minimum constraint rmsd      : %4.2f Å\n',mcrmsd);

% 
% oname=sprintf('%s_to_%s_%s.pdb',initial_PDB,target_PDB,runtype);
% tic,
% extract_frame(trj2,oname,sel);
% toc,

set(hf,'Pointer','arrow');

function correspondence=align_template_target(snum_template,snum_target)
% returns Calpha coordinates for a target structure and a correspondence
% table between Ca coordinate arrays of template and target structure based
% on sequence alignment with MUSCLE

global model

correspondence=[];

rindices1=model.coarse(snum_template).indices;
rindices2=model.coarse(snum_target).indices;

cid1=model.chain_ids(snum_template);
cid1=cid1{1};
cid2=model.chain_ids(snum_target);
cid2=cid2{1};
if length(cid1)~=length(cid2),
    add_msg_board('ERROR: Different number of chains in template and moving structure.');
    add_msg_board('Deactivate "whole structure".');
    return
end;
sel1=zeros(2000,4);
sel2=sel1;
psel=0;
for k=1:length(cid1),
    seqs{1}=model.structures{snum_template}(cid1(k)).sequence;
    seqs{2}=model.structures{snum_target}(cid2(k)).sequence;
    sindices=[snum_template,cid1(k);snum_target,cid2(k)];
    [message,inname]=align_sequences(seqs,sindices,true);
    if message.error,
        add_msg_board('ERROR: MUSCLE sequence alignment failed.');
        add_msg_board(message.text);
        add_msg_board('Deactivate "whole structure".');
        return
    end;
    alignment=get_multiple_clustal(inname);
    [asel1,asel2]=select_aligned(alignment,sindices);
    [msel,nsel]=size(asel1);
    sel1(psel+1:psel+msel,:)=asel1;
    sel2(psel+1:psel+msel,:)=asel2;
    psel=psel+msel;
end;
sel1=sel1(1:psel,:);
sel2=sel2(1:psel,:);

[mc1,nc1]=size(rindices1);
[mc2,nc2]=size(rindices2);
correspondence=zeros(2,mc1);
poi=0;
for k=1:psel,
    cindices=sel1(k,:);
    diff=repmat(cindices,mc1,1)-rindices1; % check if residue is in Calpha coordinate array of template structure
    [mdiff,ndiff]=size(diff);
    diff(:,3)=zeros(mdiff,1); % mask chain model number
    adiff=sum(abs(diff),2);
    cpoi1=find(adiff==0);
    if ~isempty(cpoi1), % if index into Calpha coordinate array was found for template
        cindices=sel2(k,:);
        diff=repmat(cindices,mc2,1)-rindices2; % check if residue is in Calpha coordinate array of target structure
        [mdiff,ndiff]=size(diff);
        diff(:,3)=zeros(mdiff,1); % mask chain model number
        adiff=sum(abs(diff),2);
        cpoi2=find(adiff==0);
        if ~isempty(cpoi2), % corresponding residues exist in both Calpha coordinate arrays
            poi=poi+1;
            correspondence(1,poi)=cpoi1;
            correspondence(2,poi)=cpoi2;
        end;
    end;
end;

correspondence=correspondence(:,1:poi);

function Ca_coor = coarse_residues_and_store(snum,modnum)

global model

if nargin<2,
    modnum=[];
end;

[Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(snum,modnum);
model.coarse(snum).Ca_coor=Ca_coor;
model.coarse(snum).indices=rindices;
model.coarse(snum).masses=masses;
model.coarse(snum).Bfactors=Bfactors;
model.coarse(snum).restypes=restypes;

function [sel1,sel2]=select_aligned(alignment,sindices)
% select only residues that are aligned

global model

seq1=alignment(1).sequence;
seq2=alignment(2).sequence;
sel1=zeros(length(seq1),4);
sel2=sel1;
poi1=0;
poi2=0;
pois=0;
for k=1:length(seq1),
    match1=false;
    match2=false;
    if char(seq1(k))~='-',
        poi1=poi1+1;
        match1=true;
    end;
    if char(seq2(k))~='-',
        poi2=poi2+1;
        match2=true;
    end;
    if match1 && match2
        tag1=sprintf('%i',poi1);
        id1=tag2id(tag1,model.structures{sindices(1,1)}(sindices(1,2)).residues{1}.residue_tags);
        tag2=sprintf('%i',poi2);
        id2=tag2id(tag2,model.structures{sindices(2,1)}(sindices(2,2)).residues{1}.residue_tags);
        if ~isempty(id1) && ~isempty(id2),
            pois=pois+1;
            sel1(pois,:)=[sindices(1,:) 1 id1];
            sel2(pois,:)=[sindices(2,:) 1 id2];
        end;
    end;
end;
sel1=sel1(1:pois,:);
sel2=sel2(1:pois,:);


function extract_frame(trj,oname,sel)

global general

mypath=pwd;
cd(general.pdb_files);

fid=fopen(trj);
if fid==-1,
    return;
end;

ofile=fopen(oname,'wt');
if fid==-1,
    return;
end;
fprintf(ofile,'REMARK Frame for structural transition extracted by MMM\n');

nl=0;
frame=0;
wmode=false;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    while length(tline)<6, % catches too short end line of NMsim files
        tline=[tline ' '];
    end;
    record=tline(1:6);
    switch record
        case 'MODEL '
            frame=frame+1;
            if frame==sel, wmode=true; end;
        case 'ENDMDL'
            wmode=false;
    end;
    if wmode  && ~strcmpi(record,'MODEL '),
        fprintf(ofile,'%s\n',tline);
    end;
end;

fprintf(ofile,'END    \n');
fclose(ofile);
fclose(fid);

cd(mypath);

function [DEER,cancelled]=process_DEER_restraints(restraints)

global model
global hMain

cancelled=false;

if ~isfield(restraints,'DEER'),
    DEER=[];
    return;
end;

snum=model.current_structure;

if isfield(restraints,'PDB'),
    if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
        button = questdlg(sprintf('Constraint file specifies template %s, while current template is %s. Do you want to continue?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
        if strcmp(button,'No'),
            cancelled=true;
            DEER=[];
            return
        end;
    end;
end;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

% check whether sites are already labelled and whether all restraint sites
% do exist
% identity of the label is checked
% labeling temperature is NOT checked
T_list=zeros(1,200);
if ~isempty(labels),
    lindices=zeros(length(labels),4);
    for k=1:length(labels),
        cindices=labels(k).indices;
        if ~isempty(cindices),
            lindices(k,:)=cindices;
        end;
    end;
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        ind1=resolve_address(adr1);
        if isempty(ind1),
            add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            add_msg_board('Processing of DEER constraints cancelled');
            cancelled=true;
            DEER=[];
            return;
        end;
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,adr1));
            end;
        end;
        adr2=restraints.DEER(k).adr2;
        ind2=resolve_address(adr2);
        if isempty(ind2),
            add_msg_board(sprintf('ERROR: Constraint %i has second label at site %s',k,adr2));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            add_msg_board('Processing of DEER constraints cancelled');
            cancelled=true;
            DEER=[];
            return;
        end;
        found=false;
        for l=1:length(labels),
            diff=ind2-lindices(l,:);
            if sum(abs(diff))==0  && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr2;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,adr2));
            end;
        end;
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            label_list{poi}=restraints.DEER(k).label;
            T_list(poi)=restraints.DEER(k).T;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
        adr2=restraints.DEER(k).adr2;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr2;
            label_list{poi}=restraints.DEER(k).label;
            T_list(poi)=restraints.DEER(k).T;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr2));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

labels=label_information(model.sites);
cindices=model.coarse(snum).indices;
[mn,nn]=size(cindices);

for k=1:length(restraints.DEER),
    cindices1=resolve_address(sprintf('%s.CA',restraints.DEER(k).adr1));
    if isempty(cindices1) || length(cindices1)~=5,
        add_msg_board(sprintf('Warning: DEER constraint site %s does not exist in template. Constraint %i will be ignored.',constraints.direct(k).adr1,k));
        continue
    else
       [message,xyz]=get_atom(cindices1,'coor');
       if message.error,
            add_msg_board(sprintf('Warning: Failed to find Calpha coordinates for DEER constraint site %s (Errror: %s). Constraint %i will be ignored.',constraints.direct(k).adr1,message.text,k));
       else
           restraints.DEER(k).Ca_xyz1=xyz;
       end;
    end;
    cindices2=resolve_address(sprintf('%s.CA',restraints.DEER(k).adr2));
    if isempty(cindices1) || length(cindices2)~=5,
        add_msg_board(sprintf('Warning: DEER constraint site %s does not exist in template. Constraint %i will be ignored.',constraints.direct(k).adr2,k));
        continue
    else
       [message,xyz]=get_atom(cindices2,'coor');
       if message.error,
            add_msg_board(sprintf('Warning: Failed to find Calpha coordinates for DEER constraint site %s (Errror: %s). Constraint %i will be ignored.',constraints.direct(k).adr2,message.text,k));
       else
           restraints.DEER(k).Ca_xyz2=xyz;
       end;
    end;
    adr1=restraints.DEER(k).adr1;
    ind1=resolve_address(adr1);
    adr2=restraints.DEER(k).adr2;
    ind2=resolve_address(adr2);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Specified residue %s does not exist in template structure.',adr1));
        cancelled=true;
        DEER=[];
        return;
    end;
    net1=0;
    for kk=1:mn,
        match=sum(abs(ind1-cindices(kk,:)));
        if match==0,
           net1=kk;
           break;
        end;
    end;
    restraints.DEER(k).net1=net1;
    net2=0;
    for kk=1:mn,
        match=sum(abs(ind2-cindices(kk,:)));
        if match==0,
           net2=kk;
           break;
        end;
    end;
    restraints.DEER(k).net2=net2;
    if net1==0 || net2==0,
        add_msg_board('ERROR: Constrained C_alpha-C_alpha distance does not exist in network model.');
        cancelled=true;
        DEER=[];
        return;
    end;
    DEER(k).r=restraints.DEER(k).r;
    DEER(k).sigr=restraints.DEER(k).sigr;
    DEER(k).indices=[ind1;ind2];
    DEER(k).adr1=adr1;
    DEER(k).adr2=adr2;
    DEER(k).Ca_xyz1=restraints.DEER(k).Ca_xyz1;
    DEER(k).Ca_xyz2=restraints.DEER(k).Ca_xyz2;
    DEER(k).net1=restraints.DEER(k).net1;
    DEER(k).net2=restraints.DEER(k).net2;
    f1=false;
    f2=false;
    for l=1:length(labels),
        diff1=ind1-labels(l).indices;
        if sum(abs(diff1))==0,
            f1=true;
            DEER(k).xyz1=labels(l).xyz;
            DEER(k).rmsd1=labels(l).rmsd;
        end;
        diff2=ind2-labels(l).indices;
        if sum(abs(diff2))==0,
            f2=true;
            DEER(k).xyz2=labels(l).xyz;
            DEER(k).rmsd2=labels(l).rmsd;
        end;
    end;
    % translate N-O midpoint distance to Calpha-Calpha distance constraint
    r_NO=norm(DEER(k).xyz2-DEER(k).xyz1);
    r_DEER=10*DEER(k).r;
    q=r_NO^2-r_DEER^2;
    r_CA=norm(DEER(k).Ca_xyz2-DEER(k).Ca_xyz1);
    p=2*r_CA;
    dr=-p/2+sqrt(p^2/4-q);
    DEER(k).r_CA=r_CA+dr;
    if ~f1 || ~f2,
        add_msg_board('ERROR: Automatic rotamer computation error.');
        add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
        cancelled=true;
        DEER=[];
        return;
    end;
end;

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
        end;
    end;
end;

function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm
