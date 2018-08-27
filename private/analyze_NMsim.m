function analyze_NMsim(initial,target,restraints,trj,runtype)

targeted=true; % if true, frame with best fit of target (rather than constraints) is selected
converted=false; % if true, trajectory already converted to Calpha .mat file was read

global general
if nargin<5,
    runtype='NMsim';
end;

if nargin<1 || isempty(initial),
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
    cd(mypath);
end;

if nargin<2 || isempty(target),
    mypath=pwd;
    cd(general.pdb_files);
    [FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB target file');
    if isequal(FileName,0) || isequal(PathName,0),
        add_msg_board('Loading of PDB file canceled by user.');
        message.error=5;
        message.text='Cancelled.';
        return;
    end;
    reset_user_paths(PathName);
    general.pdb_files=PathName;
    target=fullfile(PathName,FileName);
    target_PDB=FileName(1:4);
    cd(mypath);
end;

my_path=pwd;
cd(general.restraint_files);

if nargin<3 || isempty(restraints),
    mypath=pwd;
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
    cd(mypath);
end;

if nargin<4 || isempty(trj),
    mypath=pwd;
    cd(general.pdb_files);
    [FileName,PathName] = uigetfile({'*.pdb;*.gz;*.mat'},'Select NMsim trajectory file');
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
end;

set(gcf,'Pointer','watch');
drawnow;
[message,target_num]=add_pdb(target);

[message,initial_num]=add_pdb(initial);

template_coor = coarse_residues_and_store(initial_num);
target_coor = coarse_residues_and_store(target_num);

correspondence=align_template_target(initial_num,target_num);
template_coor_m=template_coor(correspondence(1,:),:);
target_coor_m=target_coor(correspondence(2,:),:);

rmsd0=rmsd_superimpose(template_coor_m,target_coor_m);

fprintf(1,'Template (%s)/target (%s) Calpha rmsd: %4.2f Å\n',initial_PDB,target_PDB,rmsd0);

restraint_data=rd_restraints(restraints);

[direct,cancelled]=process_direct_restraints(restraint_data,initial_num);

if cancelled || isempty(direct),
    add_msg_board('Processing of constraints failed.');
    message.error=1;
    message.text='Constraint file erro.';
    return;
end;

if ~converted,
    tic,
    [coor,residues,frames]=rd_pdb_trajectory_coarse(trj);
    save(trjout,'coor','residues','frames');
    toc,
else
    load(trj);
end;
[m,n]=size(coor);

fprintf(1,'Number of residues: %i\n',residues);
fprintf(1,'Number of frames: %i\n',frames);
fprintf(1,'Size of coordinate array: [%i,%i]\n',m,n);

rmsd2template=zeros(1,frames);
rmsd2target=zeros(1,frames);
crmsd=zeros(1,frames);
brmsd=zeros(1,frames);
capproach=1e6;

[mc,nc]=size(direct);
tic,
for k=1:frames,
    frame=coor((k-1)*residues+1:k*residues,:);
    rmsd2template(k)=rmsd_superimpose(template_coor,frame);
    frame_m=frame(correspondence(1,:),:);
    rmsd2target(k)=rmsd_superimpose(target_coor_m,frame_m);
    if rmsd2target(k)<capproach,
        capproach=rmsd2target(k);
    end;
    brmsd(k)=capproach;
    if mc>0, % direct constraints exist
        fom0=0;
        for kk=1:mc,
            vec=frame(direct(kk,2),:)-frame(direct(kk,1),:);
            rp=norm(vec);
            fom0=fom0+(10*direct(kk,4)-rp)^2;    
        end;
    end;
    crmsd(k)=sqrt(fom0/mc);
end;
toc,
[mrmsdt,best]=min(rmsd2target);
[mcrmsd,sel]=min(crmsd);

figure(2); clf;
plot(1:frames,brmsd,'b');
set(gca,'FontSize',14);
xlabel('Frame');
ylabel('min. rmsd to target (Å)');

figure(1); clf;
plot(1:frames,rmsd2template,'b');
hold on;
plot(1:frames,rmsd2target,'r');
plot(1:frames,crmsd,'g');
plot([1,frames],[rmsd0,rmsd0],'k:');
set(gca,'FontSize',14);
xlabel('Frame');
ylabel('rmsd (Å)');

cc=corrcoef(crmsd,rmsd2target);

fprintf(1,'Minimum rmsd to template: %4.2f Å\n',min(rmsd2template));
fprintf(1,'Maximum rmsd to template: %4.2f Å\n',max(rmsd2template));
fprintf(1,'Minimum rmsd to target  : %4.2f Å\n',min(rmsd2target));
fprintf(1,'Maximum rmsd to target  : %4.2f Å\n',max(rmsd2target));
fprintf(1,'Best coverage           : %4.3f at frame %i\n',(rmsd0-mrmsdt)/rmsd0,best);
fprintf(1,'Select model rmsd       : %4.2f Å at frame %i\n',rmsd2target(sel),sel);
fprintf(1,'Minimum constraint rmsd : %4.2f Å\n',mcrmsd);
fprintf(1,'Selected frame coverage : %4.3f\n',(rmsd0-rmsd2target(sel))/rmsd0);
fprintf(1,'Corr.coeff constraint/target rmsd: %6.4f\n',cc(1,2));

if targeted,
    sel=best;
end;

oname=sprintf('%s_to_%s_%s.pdb',initial_PDB,target_PDB,runtype);
tic,
extract_frame(trj2,oname,sel);
toc,

set(gcf,'Pointer','arrow');

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

function [direct,cancelled]=process_direct_restraints(restraints,snum)

global model

model.current_structure=snum;

cancelled=false;
if ~isfield(restraints,'direct'),
    direct=[];
    return;
end;

md=length(restraints.direct);
direct=zeros(md,5);

cindices=model.coarse(snum).indices;
[mn,nn]=size(cindices);
for k=1:md,
    adr1=restraints.direct(k).adr1;
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Specified residue %s does not exist in template structure.',adr1));
        cancelled=true;
        direct=[];
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
    direct(k,1)=net1;
    adr2=restraints.direct(k).adr2;
    ind2=resolve_address(adr2);
    net2=0;
    for kk=1:mn,
        match=sum(abs(ind2-cindices(kk,:)));
        if match==0,
           net2=kk;
           break;
        end;
    end;
    direct(k,2)=net2;    
    direct(k,4)=restraints.direct(k).r;
    direct(k,5)=restraints.direct(k).sigr;
    if net1==0 || net2==0,
        add_msg_board('ERROR: Constrained C_alpha-C_alpha distance does not exist in network model.');
        cancelled=true;
        direct=[];
        return;
    end;
end;

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