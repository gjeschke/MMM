function [ merit, f_ENM, f_hom, message ] = analyze_transition_fit( initial, final, ENM_fit, hom_fit )
%ANALYZE_TRANSITION_FIT Determines fit qualities and agreement between
%different fits for a structural transition
%   initial     pathname of PDB file for initial state
%   final       pathname of PDB file for final state
%   ENM_fit     pathname of PDB file for final state fitted by elastic
%               network model
%   hom_fit     pathname of PDB file for final state fitted by homology
%               modeling
%   merit       number of merit for the agreement
%   f_ENM       fractional coverage of ENM fit
%   f_hom       fractional coverage of homology model
%   message     error message with fields
%               .error: 0 no error, otherwise error code
%               .text: 'Initial structure missing.', error code 1
%                      'Final structure missing.', error code 2
%                      'ENM fit missing.', error code 3
%                      'Homology fit missing.', error code 4
%                      'Cancelled.', error code 5
%                      'Residue number mismatch.', error code 6
%
% G. Jeschke, 3.9.2012

global general

ENM_weight=0.5;

automatic = true;

initial_PDB='';
final_PDB='';
fit_path='E:/projects/homology_modeling/ensembles';
log_path='E:/projects/homology_modeling/ensembles/Modeller_ensemble_50_logfiles';
log_suffix='ensemble_50.log';
fit_type='50_label_ensemble';

merit = [];
f_ENM = [];
f_hom = [];
message.error=0;
message.text='No error.';

if nargin<1,
    [FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file for initial state');
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
end;



if nargin<2,
    [FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file for final state');
    if isequal(FileName,0) || isequal(PathName,0),
        add_msg_board('Loading of PDB file canceled by user.');
        message.error=5;
        message.text='Cancelled.';
        return;
    end;
    reset_user_paths(PathName);
    general.pdb_files=PathName;
    final=fullfile(PathName,FileName);
    final_PDB=FileName(1:4);
end;

ENM_file=sprintf('%s_to_%s_%s_edENM.pdb',initial_PDB,final_PDB,fit_type);
hom_file=sprintf('%s_to_%s_%s_homology.pdb',initial_PDB,final_PDB,fit_type);
log_file=sprintf('%s_%s',final_PDB,log_suffix);
log_ENM=fullfile(log_path,log_file);

ENM_representative=scan_log(log_ENM);

if isempty(ENM_representative),
    return;
else
    fprintf(1,'Representative ENM structure is #%i.\n',ENM_representative);
end;


if nargin<3 
    if automatic,
        ENM_fit = fullfile(fit_path,ENM_file);
    else
        [FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file for ENM fit');
        if isequal(FileName,0) || isequal(PathName,0),
            add_msg_board('Loading of PDB file canceled by user.');
            message.error=5;
            message.text='Cancelled.';
            return;
        end;
        reset_user_paths(PathName);
        general.pdb_files=PathName;
        ENM_fit=fullfile(PathName,FileName);
    end;
end;

if nargin<4,
    if automatic,
        hom_fit = fullfile(fit_path,hom_file);
    else
        [FileName,PathName,FilterIndex] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file for final state');
        if isequal(FileName,0) || isequal(PathName,0),
            add_msg_board('Loading of PDB file canceled by user.');
            message.error=5;
            message.text='Cancelled.';
            return;
        end;
        reset_user_paths(PathName);
        general.pdb_files=PathName;
        hom_fit=fullfile(PathName,FileName);
    end;
end;

[message,initial_num]=add_pdb(initial);
if message.error,
    message.error=1;
    message.text='Initial structure missing.';
    return
end;
[message,final_num]=add_pdb(final);
if message.error,
    message.error=2;
    message.text='Final structure missing.';
    return
end;
disp(ENM_fit);
[message,ENM_num]=add_pdb(ENM_fit);
if message.error,
    message.error=3;
    message.text='ENM fit missing.';
    return
end;
[message,hom_num]=add_pdb(hom_fit);
if message.error,
    message.error=4;
    message.text='Homology fit missing.';
    return
end;

initial_coor = coarse_residues_and_store(initial_num);
if isempty(initial_coor),
    add_msg_board('Initial structure could not be coarse grained');
    message.error=1;
    message.text='Initial structure missing.';
    return
end

final_coor = coarse_residues_and_store(final_num);
if isempty(final_coor),
    add_msg_board('Final structure could not be coarse grained');
    message.error=2;
    message.text='Final structure missing.';
    return
end

ENM_coor = coarse_residues_and_store(ENM_num,ENM_representative);
if isempty(ENM_coor),
    add_msg_board('ENM fit could not be coarse grained');
    message.error=3;
    message.text='ENM fit missing.';
    return
end

hom_coor = coarse_residues_and_store(hom_num);
if isempty(hom_coor),
    add_msg_board('Homology fit could not be coarse grained');
    message.error=4;
    message.text='Homology fit missing.';
    return
end

isize=size(initial_coor);
fsize=size(final_coor);
esize=size(ENM_coor);
hsize=size(hom_coor);

if isize(2)~=3,
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Calpha coordinates of initial structure are improper.');
    return
end;

if fsize(2)~=3,
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Calpha coordinates of final structure are improper.');
    return
end;

if esize(2)~=3,
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Calpha coordinates of ENM fit are improper.');
    return
end;

if hsize(2)~=3,
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Calpha coordinates of homology fit are improper.');
    return
end;

correspondence=align_template_target(initial_num,final_num);

if isempty(correspondence),
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Final structure could not be aligned to initial structure.');
    return
else
    final_coor_0=final_coor;
    correspondence_E=align_template_target(initial_num,ENM_num);
    initial_coor_E=initial_coor(correspondence_E(1,:),:);
    ENM_coor_i=ENM_coor(correspondence_E(2,:),:);
    correspondence_h=align_template_target(initial_num,hom_num);
    initial_coor_h=initial_coor(correspondence_h(1,:),:);
    hom_coor_i=hom_coor(correspondence_h(2,:),:);
    correspondence_E=align_template_target(final_num,ENM_num);
    final_coor_E=final_coor(correspondence_E(1,:),:);
    ENM_coor_f=ENM_coor(correspondence_E(2,:),:);
    correspondence_h=align_template_target(final_num,hom_num);
    final_coor_h=final_coor(correspondence_h(1,:),:);
    initial_coor=initial_coor(correspondence(1,:),:);
    final_coor=final_coor(correspondence(2,:),:);
    hom_coor_f=hom_coor(correspondence_h(2,:),:);
    correspondence_eh=align_template_target(ENM_num,hom_num);
    ENM_coor_s=ENM_coor(correspondence_eh(1,:),:);
    hom_coor_s=hom_coor(correspondence_eh(2,:),:);
    mixed_coor=ENM_weight*ENM_coor_s+(1-ENM_weight)*hom_coor_s;
    correspondence_mh=correspondence_eh;
    [~,mm]=size(correspondence_eh);
    poi=0;
    for k=1:mm,
        indi=correspondence_eh(2,k);
        ipoi = find(correspondence_h(2,:)==indi,1);
        if ~isempty(ipoi),
            poi=poi+1;
            correspondence_mh(1,poi)=correspondence_h(1,ipoi);
            correspondence_mh(2,poi)=k;
        end;
    end;
    correspondence_mh=correspondence_mh(:,1:poi);
    final_coor_mh=final_coor_0(correspondence_mh(1,:),:);
    mixed_coor_mh=mixed_coor(correspondence_mh(2,:),:);
    correspondence_mE=correspondence_E;
    [~,mm]=size(correspondence_mE);
    poi=0;
    for k=1:mm,
        indi=correspondence_eh(1,k);
        ipoi = find(correspondence_E(2,:)==indi,1);
        if ~isempty(ipoi),
            poi=poi+1;
            correspondence_mE(1,poi)=correspondence_E(1,ipoi);
            correspondence_mE(2,poi)=k;
        end;
    end;
    correspondence_mE=correspondence_mE(:,1:poi);
    final_coor_mE=final_coor_0(correspondence_mE(1,:),:);
    mixed_coor_mE=mixed_coor(correspondence_mE(2,:),:);
    isize=size(initial_coor);
    isize_E=size(initial_coor_E);
    isize_h=size(initial_coor_h);
    fsize=size(final_coor);
    fsize_E=size(final_coor_E);
    fsize_h=size(final_coor_h);
    esize=size(ENM_coor);
    hsize=size(hom_coor);
    esize_i=size(ENM_coor_i);
    hsize_i=size(hom_coor_i);
    esize_f=size(ENM_coor_f);
    hsize_f=size(hom_coor_f);
    if isize(1)~=fsize(1),
        message.error = 6;
        message.text = 'Residue number mismatch.';
        add_msg_board('Residue number of final structure does not match residue number of initial structure.');
        return
    end;
end;

if esize_i(1)~=isize_E(1),
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Residue number of ENM fit does not match residue number of initial structure.');
    return
end;

if hsize_i(1)~=isize_h(1),
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Residue number of homology fit does not match residue number of initial structure.');
    return
end;

% at this point we have four coordinate arrays of the same size with
% aligned residues

rmsd0=rmsd_superimpose(final_coor,initial_coor);

add_msg_board(sprintf('Calpha r.m.s.d. between initial and final state is %6.2f Å',rmsd0));

rmsd_ENM = rmsd_superimpose(final_coor_E,ENM_coor_f);
change_ENM = rmsd_superimpose(initial_coor_E,ENM_coor_i);

f_ENM=(rmsd0-rmsd_ENM)/rmsd0;

[mm,~]=size(final_coor_E);
add_msg_board(sprintf('Fractional coverage of ENM fit is %6.3f for %i residues.',f_ENM,mm));
add_msg_board(sprintf('Calpha r.m.s.d. between ENM fit and initial state is %6.2f Å',change_ENM));

rmsd_hom = rmsd_superimpose(final_coor_h,hom_coor_f);
change_hom = rmsd_superimpose(initial_coor_h,hom_coor_i);

f_hom=(rmsd0-rmsd_hom)/rmsd0;

[mm,~]=size(final_coor_h);
add_msg_board(sprintf('Fractional coverage of homology fit is %6.3f for %i residues.',f_hom,mm));
add_msg_board(sprintf('Calpha r.m.s.d. between homology fit and initial state is %6.2f Å',change_hom));

[mm,~]=size(ENM_coor_s);
rmsd_fits = rmsd_superimpose(ENM_coor_s,hom_coor_s);
add_msg_board(sprintf('Calpha r.m.s.d. between ENM fit and homology fit is %6.2f Å for %i residues',rmsd_fits,mm));

ra=change_ENM;
rb=change_hom;
rc=rmsd_fits;
s=(ra+rb+rc)/2;
% thbc0=2*asin(sqrt((s-rb)*(s-rc)/(rb*rc))); % Halbwinkelsätze
% thac0=2*asin(sqrt((s-ra)*(s-rc)/(ra*rc)));
thab0=2*asin(sqrt((s-ra)*(s-rb)/(ra*rb)));

merit=cos(thab0);
add_msg_board(sprintf('Opening angle: %4.1f, Figure of merit: %5.3f',180*thab0/pi,merit));

x0=0;
y0=0;
x1=0;
y1=ra;
x2=sin(thab0)*rb;
y2=cos(thab0)*rb;

rmsd0=rmsd_superimpose(final_coor,initial_coor);

[mm,~]=size(final_coor_mh);
rmsd_fits = rmsd_superimpose(final_coor_mh,mixed_coor_mh);
add_msg_board(sprintf('Calpha r.m.s.d. between average fit and final structure is %6.2f Å for % i residues.',rmsd_fits,mm));
f_mixed=(rmsd0-rmsd_fits)/rmsd0;
add_msg_board(sprintf('Fractional coverage of average fit is %5.3f',f_mixed));

rmsd_fits = rmsd_superimpose(final_coor_mE,mixed_coor_mE);
add_msg_board(sprintf('Calpha r.m.s.d. between average fit and final structure is %6.2f Å',rmsd_fits));
f_mixed=(rmsd0-rmsd_fits)/rmsd0;
add_msg_board(sprintf('Fractional coverage of average fit is %5.3f',f_mixed));

fprintf(1,'\n');
for kk=0:20,
    ENM_weight=kk/20;
    mixed_coor=ENM_weight*ENM_coor_s+(1-ENM_weight)*hom_coor_s;
    mixed_coor_mh=mixed_coor(correspondence_mh(2,:),:);
    rmsd_fits = rmsd_superimpose(final_coor_mh,mixed_coor_mh);
    f_mixed=(rmsd0-rmsd_fits)/rmsd0;
    fprintf(1,'%5.3f ',f_mixed);
end;
fprintf(1,'\n');

figure(1); clf;
plot(x0,y0,'k.');
hold on;
plot(x1,y1,'b.');
plot(x2,y2,'r.');
plot([x0,x1],[y0,y1],'k');
plot([x0,x2],[y0,y2],'k');
plot([x1,x2],[y1,y2],'k');
axis([-0.1,max([x0,x1,x2])+0.1,-0.1,max([y0,y1,y2])+0.1]);
axis equal

% Distance difference matrix analysis

mymap=ones(65,3);
for k=1:32;
    nsat=(k-1)/32;
    mymap(k,1)=nsat;
    mymap(k,2)=nsat;
    mymap(66-k,2)=nsat;
    mymap(66-k,3)=nsat;
end;
simplified=false;

ddm_x=get_ddm(initial_coor,final_coor,simplified);

figure(2); clf;
title('x-ray difference distance matrix');
if simplified,
    colormap(gray);
    mymap=colormap;
else
    colormap(mymap);
end;
[m,n]=size(mymap);
image(m*ddm_x);
axis equal
axis off
axis xy

ddm_E=get_ddm(initial_coor_E,ENM_coor_i,simplified);

figure(3); clf;
title('ENM difference distance matrix');
if simplified,
    colormap(gray);
    mymap=colormap;
else
    colormap(mymap);
end;
[m,n]=size(mymap);
image(m*ddm_E);
axis equal
axis off
axis xy

ddm_h=get_ddm(initial_coor_h,hom_coor_i,simplified);

figure(4); clf;
title('Modeller difference distance matrix');
if simplified,
    colormap(gray);
    mymap=colormap;
else
    colormap(mymap);
end;
mymap=colormap;
[m,n]=size(mymap);
image(m*ddm_h);
axis equal
axis xy
axis off

ddm_E=2*ddm_E-1;
ddm_h=2*ddm_h-1;
entropy_E=matrix_entropy(ddm_E);
entropy_h=matrix_entropy(ddm_h);
add_msg_board(sprintf('Entropy of ENM DDM     : %6.4f.',entropy_E));
add_msg_board(sprintf('Entropy of Modeller DDM: %6.4f.',entropy_h));

[mE,~]=size(ddm_E);
[mh,~]=size(ddm_h);

if mE==mh,
    ddm_consensus=zeros(mE);
    for k=1:mE-1,
        for kk=k+1:mE,
            if sign(ddm_E(k,kk))==sign(ddm_h(k,kk)),
                elm=(ddm_E(k,kk)+ddm_h(k,kk))/2;
            else
                elm=0;
            end;
            ddm_consensus(k,kk)=elm;
            ddm_consensus(kk,k)=elm;
        end;
    end;
    ddm_consensus=ddm_consensus/max(max(abs(ddm_consensus)));
%     threshold=1/5;
%     ddm_consensus(abs(ddm_consensus)<threshold)=0;
    merit_E=sum(sum(ddm_consensus.*ddm_E))/mE^2;
    merit_h=sum(sum(ddm_consensus.*ddm_h))/mE^2;
    dissimilarity=matrix_dissimilarity(ddm_E,ddm_h);
    add_msg_board(sprintf('Dissimilarity of ENM and Modeller DDM: %6.4f.',dissimilarity));
    add_msg_board(sprintf('Agreement of consensus DDM with ENM DDM     : %6.4f.',merit_E));
    add_msg_board(sprintf('Agreement of consensus DDM with Modeller DDM: %6.4f.',merit_h));
    
    ddm_consensus=(ddm_consensus+1)/2;

    figure(5); clf;
    title('Consensus difference distance matrix');
    if simplified,
        colormap(gray);
        mymap=colormap;
    else
        colormap(mymap);
    end;
    [m,n]=size(mymap);
    image(m*ddm_consensus);
    axis equal
    axis off
    axis xy
else
    figure(5);
    close(5);
end;


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

function ENM_representative=scan_log(log_ENM)

ENM_representative=[];

key='The most typical result was obtained in trial';

fid=fopen(log_ENM);
if fid==-1,
    fprintf(2,'%s\n',log_ENM);
    fprintf(2,'ERROR: ENM logfile could not be opened.\n');
    return;
end;

nl=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    poi=strfind(tline,key);
    if ~isempty(poi),
        rem=tline(poi+length(key):end);
        ENM_representative=str2double(rem);
    end;
end;

function ddm_x=get_ddm(coor1,coor2,simplified)

filter=true;

xi_dmat=coor2dmat(coor1);
xf_dmat=coor2dmat(coor2);
ddm_x=xi_dmat-xf_dmat;
if filter,
    threshold=max(max(abs(ddm_x)))/5;
    [m,n]=size(ddm_x);
    for k=1:m-1,
        for kk=k+1:m,
            if abs(ddm_x(k,kk))<threshold,
                ddm_x(k,kk)=0;
                ddm_x(kk,k)=0;
            end;
        end;
    end;
end;
if simplified,
    [m,n]=size(ddm_x);
    for k=1:m-1,
        for kk=k+1:m,
            if abs(ddm_x(k,kk))<0.5,
                ddm_x(k,kk)=0;
                ddm_x(kk,k)=0;
            elseif abs(ddm_x(k,kk))<1.0,
                ddm_x(k,kk)=0.5;            
                ddm_x(kk,k)=0.5;
            else
                ddm_x(k,kk)=1;
                ddm_x(kk,k)=1;
            end;
        end;
    end;
else
    ddm_x=ddm_x/max(max(abs(ddm_x)));
    ddm_x=(ddm_x+1)/2;
end;

function entropy=matrix_entropy(ddm)

% Normalization

ddm=abs(ddm)/sum(sum(abs(ddm)));
lddm=ddm;
lddm(lddm>0)=log(lddm(lddm>0));
ddm=lddm.*ddm;
entropy=sum(sum(ddm));

function dissimilarity=matrix_dissimilarity(ddm1,ddm2)

[m1,~]=size(ddm1);
[m2,~]=size(ddm2);

if m1~=m2,
    similarity=[];
else
    ddm1=ddm1/max(max(abs(ddm1)));
    ddm2=ddm2/max(max(abs(ddm2)));
    dissimilarity=sum(sum(abs(ddm1-ddm2)))/(sum(sum(abs(ddm1)))+sum(sum(abs(ddm2))));
end;

