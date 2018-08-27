function [ merit, message ] = analyze_model_fom( initial_num, ENM_num, hom_num )
%ANALYZE_MODEL_FOM Determines figure of merit for a structural transition
%                  fit, if the initial structure and two models are given
%   initial_num (optional) structure number for initial state
%   ENM_num     (optional) structure number for first model
%   hom_num     (optional) structure number for second model
%               network model
%   merit       figure of merit for model agreement
%   message     error message with fields
%               .error: 0 no error, otherwise error code
%               .text: 'Initial structure missing.', error code 1
%                      'First model missing.', error code 3
%                      'Second model missing.', error code 4
%                      'Cancelled.', error code 5
%                      'Residue number mismatch.', error code 6
%
% G. Jeschke, 14.1.2013

global general

merit = [];

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
    [message,initial_num]=add_pdb(initial);
    if message.error,
        message.error=1;
        message.text='Initial structure missing.';
        return
    end;
end;

if nargin<2,
    [FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file for first model');
    if isequal(FileName,0) || isequal(PathName,0),
        add_msg_board('Loading of PDB file canceled by user.');
        message.error=5;
        message.text='Cancelled.';
        return;
    end;
    reset_user_paths(PathName);
    general.pdb_files=PathName;
    final=fullfile(PathName,FileName);
    [message,ENM_num]=add_pdb(final);
    if message.error,
        message.error=3;
        message.text='First model missing.';
        return
    end;
end;

if nargin<3,
    [FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file for second model');
    if isequal(FileName,0) || isequal(PathName,0),
        add_msg_board('Loading of PDB file canceled by user.');
        message.error=5;
        message.text='Cancelled.';
        return;
    end;
    reset_user_paths(PathName);
    general.pdb_files=PathName;
    final=fullfile(PathName,FileName);
    [message,hom_num]=add_pdb(final);
    if message.error,
        message.error=3;
        message.text='ENM fit missing.';
        return
    end;
end;

initial_coor = coarse_residues_and_store(initial_num);
if isempty(initial_coor),
    add_msg_board('Initial structure could not be coarse grained');
    message.error=1;
    message.text='Initial structure missing.';
    return
end

ENM_coor = coarse_residues_and_store(ENM_num);
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
esize=size(ENM_coor);
hsize=size(hom_coor);

if isize(2)~=3,
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Calpha coordinates of initial structure are improper.');
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

correspondence_E=align_template_target(initial_num,ENM_num);
initial_coor_E=initial_coor(correspondence_E(1,:),:);
ENM_coor_i=ENM_coor(correspondence_E(2,:),:);
correspondence_h=align_template_target(initial_num,hom_num);
initial_coor_h=initial_coor(correspondence_h(1,:),:);
hom_coor_i=hom_coor(correspondence_h(2,:),:);
correspondence_M=align_template_target(ENM_num,hom_num);
ENM_coor_M=ENM_coor(correspondence_M(1,:),:);
hom_coor_M=hom_coor(correspondence_M(2,:),:);
isize_E=size(initial_coor_E);
isize_h=size(initial_coor_h);
esize_i=size(ENM_coor_i);
hsize_i=size(hom_coor_i);
esize_M=size(ENM_coor_M);
hsize_M=size(hom_coor_M);
if isize_E(1)~=esize_i(1),
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Residue number of first model does not match residue number of initial structure.');
    return
end;
if isize_h(1)~=hsize_i(1),
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Residue number of second model does not match residue number of initial structure.');
    return
end;
if esize_M(1)~=hsize_M(1),
    message.error = 6;
    message.text = 'Residue number mismatch.';
    add_msg_board('Residue number of first model does not match residue number of second model.');
    return
end;

% at this point we have three pairs of two coordinate arrays of the same size with
% aligned residues

rmsd_iE=rmsd_superimpose(initial_coor_E,ENM_coor_i);

add_msg_board(sprintf('Calpha r.m.s.d. between initial structure and first model is %6.2f Å',rmsd_iE));

rmsd_ih=rmsd_superimpose(initial_coor_h,hom_coor_i);

add_msg_board(sprintf('Calpha r.m.s.d. between initial structure and second model is %6.2f Å',rmsd_ih));

rmsd_Eh=rmsd_superimpose(ENM_coor_M,hom_coor_M);

add_msg_board(sprintf('Calpha r.m.s.d. between first and second model is %6.2f Å',rmsd_Eh));

s=(rmsd_iE+rmsd_ih+rmsd_Eh)/2;
th=2*asin(sqrt((s-rmsd_iE)*(s-rmsd_ih)/(rmsd_iE*rmsd_ih)));
arg=th^-10;
merit=3.3103*arg-3.1506;

add_msg_board(sprintf('Figure of merit for the model pair is %6.3f',merit));


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

