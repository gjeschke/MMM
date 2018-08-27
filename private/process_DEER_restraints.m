function [DEER,cancelled]=process_DEER_restraints(restraints)
% [DEER,cancelled]=process_DEER_restraints(restraints)
%
% Processes the DEER constraints in a list obtained with rd_restraints.m
% checks, whether rotamers for the requested label type and temperature
% were already computed, computes them, if not
% rotamers are not actually attached
% residue indices into the template structure, template label coordinates
% (mean N-O midpoint coordinates), and label position r.m.s.d. are computed
% and stored
%
% restraints    constraint structure (obtained with rd_restraints.m)
%
% DEER          array of result structure, length is number of constraints,
%               fields are:
%               .adr1   address of the first site
%               .ind1   residue indices of the first site
%               .adr2   address of the second site
%               .ind2   residue indices of the second site
%               .r      constraint on mean distance
%               .sigr   uncertainty of constraint
%               .xyz1   Cartesian coordinates of mean N-O midpoint of first
%                       site
%               .xyz2   Cartesian coordinates of mean N-O midpoint of
%                       second site
%               .rmsd1  r.m.s.d of mean N-O midpoint coordinate of first
%                       site
%               .rmsd2  r.m.s.d. of mean N-O midpoint coordinate of
%                       second site
% cancelled     flag that reporst on interactive cancellation of
%               processing, true cancelled by user, false not cancelled
%
% G. Jeschke 2010-2012

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

for k=1:length(restraints.DEER),
    adr1=restraints.DEER(k).adr1;
    ind1=resolve_address(adr1);
    adr2=restraints.DEER(k).adr2;
    ind2=resolve_address(adr2);
    DEER(k).r=restraints.DEER(k).r;
    DEER(k).sigr=restraints.DEER(k).sigr;
    DEER(k).ind1=ind1;
    DEER(k).ind2=ind2;
    DEER(k).adr1=adr1;
    DEER(k).adr2=adr2;
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
