function fname=prepare_Modeller_fit(rfile,template,target)
% function fname=prepare_Modeller_fit(rfile)
%
% prepares a PDB template with spin labels, an alignment file, and a Python
% script for Modeller for restraint-driven homology modelling, restraints
% are provided in MMM restraint file format
%
% rfile     restraint file name (MMM format)
% template  template identifier, used for template PDB filename
% target    target identifier, used for output file name
%
% fname     name of the Modeller Python script file
%
% G. Jeschke, 2011

global model
global hMain

fname=['model_' template '_to_' target '.py'];

restraints=rd_restraints(rfile);

if ~isfield(restraints,'DEER'),
    DEER=[];
    return;
end;

snum=model.current_structure;

if isfield(restraints,'PDB'),
    if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
        button = questdlg(sprintf('Restraint file specifies template %s, while current template is %s. Do you want to continue?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
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
    DEER(k).indices=[ind1;ind2];
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

