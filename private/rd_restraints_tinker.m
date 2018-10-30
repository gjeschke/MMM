function [restraints,failed] = rd_restraints_tinker(fname)
% function [restraints,failed] = rd_restraints_tinker(fname)
%
% reads TINKER job definition and restraints from an ASCII file
% restraints are returned in a structure whose fields correspond to
% different restraint types
%
% the restraints are pre-processed with respect to the current structure in
% MMM
%
% currently implemented types are:
%
% TINKER    tinker command and parameter, current options:
%           minimize thr   % thr is the rms gradient convergence threshold
%           optimize thr   % thr is the rms gradient convergence threshold
% INACTIVE  block key, residue ranges that are inactive during minimization
%           useful for defining rigid bodies
%           each subsequent line defines a continous segment of inactive
%           residues
% RANDOMIZE defines the restraint range for an individual conformation as a
%           fraction of the experimental restraint width for the ensemble,
%           argument N is the number of intervals within +/- the standard
%           deviation, mean values conform to a Gaussian distribution
% DEER      distance distribution restraints between spin labels (block
%           key), restraints are translated to Calpha-Calpha restraints
% UNITS         can be NM (nanometers, default) or A (Angstroem)
%
% G. Jeschke, 24.10.2018

global model

failed = true;

fid=fopen(fname);
if fid==-1
    add_msg_board('ERROR: Restraint file does not exist');
    return;
end

clear restraints
restraints.DEER(1).r = [];
restraints.DEER(1).sigr = [];
restraints.DEER(1).label = '';
restraints.DEER(1).CA1 = [];
restraints.DEER(1).CA2 = [];
restraints.DEER(1).r0 = [];
restraints.DEER(1).sigr0 = [];
restraints.DEER(1).rCA0 = [];
restraints.DEER(1).rCA = [];
restraints.DEER(1).sigrCA = [];

restraints.randomize = 0;

restraints.PDB = '';

restraints.inactive(1).initial = [];
restraints.inactive(1).final = [];

DEER_poi=0;
inactive_poi=0;

label_adr = ':';
label_types = ':';
label_coor = zeros(1000,4);
label_indices = zeros(1000,4);
mode=0;
scale_units=1;

nl = 0;
while 1
    tline = fgetl(fid);
    nl = nl + 1;
    if ~ischar(tline) || mode<0, break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            if k(1)>1
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end
        end
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#')
            switch upper(char(args(2)))
                case 'PDB'
                    mode=0;
                    restraints.current = char(args(3));
                    restraints.PDB = char(args(4));
                case 'RANDOMIZE'
                    mode=0;
                    restraints.randomize = str2double(char(args(3)));
                case 'TINKER'
                    mode=0;
                    restraints.tinker.cmd = char(args(3));
                    if length(args) > 3
                        restraints.tinker.thr = str2double(char(args(4)));
                    else
                        restraints.tinker.thr = 0.1;
                    end
                case 'DEER'
                    mode=1;
                    label = char(args(3));
                    label = check_label(label);
                    if isempty(label)
                        add_msg_board(sprintf('Warning: Label %s is unknown in line %i. Reverting to MTSL.',char(args(3))),nl);
                        label = 'R1A';
                    end
               case 'INACTIVE'
                    mode=2;
               case 'END'
                    mode=-1;
                otherwise
                    mode=0;
                    add_msg_board(sprintf('Warning: Unknown restraint mode %s',upper(char(args(2)))));
            end
        elseif mode>0 && ~strncmp(strtrim(char(args(1))),'%',1)
% %             keyboard
            switch mode
                case 1 % DEER
                    DEER_poi = DEER_poi + 1;
                    restraints.DEER(DEER_poi).r = 10*scale_units*str2double(char(args(3)));
                    restraints.DEER(DEER_poi).sigr = 10*scale_units*str2double(char(args(4)));
                    restraints.DEER(DEER_poi).label = label;
                    restraints.DEER(DEER_poi).adr1 = char(args(1));
                    restraints.DEER(DEER_poi).adr2 = char(args(2));
                    [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(char(args(1)),label,label_adr,label_types,label_coor,label_indices);
                    if isempty(ecoor)
                        add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(1))));
                        fclose(fid);
                        return
                    end
                    restraints.DEER(DEER_poi).ecoor1 = ecoor;
                    restraints.DEER(DEER_poi).indices1 = indices;
                    [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(char(args(2)),label,label_adr,label_types,label_coor,label_indices);
                    if isempty(ecoor)
                        add_msg_board(sprintf('ERROR: Residue %s could not be labeled. Aborting.',char(args(2))));
                        fclose(fid);
                        return
                    end
                    restraints.DEER(DEER_poi).ecoor2 = ecoor;
                    restraints.DEER(DEER_poi).indices2 = indices;
                    restraints.DEER(DEER_poi).r0 = norm(restraints.DEER(DEER_poi).ecoor1(1:3)-restraints.DEER(DEER_poi).ecoor2(1:3));
                    restraints.DEER(DEER_poi).sigr0 = sqrt(restraints.DEER(DEER_poi).ecoor1(4)^2+restraints.DEER(DEER_poi).ecoor2(4)^2);
                    adr1 = sprintf('%s.CA',char(args(1)));
                    [indices,message]=resolve_address(adr1);
                    if message.error
                        add_msg_board(sprintf('ERROR: Calpha atom of residue %s does not exist. Aborting.',adr1));
                        fclose(fid);
                        return
                    end
                    [message,coor1]=get_atom(indices,'coor');
                    if message.error || isempty(coor1)
                        add_msg_board(sprintf('ERROR: Coordinates of Calpha atom of residue %s could not be retrieved. Aborting.',adr1));
                        fclose(fid);
                        return
                    end
                    restraints.DEER(DEER_poi).CA1 = coor1;
                    adr2 = sprintf('%s.CA',char(args(2)));
                    [indices,message]=resolve_address(adr2);
                    if message.error
                        add_msg_board(sprintf('ERROR: Calpha atom of residue %s does not exist. Aborting.',adr1));
                        fclose(fid);
                        return
                    end
                    [message,coor2]=get_atom(indices,'coor');
                    if message.error || isempty(coor1)
                        add_msg_board(sprintf('ERROR: Coordinates of Calpha atom of residue %s could not be retrieved. Aborting.',adr2));
                        fclose(fid);
                        return
                    end
                    restraints.DEER(DEER_poi).CA2 = coor2;
                    restraints.DEER(DEER_poi).rCA0 = norm(coor1-coor2);
                    corr = restraints.DEER(DEER_poi).rCA0 - restraints.DEER(DEER_poi).r0;
                    restraints.DEER(DEER_poi).rCA = restraints.DEER(DEER_poi).r + corr;
                    if restraints.DEER(DEER_poi).sigr > restraints.DEER(DEER_poi).sigr0
                        restraints.DEER(DEER_poi).sigrCA = sqrt(restraints.DEER(DEER_poi).sigr^2 - restraints.DEER(DEER_poi).sigr0^2);
                    else
                        restraints.DEER(DEER_poi).sigrCA = 0.1;
                    end
                case 2 % INACTIVE
                    inactive_poi = inactive_poi + 1;
                    adr1 = char(args(1));
                    [indices,message]=resolve_address(adr1);
                    if message.error || length(indices) ~= 4
                        add_msg_board(sprintf('ERROR: Residue %s does not exist. Aborting.',adr1));
                        fclose(fid);
                        return
                    end
                    restraints.inactive(inactive_poi).initial = indices;    
                    adr2 = char(args(2));
                    [indices,message]=resolve_address(adr2);
                    if message.error || length(indices) ~= 4
                        add_msg_board(sprintf('ERROR: Residue %s does not exist. Aborting.',adr2));
                        fclose(fid);
                        return
                    end
                    restraints.inactive(inactive_poi).final = indices;   
            end
        end
    end
end

    
fclose(fid);

if restraints.randomize
    for kr = 1:length(restraints.DEER)
        rmean0 = restraints.DEER(kr).r;
        sigr0 = restraints.DEER(kr).sigr;
        rmean = rmean0 + sigr0*randn/sqrt(2);
        sigr = sigr0/restraints.randomize;
        restraints.DEER(kr).r_rand = rmean;
        restraints.DEER(kr).sigr_rand = sigr;
        corr = restraints.DEER(kr).rCA0 - restraints.DEER(kr).r0;
        restraints.DEER(kr).rCA = rmean + corr;
        if sigr > restraints.DEER(kr).sigr0
            restraints.DEER(kr).sigrCA = sqrt(sigr^2 - restraints.DEER(kr).sigr0^2);
        else
            restraints.DEER(kr).sigrCA = 0.1;
        end
    end
end

failed = false;

function label_std = check_label(label)
% checks if a rotamer library exists for a requested spin label and returns
% the MMM-internal three-letter code for this label
% an empty string is returned, if the label does not exist

global rotamer_libraries

label_std = '';

for k = 1:length(rotamer_libraries)
    if strcmpi(label,rotamer_libraries(k).label) || strcmpi(label,rotamer_libraries(k).tc)
        label_std = rotamer_libraries(k).tc;
    end
end

function [ecoor,indices,label_adr,label_types,label_coor,label_indices] = get_label_coor(adr,label,label_adr,label_types,label_coor,label_indices)
% returns chain number, label coordinate, and label position uncertainty as
% well as labeled residue indices
% (ecoor empty, if not defined) 

global hMain
global model

ml = sum(label_coor(:,1) ~= 0);

labnum = [];
labnums = tag2ids(adr,label_adr);
if ~isempty(labnum)
    for k = 1:length(labnums)
        type = id2tag(labnums(k),label_types);
        if strcmpi(label,type)
            labnum = labnums(k);
        end
    end
end

if ~isempty(labnum)
    ecoor = label_coor(labnum,:);
    indices = label_indices(labnum,:);
else
    [indices,message] = resolve_address(adr);
    if message.error % this could be a label site in a stemloop library
        return
    else
        command=sprintf('rotamers %s %s %i',adr,label,298);
        hMain.store_undo=false;
        hMain.dynamic_rotamers=false;
        cmd(hMain,command);
        labels = label_information(model.sites);
        found = false;
        ecoor = zeros(1,4);
        for k = 1:length(labels)
            if abs(sum(indices - labels(k).indices)) == 0
                found = true;
                ecoor(1:3) = labels(k).xyz;
                ecoor(4) = labels(k).rmsd;
            end
        end
        if found
            label_coor(ml+1,:) = ecoor;
            label_indices(ml+1,:) = indices;
            label_adr = [label_adr adr ':'];
            label_types = [label_types label ':'];
        else
            ecoor = [];
        end
    end
end
    

function labels = label_information(sites)

global model
global label_defs

poi=0;
for k0=1:length(sites)
    for k1=1:length(sites{k0})
        for k=1:length(sites{k0}(k1).residue)
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
        end
    end
end

function [krb,kp] = identify(indices,label,restraints)

krb = 0;
kp = 0;
for kr = 1:length(restraints.pindices)
    cindices = restraints.pindices{kr};
    if min(abs(indices(2)-restraints.rb(kr).chains)) == 0 % chain is in this rigid body
        krb = kr;
        [mc,~] = size(cindices);
        for kpc = 1:mc
            if abs(sum(cindices(kpc,:)-indices)) == 0 % same residue indices
                clabels = restraints.plabels{kr};
                if strcmpi(label,clabels{kpc}) % same spin label
                    kp = kpc;
                end
            end
        end
    end
end

function [rmsd,xyz]=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
xyz = [xmean,ymean,zmean];
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1)); 



