function [PDBind,Tink2PDB,snum,msg] = prepare_tinker_input(basname,molecule,selected,options)
% [PDBind,Tink2PDB,msg] = prepare_tinker_input(basname,molecule,selected,options)
%
% Writes Tinker input .xyz and .key files
%
% basname   name of the input molecular system file, if empty, a default
%           name tinker_mol is used
% molecule  cell array of index arrays of objects that make up the whole 
%           input molecule, defaults to current structure
% selected  cell array of index arrays of objects that are active or 
%           inactive during the Tinker computation, defaults to the current 
%           selection
% options   Tinker computation options, contain also restraints, see
%           wr_tinker_key source code for possibilities and defaults
%
% selection list of the selected atoms (indices) in the Tinker .xyz file
% PDBind    index array of location indices corresponding to the lines of
%           the PDB input file
% Tink2PDB  vector pointers into PDBind corresponding to the lines in the 
%           Tinker xyz file, this is required as Tinker pdbxyz adds and 
%           removes atoms
% snum      number of the structure corresponding to the input PDB file
% msg       message on success, msg.error is error code (0 for success) and
%           msg.text the error message
%
% G. Jeschke, 12.7.2017

global model

if isfield(model,'selected')
    old_sel = model.selected;
else
    old_sel = [];
end

msg.error = 0;
msg.text = 'Tinker input files written';

if ~exist('basname','var') || isempty(basname)
    basname = 'tinker_mol';
end


if ~exist('molecule','var') || isempty(molecule)
    molecule{1} = model.current_structure;
end

if ~exist('selected','var') || isempty(selected)
    selected = model.selected;
end

if isfield(model,'selected')
    curr_sel = model.selected;
else
    curr_sel = [];
end

if ~exist('options','var') || isempty(options)
    options.active = true;
end

if ~isfield(options,'forcefield')
    options.forcefield = 'amber99';
end

tinker_path = get_tinker_path;

if isempty(tinker_path)
    msg.error = 1;
    msg.text = 'Tinker distribution (file amber99.prm) not found on Matlab path.';
    return
end

delete(fullfile(tinker_path,strcat(basname,'*.*')));

pdbfile = fullfile(tinker_path,strcat(basname,'.pdb'));

model.selected = selected;

sel_atoms = sorted_selection;

if ~isempty(sel_atoms)
    % check whether there are multiple locations of an atom, which would poison
    % Tinker
    multloc = sum(abs(sel_atoms(:,6)-1));
    if multloc
        msg.error = 2;
        msg.text = 'Selection contains multiple locations of an atom.';
        return
    end
end

model.selected = molecule;

[message,PDBind0,PDBind] = wr_pdb_selected(pdbfile,'TINK');

[msg,snum] = add_pdb(pdbfile);

multichain = false;
indi = PDBind0(1,1:2);
[nat,~] = size(PDBind0);
for kk = 2:nat
    if sum(abs(indi-PDBind0(kk,1:2)))
        multichain = true;
    end
end

if message.error
    add_msg_board(strcat('ERROR (wr_pdb_selected): ',message.text));
    return
end

% check whether there are multiple locations of an atom, which would poison
% Tinker
multloc = sum(abs(PDBind0(:,6)-1));
if multloc
    msg.error = 3;
    msg.text = 'Input molecule contains multiple locations of an atom.';
    return
end

if ~isempty(curr_sel)
    model.selected = curr_sel;
else
    model = rmfield(model,'selected');
end

dospath=which('pdbxyz.exe');
if isempty(dospath)
    add_msg_board('Tinker module pdbxyz not found on Matlab path.');
    return
end
[modpath, modcmd] = fileparts(dospath);
add_msg_board('Now calling Tinker pdbxyz');

cmd = sprintf('pdbxyz ..%s%s ..%sparams%s%s.prm',filesep,basname,filesep,filesep,options.forcefield);

if multichain
    cmd = sprintf('pdbxyz ..%s%s ALL ..%sparams%s%s.prm',filesep,basname,filesep,filesep,options.forcefield);
end

my_dir = pwd;
cd(modpath);

tic,
[s, w] = dos(cmd);
runtime=toc;
add_msg_board(sprintf('Tinker pdbxyz was running %5.1f s\n',runtime));
if s~=0
    rem=w;
    while ~isempty(rem)
        [token,rem]=strtok(rem,newline);
        if ~isempty(token)
            add_msg_board(token);
        end
    end
    add_msg_board('ERROR: Tinker pdbxyz did not run successfully.');
    cd(my_dir);
    return
end

ecoor = read_tinker_xyz(basname);
[nt,~] = size(ecoor);
tcoor = ecoor(:,2:4);

[ns,~] = size(sel_atoms);
selection = zeros(1,ns);
for k = 1:ns
    [msg,ccoor] = get_location(sel_atoms(k,:),'xyz');
    if msg.error
        add_msg_board('ERROR: Coordinates of selected atom could not be retrieved.');
        cd(my_dir);
        return
    end
    comp = repmat(ccoor,nt,1);
    match = sum(abs(tcoor-comp),2);
    poi = find(match < 1e-3);
    if length(poi) > 1 % shouldn't happen, but better safe than sorry
        msg.error = 4;
        msg.text = 'Several lines in Tinker xyz file match the same selected atom.';
        cd(my_dir);
        return
    end
    if isempty(poi)
        msg.error = -1;
        msg.text = 'A selected atom is missing in the molecule';
    else
        selection(k) = poi;
    end
end

Tink2PDB = zeros(1,nt);
[np,~] = size(PDBind0);
for k = 1:np
    [msg,ccoor] = get_location(PDBind0(k,:),'xyz');
    if msg.error % should not happen
        add_msg_board('ERROR: Coordinates of an atom could not be retrieved.');
        cd(my_dir);
        return
    end
    comp = repmat(ccoor,nt,1);
    match = sum(abs(tcoor-comp),2);
    poi = find(match < 1e-3);
    if length(poi) > 1 % shouldn't happen, but better safe than sorry
        msg.error = 5;
        msg.text = 'Several lines in Tinker xyz file match the same atom of the PDB file.';
        cd(my_dir);
        return
    end
    if ~isempty(poi)
        Tink2PDB(poi) = k;
    end
end

if isfield(options,'restraints')
    options.restrain_d = zeros(length(options.restraints),5);
    for kr = 1:length(options.restraints)
        ccoor = options.restraints(kr).CA1;
        comp = repmat(ccoor,nt,1);
        match = sum(abs(tcoor-comp),2);
        poi1 = find(match < 1e-3);
        ccoor = options.restraints(kr).CA2;
        comp = repmat(ccoor,nt,1);
        match = sum(abs(tcoor-comp),2);
        poi2 = find(match < 1e-3);
        if length(poi1) > 1 % shouldn't happen, but better safe than sorry
            msg.error = 5;
            msg.text = 'Several lines in Tinker xyz file match the same CA1 atom of a restraint.';
            cd(my_dir);
            return
        end
        if length(poi2) > 1 % shouldn't happen, but better safe than sorry
            msg.error = 5;
            msg.text = 'Several lines in Tinker xyz file match the same CA2 atom of a restraint.';
            cd(my_dir);
            return
        end
        if isempty(poi1) % shouldn't happen, but better safe than sorry
            msg.error = 6;
            msg.text = 'CA1 atom of a restraint not found.';
            cd(my_dir);
            return
        end
        if isempty(poi2) % shouldn't happen, but better safe than sorry
            msg.error = 6;
            msg.text = 'CA2 atom of a restraint not found.';
            cd(my_dir);
            return
        end
        options.restrain_d(kr,1) = poi1;
        options.restrain_d(kr,2) = poi2;
        options.restrain_d(kr,3) = 25/options.restraints(kr).sigr;
        options.restrain_d(kr,4) = options.restraints(kr).rCA - options.restraints(kr).sigrCA;
        options.restrain_d(kr,5) = options.restraints(kr).rCA + options.restraints(kr).sigrCA;
    end
end

if options.active
    inserted = find(Tink2PDB==0);
    newselection = [selection inserted];
    selection = sort(newselection);
end


msg = wr_tinker_key(basname,options,selection);

if msg.error
    add_msg_board(sprintf('ERROR (wr_tinker_key): %s',msg.text));
    cd(my_dir);
    return
end


cd(my_dir);

if ~isempty(old_sel)
    model.selected = old_sel;
elseif isfield(model,'selected')
    model = rmfield(model,'selected');
end


