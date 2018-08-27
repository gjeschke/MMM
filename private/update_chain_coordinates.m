function update_chain_coordinates(cartfile)
% function update_chain_coordinates(cartfile)
%
% updates the coordinates of the current chain (use hierachy window to
% select it) with coordinates from a Cartesian file cartfile, which can be
% in .xyz or in Chem3D's .cc1 format
% the number and sequence of atoms in the Cartesian coordinate file must
% exactly match the one in the PDB structure from which the chain was read
% there must be exactly on chain model
% if the number does not match or if there is more than one chain model, 
% the function aborts without changing coordinates
%
% useful if some optimization program has lost biomolecule information
% (residues or atom types)
%
% if the current chain is displayed, graphics will no longer be consistent
% with coordinates after this call, the caller is responsible for redrawing
%
% use with Save selection as PDB... in the File menu to generate a PDB file
% with updated coordinates (after selecting this chain and possibly others)
%
% cartfile  name of the Cartesian coordinate file including extension
%
% G. Jeschke, 25.9.2017

global model

% determine indices
sind = model.current_structure;
cind = tag2id(model.current_chain,model.chain_tags{sind});

% check old coordinate array size
csets = length(model.structures{sind}(cind));
if csets ~=1
    add_msg_board('ERROR: The current chain must have exactly one model to replace coordinates. Aborting.');
    return
end

old_coor = model.structures{sind}(cind).xyz{1};
[ao,~] = size(old_coor);
add_msg_board(sprintf('%i atom coordinates are to be replaced.',ao));

[~,~,ext] = fileparts(cartfile);

if isempty(ext)
    add_msg_board('ERROR: The Cartesian coordinate file must be supplied with extension to determine its type. Aborting.');
    return
end

switch ext
    case '.xyz'
        ecoor = read_xyz(cartfile,false);
    case '.cc1'
        ecoor = read_xyz(cartfile,true);
    otherwise
        add_msg_board(sprintf('ERROR: Unknown extension %s of the Cartesian coordinate file. Aborting.',ext));
        return
end

[an,~] = size(ecoor);

if an ~= ao
    add_msg_board(sprintf('ERROR: Number %i of atoms in Cartesian coordinate file does not match number %i in the structure. Aborting.',an,ao));
    return
end

model.structures{sind}(cind).xyz{1} = ecoor(:,2:4); 