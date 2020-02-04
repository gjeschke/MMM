function ensemble_standard_view(inname,outname,selection)
% ensemble_standard_view(fname,selection)
%
% transforms coordinates of conformers in an ensemble to a standard frame
% where the conformers are superimposed with minimum rmsd and the Cartesian
% coordinate axes are principal axes of the tensor of inertia of the whole
% superimposed ensemble
%
% inname    file name of the original ensemble PDB file
% outname   file name of the output file
% selection (optional) selection of the objects to be superimposed in MMM
%           address format, defaults to all atoms, should not contain the
%           structure tag

%
% G. Jeschke, 16.1.2020

global model


% remove extension from input file name, if present
s = strfind(inname,'.pdb');
if ~isempty(s)
    inname = inname(1:s-1);
end

% generate output file name, if not provided
if ~exist('outname','var') || isempty(outname)
    outname = strcat(inname,'_std_view');
end
% add extension to output file name
if ~contains(outname,'.pdb')
    outname = strcat(outname,'.pdb');
end

% select all atoms, if no selection is provided
if ~exist('selection','var') || isempty(selection)
    selection = '(:){1}:.:';
end

% load the input file
[message,snum]=add_pdb(inname);
if message.error
    add_msg_board(sprintf('ERROR (ensemble_standard_view) reading input file: %s',message.text));
    return
end

% determine structure tag
stag = mk_address(snum);

% determine atom indices
sel_indices = resolve_address(strcat(stag,selection));

% number of atoms
[m,~] = size(sel_indices);
if m < 4
    add_msg_board(sprintf('ERROR (ensemble_standard_view): less than 4 (%i) atoms were selected',m));
    return
end

% number of models in the ensemble
n = length(model.structures{snum}(1).xyz);

% extract atom ccordinates of selection for all models
coor = cell(1,n);
for k = 1:n
    sel_indices(:,3) = k;
    [~,coor0] = get_object(sel_indices,'xyz');
    coor0 = cell2mat(coor0);
    coor0 = reshape(coor0,[3,m]);
    coor{k} = coor0';
end

% make the pairwise rmsd matrix and the transformation matrices
pair_rmsd = zeros(n);
trafo = cell(n);
for k1 = 1:n
    for k2 = 1:n
        [rms,~,transmat] = rmsd_superimpose(coor{k1},coor{k2});
        pair_rmsd(k1,k2) = rms;
        trafo{k1,k2} = transmat;
    end
end

% determine the central model
[~,central] = min(sum(pair_rmsd));

% superimpose the selected atoms of all models to the one of the central
% model
for k = 1:n
    for kc = 1:length(model.structures{snum})
        coor0 = model.structures{snum}(kc).xyz{k};
        coor1 = affine_trafo_coor(coor0,trafo{central,k});
        model.structures{snum}(kc).xyz{k} = coor1;
    end
end

% determine the tensor of inertia of the whole ensemble
all_coor = zeros(1000000,3);
poi = 0;
for km = 1:n
    for kc = 1:length(model.structures{snum})
        coor = model.structures{snum}(kc).xyz{k};
        [a,~] = size(coor);
        all_coor(poi+1:poi+a,:) = coor;
        poi = poi + a;
    end
end
all_coor = all_coor(1:poi,:);
inertia = inertia_tensor(all_coor);
% diagonalize the tensor of inertia and make the affine transformation
% matrix towards its eigenframe
[evec,~] = eig(inertia); 
trafo = eye(4);
trafo(1:3,1:3) = evec';

% transform the ensemble to the eigenframe of its tensor of inertia
for k = 1:n
    for kc = 1:length(model.structures{snum})
        coor0 = model.structures{snum}(kc).xyz{k};
        coor1 = affine_trafo_coor(coor0,trafo);
        model.structures{snum}(kc).xyz{k} = coor1;
    end
end

% write output PDB file
model.current_structure = snum;
idCode = model.info{snum}.idCode;
if isempty(idCode)
    idCode = 'MMM0';
end

wr_pdb(outname,idCode);


