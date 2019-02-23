function stitch_chain(indices,nindices,chain_id,init)
% stitch_chain(indices,newindices,chain_id,init)
%
% Makes a new chain from segments, caller is responsible for the segments
% making up a continuous chain with continuous residue numbers
%
% indices   array (m,3) of chain model indices for the segments
% nindices  vector (1,3) of indices for the new chain model, the structure
%           must already exist
% chain_id  chain tag, needed only if init = true
% init      initialization flag, set true for first chain model, optional,
%           defaults to true
%
% G. Jeschke, 22.05.2018

global model

maxlen = 5000; % maximum chain length

if ~exist('init','var') || isempty(init)
    init = true;
end

snum = nindices(1);
cnum = nindices(2);
mnum = nindices(3);

if init
    model.chain_tags{snum} = strcat(model.chain_tags{snum},sprintf('%s:',chain_id));
    model.chain_ids{snum} = [model.chain_ids{snum} model.chain_ids{snum}(end)+1];
    model.structures{snum}(cnum) = model.structures{indices(1,1)}(indices(1,2));
    model.structures{snum}(cnum).name = chain_id;
    model.structures{snum}(cnum).modified = 0;
    model.structures{snum}(cnum).nonstandard = 0;
    model.structures{snum}(cnum).helices = 0;
    model.structures{snum}(cnum).strands = 0;
    model.structures{snum}(cnum).maxconn = 0;
    maxconn = 0;
    model.structures{snum}(cnum).resnum = 0;
    model.structures{snum}(cnum).sequence = char(double('?')*ones(1,maxlen));
    model.structures{snum}(cnum).restags = ':';
%     model.structures{snum}(cnum).helix = struct([]);
%     model.structures{snum}(cnum).strand = struct([]);
    model.structures{snum}(cnum).seqres = [];
    model.structures{snum}(cnum).header = '';
    model.structures{snum}(cnum).loop_defs = cell(1,0);
    model.structures{snum}(cnum).helix_defs = cell(1,0);
    model.structures{snum}(cnum).sheet_defs = cell(1,0);
end

[segments,~] = size(indices);
% determine total number of atoms and residues
atnum = 0;
resnum = 0;
nconn = 0;
hpoi = 0;
spoi = 0;
lpoi = 0;
hpoi2 = 0;
spoi2 = 0;
maxresnum = 0;
for k = 1:segments
    ssnum = indices(k,1);
    scnum = indices(k,2);
    smnum = indices(k,3);
    atnum = atnum + model.structures{ssnum}(scnum).atoms{smnum};
    resnum = resnum + length(model.structures{ssnum}(scnum).residues{smnum}.info);
    if k == 1
        residues = model.structures{ssnum}(scnum).residues{smnum};
        if isfield(residues,'secondary_graphics')
            residues = rmfield(residues,'secondary_graphics');
        end
        for ki = 1:length(residues.info)
            if ~isfield(residues.info(ki),'dssp')
                residues.info(ki).dssp = ' ';
            end
        end
    end
    if init
        model.structures{snum}(cnum).modified = model.structures{snum}(cnum).modified + model.structures{ssnum}(scnum).modified;
        model.structures{snum}(cnum).nonstandard = model.structures{snum}(cnum).nonstandard + model.structures{ssnum}(scnum).nonstandard;
        model.structures{snum}(cnum).helices = model.structures{snum}(cnum).helices + model.structures{ssnum}(scnum).helices;
        model.structures{snum}(cnum).strands = model.structures{snum}(cnum).strands + model.structures{ssnum}(scnum).strands;
        [mc,nc] = size(model.structures{ssnum}(scnum).conn);
        nconn = nconn + mc;
        if nc > maxconn 
            maxconn = nc;
        end
        model.structures{snum}(cnum).sequence = strcat(model.structures{snum}(cnum).sequence,model.structures{ssnum}(scnum).sequence);
        model.structures{snum}(cnum).restags = strcat(model.structures{snum}(cnum).restags(1:end-1),model.structures{ssnum}(scnum).restags);
        model.structures{snum}(cnum).header = strcat(model.structures{snum}(cnum).header,model.structures{ssnum}(scnum).header);
        if isfield(model.structures{ssnum}(scnum),'helix')
            nh = length(model.structures{ssnum}(scnum).helix);
            if nh > 0
                if hpoi > 0
                    model.structures{snum}(cnum).helix(hpoi+1:hpoi+nh) = model.structures{ssnum}(scnum).helix;
                else
                    model.structures{snum}(cnum).helix = model.structures{ssnum}(scnum).helix;
                end
                hpoi = hpoi + nh;
            end
        end
        if isfield(model.structures{ssnum}(scnum),'strand')
            ns = length(model.structures{ssnum}(scnum).strand);
            if ns > 0
                if spoi > 0
                    model.structures{snum}(cnum).strand(spoi+1:spoi+ns) = model.structures{ssnum}(scnum).strand;
                else
                    model.structures{snum}(cnum).strand = model.structures{ssnum}(scnum).strand;
                end
                spoi = spoi + ns;
            end
        end
        if isfield(model.structures{ssnum}(scnum),'loop_defs')
            nl = length(model.structures{ssnum}(scnum).loop_defs);
            if nl > 0
                model.structures{snum}(cnum).loop_defs(lpoi+1:lpoi+nl) = model.structures{ssnum}(scnum).loop_defs;
                lpoi = lpoi + nl;
            end
        end
        if isfield(model.structures{ssnum}(scnum),'helix_defs')
            nh = length(model.structures{ssnum}(scnum).helix_defs);
            if nh > 0
                model.structures{snum}(cnum).helix_defs(hpoi2+1:hpoi2+nh) = model.structures{ssnum}(scnum).helix_defs;
                hpoi2 = hpoi2 + nh;
            end
        end
        if isfield(model.structures{ssnum}(scnum),'sheet_defs')
            ns = length(model.structures{ssnum}(scnum).sheet_defs);
            if ns > 0
                model.structures{ssnum}(scnum).helix_defs(spoi2+1:spoi2+ns) = model.structures{ssnum}(scnum).sheet_defs;
                spoi2 = spoi2 + ns;
            end
        end
    end
end

if init
    model.structures{snum}(cnum).maxconn = maxconn;
    model.structures{snum}(cnum).conn = zeros(nconn,maxconn,'int32');
    connpoi = 0;
    model.structures{snum}(cnum).isotopes = zeros(atnum,2,'single');
    model.structures{snum}(cnum).resnum = resnum;
end
model.structures{snum}(cnum).xyz{mnum} = zeros(atnum,3);
atpoi = 0;
model.structures{snum}(cnum).Bfactor{mnum} = zeros(1,atnum);
model.structures{snum}(cnum).Btensor{mnum} = zeros(atnum,6,'int32');
model.structures{snum}(cnum).atoms{mnum} = atnum;

% second pass for filling the large arrays
for k = 1:segments
    ssnum = indices(k,1);
    scnum = indices(k,2);
    smnum = indices(k,3);
    xyz = model.structures{ssnum}(scnum).xyz{smnum};
    [nat,~] = size(xyz); 
    model.structures{snum}(cnum).xyz{mnum}(atpoi+1:atpoi+nat,:) = xyz;
    model.structures{snum}(cnum).Bfactor{mnum}(atpoi+1:atpoi+nat) = model.structures{ssnum}(scnum).Bfactor{smnum};
    model.structures{snum}(cnum).Btensor{mnum}(atpoi+1:atpoi+nat,:) = model.structures{ssnum}(scnum).Btensor{smnum};
    model.structures{snum}(cnum).atoms{mnum} = nat;
    if init
        [mc,nc] = size(model.structures{ssnum}(scnum).conn);
        model.structures{snum}(cnum).conn(connpoi+1:connpoi+mc,1:nc) = model.structures{ssnum}(scnum).conn;
        connpoi = connpoi + mc;
        model.structures{snum}(cnum).isotopes(atpoi+1:atpoi+nat,:) = model.structures{ssnum}(scnum).isotopes;
    end
    newinfo = model.structures{ssnum}(scnum).residues{smnum}.info;
    lpoi = 0;
    for kr = 1:length(newinfo)
        model.structures{snum}(cnum).sequence(newinfo(kr).number) = model.structures{ssnum}(scnum).sequence(newinfo(kr).number);
        if newinfo(kr).number > maxresnum
            maxresnum = newinfo(kr).number;
        end
        if ~isfield(newinfo(kr),'dssp')
            newinfo(kr).dssp = ' ';
        end
        atom_indices = newinfo(kr).atom_numbers;
        for ka = 1:length(atom_indices)
            cai = atom_indices{ka};
            [mat,~] = size(cai);
            cai(:,1) = zeros(size(cai(:,1))) + lpoi + atpoi;
            atom_indices{ka} = cai;
            lpoi = lpoi + mat;
        end
        newinfo(kr).atom_numbers = atom_indices;
        if k ~= segments % remove signalling of terminal nt/residue, except for the last segment
            newinfo(kr).terminal = 0;
        end
    end
    if k > 1
        residues.residue_tags = strcat(residues.residue_tags(1:end-1),model.structures{ssnum}(scnum).residues{smnum}.residue_tags);
        residues.info = [residues.info,newinfo];
    else
        residues.info = newinfo;
    end
    atpoi = atpoi + nat;
end
model.structures{snum}(cnum).residues{mnum} = residues;
model.structures{snum}(cnum).sequence = model.structures{snum}(cnum).sequence(1:maxresnum);

