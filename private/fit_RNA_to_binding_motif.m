function conformation = fit_RNA_to_binding_motif(decoys,motif,nts,rrmcoor,used,offset,libdir,chaintag,conformation,options)
% Fits nucleotide base heavy-atom coordinates of Rosetta RNA decys to
% corresponding coordinates in a binding motif determined by an NMR or
% CLIR/MS study
% use the Rosie FarFar server to generate decoys: 
% http://rosie.rosettacommons.org/rna_denovo
%
% decoys    path name for the directory with decoys to be tested
% motif     MMM indices of the B nucleotides in the binding motif, (B,3)
%           array
% nts       nucleotide (nt) numbers of the corresponding nucleotides in the
%           Rosetta decoys
% rrmcoor   coordinates of RRM heavy atoms for clash test
% used      nts that are actually used, these are absolute numbers,
%           starting from 1 for the first defined nt
% offset    nt number offset for storing the fitted conformations
% libdir    directory where the library should be stored
% chaintag  tag for the RNA chain
% conformation  number of already existing conformations, defaults to zero
% options   determines whether the models are optimized by Tinker,
%           provides optimization options, and can define an acceptance
%           threshold for the fit of the binding motif, if the acceptance
%           threshold is negative, the model is only compiled and written,
%           not fitted
%
% G. Jeschke, 19.03.2018-25.10.2019

acceptance = 2.0; % acceptance threshold, RMSD (Å) for overlap of bases of the binding nts

if ~exist('conformation','var') || isempty(conformation)
    conformation = 0;
end

if ~exist('options','var') || isempty(options)
    options.optimize = false;
else
    if isfield(options,'acceptance') && ~isempty(options.acceptance)
        acceptance = options.acceptance;
    end
end

[B,~] = size(motif);
B2 = length(nts);
if B ~= B2
    add_msg_board(sprintf('ERROR: Number %i of nucleotides in binding motif and of binding nucleotides in decoys (%i) do not match. Aborting.',B,B2));
    return
end



base_atoms = zeros(1,B);

set(gcf,'Pointer','watch');

mcoor = zeros(11*B,3); % allocate sufficient array for base coordinates
mpoi = 0;

% retrieve base heavy atom coordinates of the binding motif
for kb = 1:B
    indi = motif(kb,:);
    [msg, xyz] = get_object(indi,'xyz_base');
    if msg.error
        adr = mk_address(indi);
        add_msg_board(sprintf('ERROR: Base coordinates of nucleotide %s could not be retrieved: %s. Aborting.',adr,msg.text));
        return
    end
    [na,~] = size(xyz);
    base_atoms(kb) = na;
    mcoor(mpoi+1:mpoi+na,:) = xyz;
    mpoi = mpoi + na;
end
mcoor = mcoor(1:mpoi,:);

repname = fullfile(libdir,'report.txt');
[report,errmsg] = fopen(repname,'at');
libdir = fullfile(pwd,libdir);
mydir = pwd;
cd(decoys);

file_list = get_file_list('diverse_SLs.dat');
if ~isempty(file_list)
    for kf = 1:length(file_list)
        RNA_list(kf).name = file_list{kf};
    end
else
    RNA_list = dir('*.pdb');
end

ndecoys = length(RNA_list);

tic,
for kd = 1:ndecoys
    add_msg_board(sprintf('Processing decoy %s.',RNA_list(kd).name));
    [msg,snum] = add_pdb(RNA_list(kd).name);
    if msg.error
        add_msg_board(sprintf('ERROR: PDB file %s could not be read (%s). Aborting.',RNA_list(kd).name,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    % write copy without fitting, if requested
    if acceptance < 0
        write_unfitted_rna(snum,nts,used,offset,libdir,chaintag,conformation,options);
        conformation = conformation + 1;
        continue
    end
    stag = mk_address_parts(snum);
    dpoi = 0;
    dcoor = zeros(11*B,3); % allocate sufficient array for base coordinates
    C1p = zeros(B,3);
    for kb = 1:B
        adr = sprintf('[%s]%i',stag,nts(kb));
        [msg, xyz] = get_object(adr,'xyz_base');
        if msg.error
            add_msg_board(sprintf('ERROR: Base coordinates of nucleotide %s could not be retrieved: %s. Aborting.',adr,msg.text));
            return
        end
        [na,~] = size(xyz);
        dcoor(dpoi+1:dpoi+na,:) = xyz;
        dpoi = dpoi + na;
        adr = sprintf('%s.C1''',adr);
        [msg, coor] = get_object(adr,'coor');
        if msg.error
            add_msg_board(sprintf('ERROR: Base coordinates of nucleotide %s could not be retrieved: %s. Aborting.',adr,msg.text));
            return
        end
        C1p(kb,:) = coor;
    end
    dcoor = dcoor(1:dpoi,:);
    if dpoi ~= mpoi % should never happen
        add_msg_board(sprintf('ERROR": Numbers of base heavy atoms in binding motif (%i) and decoy (%i) do not match. Aborting.',mpoi,dpoi));
        fprintf(1,'Why does this happen in decoy %i?\n',kd);
        return
    end
    [~,dcoor1,transmat] = rmsd_superimpose(mcoor,dcoor);
    bas = 0;
    for kb = 1:B
        btcoor = mcoor(bas+1:bas+base_atoms(kb),:);
        bcoor = dcoor1(bas+1:bas+base_atoms(kb),:);
        nbcoor = base_rotation(btcoor,bcoor,C1p(kb,:));
        dcoor1(bas+1:bas+base_atoms(kb),:) = nbcoor;
        bas = bas + base_atoms(kb);
    end
    match = sqrt(sum(sum((mcoor-dcoor1).^2))/mpoi);
    fprintf(1,'%s/%s: Match = %5.2f Å at acceptance threshold %5.2f Å\n',libdir,RNA_list(kd).name,match,acceptance);
    if match < acceptance
        [clash,options] = write_fitted_rna(snum,dcoor,dcoor1,transmat,rrmcoor,nts,used,offset,libdir,chaintag,conformation,options);
        if ~clash && ~isnan(options.energy)
            conformation = conformation + 1;
            fprintf(report,'%i: Match: %5.2f Å derived from %s',conformation,match,fullfile(decoys,RNA_list(kd).name));
            if options.optimize
                fprintf(report,'; energy %12.1f kcal/mol',options.energy);
            end
            fprintf(report,'\n');
        end
    end
end
set(gcf,'Pointer','arrow');

fclose(report);

cd(mydir);


function ncoor = base_rotation(tcoor,coor,C1p)

options.MaxIter = 100;
options.TolX = 1e-3;
fi = fminbnd(@(fi) rmsd_base_rot(fi,tcoor,coor,C1p),0,2*pi-eps,options);
ncoor = base_rot(fi,coor,C1p);

function rmsd = rmsd_base_rot(fi,tcoor,coor,C1p)

[m,~] = size(tcoor);
ncoor = base_rot(fi,coor,C1p);
rmsd = sqrt(sum(sum((tcoor-ncoor).^2))/m);

function ncoor = base_rot(fi,coor,C1p)

[m,~] = size(coor);
rotax = coor(1,:) - C1p;
transmat=affine('rotn',fi,rotax);
ncoor = [coor ones(m,1)];
ncoor = ncoor*transmat';
ncoor = ncoor(:,1:3);

function write_fitted_rna(snum,dcoor,dcoor1,transmat,rrmcoor,nts,used,offset,libdir,chaintag,conformation,fit_options)

global model

fit_options.energy = 0;

[stag,ctag] = mk_address_parts([snum 1]);
adr0 = sprintf('[%s](%s)',stag,ctag);
nti = zeros(1,length(nts));
for k = 1:length(nts)
    adr = sprintf('%s{1}%i',adr0,nts(k));
    indi = resolve_address(adr);
    nti(k) = indi(4);
end

fname = sprintf('F_%i',conformation+1);
fname2 = sprintf('S_%i',conformation+1);
rnacoor = model.structures{snum}(1).xyz{1};
[msg,rcoor] = get_chain_model([snum 1 1],'rcoor');

[m,~] = size(rnacoor);
rnacoor1 = [rnacoor ones(m,1)];
rnacoor1 = rnacoor1*transmat';
rnacoor1 = rnacoor1(:,1:3);
[m1,~] = size(rcoor);
rcoor1 = [rcoor(:,2:4) ones(m1,1)];
rcoor1 = rcoor1*transmat';
rcoor1 = rcoor1(:,1:3);
% replace coordinates of the rotated bases
[m2,~] = size(dcoor);
poilist = zeros(1,m2);
for k = 1:m2
    diff = sqrt(sum((rnacoor - repmat(dcoor(k,:),m,1)).^2,2));
    [~,poi] = min(diff);
    poilist(k) = poi;
    rnacoor1(poi,:) = dcoor1(k,:);
end
rcoor2 = rcoor1;
ppoi = 0;
for k = 1:m1
    if min(abs(nti-rcoor(k,1))) > 0 % not a base pair in the binding motif
        ppoi = ppoi + 1;
        rcoor2(ppoi,:) = rcoor1(k,:);
    end
end
rcoor2 = rcoor2(1:ppoi,:);
clash = false;
environ = rrmcoor;
options.clash_thr = 1.2;
options.max_shift = 3;

[rcoor2,min_dist,cycles,max_shift] = clash_repair(rcoor2,environ,dcoor1,options);
fprintf(1,'%i cycles for %s: min. dist. %4.2f Å, max. shift %4.2f Å\n',cycles,fname,min_dist,max_shift);
if min_dist + 0.001 < options.clash_thr || max_shift > options.max_shift
    clash = true;
end
% if strcmp(fname,'fit_m1_R_000036.pdb')
%     clash = false;
%     [rnacoor2,min_dist,cycles,max_shift] = clash_repair(rnacoor2,environ,dcoor1,options);
% end
if ~clash
    ppoi = 0;
    for k = 1:m1
        if min(abs(nti-rcoor(k,1))) > 0 % not a base pair in the binding motif
            ppoi = ppoi + 1;
            ccoor = rcoor2(ppoi,:);
            ccoor0 = rcoor(k,2:4);
            diff = sqrt(sum((rnacoor - repmat(ccoor0,m,1)).^2,2));
            [midi,poi] = min(diff);
            if midi < eps
                rnacoor1(poi,:) = ccoor;
            end
        end
    end
    model.structures{snum}(1).xyz{1} = rnacoor1;
    if isfield(model,'selected')
        model = rmfield(model,'selected');
    end
    ctag = chaintag(2);
    model.current_chain = ctag;
    model.chain_tags{snum} = [':' ctag ':'];
    model.structures{snum}.name = ctag;
    % renumber the nts
    id = 1;
    tag = id2tag(id,model.structures{snum}(1).residues{1}.residue_tags);
    newtags = ':';
    while ~isempty(tag)
        oldnum = str2double(tag);
        newnum = oldnum + offset;
        model.structures{snum}(1).residues{1}.info(id).number = newnum;
        newtags = [newtags sprintf('%i',newnum) ':'];
        id = id + 1;
        tag = id2tag(id,model.structures{snum}(1).residues{1}.residue_tags);
    end
    model.structures{snum}(1).residues{1}.residue_tags = newtags;
    selection = ones(1,4);
    if fit_options.optimize
        structure{1} = fit_options.RRM_indices;
        structure{2} = snum;
        selected{1} = snum;
        [~,snum,energy] = optimize_by_tinker('RNA_fit_tinker',structure,selected,fit_options);
        fit_options.energy = energy;
        selection(2) = 2;
        if isnan(energy)
            return
        end
    end
    % select the residues to be stored
    selection(1) = snum;
    select_nt = used(1):used(2);
    poi = 0;
    for k = select_nt
        selection(4) = k;
        poi = poi + 1;
        model.selected{poi} = selection;
    end
    pdbid = sprintf('RNA%s',ctag);
    wr_pdb_selected(fullfile(libdir,fname),pdbid);
    % now write a structure without the binding motif
    if isfield(model,'selected')
        model = rmfield(model,'selected');
    end
    poi = 0;
    for k = select_nt
        if min(abs(nts-k-fit_options.offset)) ~= 0
            selection(4) = k;
            poi = poi + 1;
            model.selected{poi} = selection;
        end
    end
    wr_pdb_selected(fullfile(libdir,fname2),pdbid);
end

function write_unfitted_rna(snum,nts,used,offset,libdir,chaintag,conformation,fit_options)

global model

fname = sprintf('F_%i',conformation+1);
fname2 = sprintf('S_%i',conformation+1);

if isfield(model,'selected')
    model = rmfield(model,'selected');
end
ctag = chaintag(2);
model.current_chain = ctag;
model.chain_tags{snum} = [':' ctag ':'];
model.structures{snum}.name = ctag;
% renumber the nts
id = 1;
tag = id2tag(id,model.structures{snum}(1).residues{1}.residue_tags);
newtags = ':';
while ~isempty(tag)
    oldnum = str2double(tag);
    newnum = oldnum + offset;
    model.structures{snum}(1).residues{1}.info(id).number = newnum;
    newtags = [newtags sprintf('%i',newnum) ':'];
    id = id + 1;
    tag = id2tag(id,model.structures{snum}(1).residues{1}.residue_tags);
end
model.structures{snum}(1).residues{1}.residue_tags = newtags;
selection = ones(1,4);
% select the residues to be stored
selection(1) = snum;
select_nt = used(1):used(2);
poi = 0;
for k = select_nt
    selection(4) = k;
    poi = poi + 1;
    model.selected{poi} = selection;
end
pdbid = sprintf('RNA%s',ctag);
wr_pdb_selected(fullfile(libdir,fname),pdbid);
% now write a structure without the binding motif
if isfield(model,'selected')
    model = rmfield(model,'selected');
end
poi = 0;
for k = select_nt
    if min(abs(nts-k-fit_options.offset)) ~= 0
        selection(4) = k;
        poi = poi + 1;
        model.selected{poi} = selection;
    end
end
wr_pdb_selected(fullfile(libdir,fname2),pdbid);