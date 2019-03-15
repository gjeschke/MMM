function loop_assembler_PTB1(fname)
% Morphs and assigns loop models to rigid-body arrangements of PTBP1 with
% an already linked EMCV-IRES RNA model
% models are tested for fulfilment of small-angle scattering restraints
%
% assumes a set of RBAs in structure PTB1
% assumes a set of RBAs with connected RRM1/RRM2 loop in structure PTF1
% assumes a set of RBAs with connected RRM2/RRM3 loop in structure PTF2
%
% G. Jeschke, 5.2.2019

global model

soln_name = sprintf('%s_linkable.dat',fname);
report_name = sprintf('%s_assembler_report.dat',fname);
list_file = sprintf('%s_scores.dat',fname);
match_file = sprintf('%s_section_matches.mat',fname);

restraint_file = 'PTBP1_restraints_190208.dat';
solutions = load(soln_name);

tic,

snum_RBA = resolve_address('[PTB1]');
snum_f1 = resolve_address('[PTF1]');
excise_chain(snum_f1,1,4,155,181,'L');
snum_f2 = resolve_address('[PTF2]');
excise_chain(snum_f2,2,4,284,336,'L');
snum_new = copy_structure(snum_RBA,'PTBP');

n_RBA = length(model.structures{snum_RBA}(1).xyz);
n_f1 = length(model.structures{snum_f1}(1).xyz);
n_f2 = length(model.structures{snum_f2}(1).xyz);

chi2_SAS = zeros(n_RBA,3);

fid = fopen(report_name,'wt');

fprintf(fid,'--- PTBP1/EMCV-IRES loop assembler report ---\n\n');
fprintf(fid,'%i RBAs with linked RNA\n',n_RBA);
fprintf(fid,'%i RBAs with RMM1/RRM2 loop model\n',n_f1);
fprintf(fid,'%i RBAs with RRM2/RRM3 loop model\n\n',n_f2);

% extract backbone coordinates of RRM1 and RRM2 for link(1,2) models
f1_RRM1 = cell(1,n_f1);
f1_RRM2 = cell(1,n_f1);
for k1 = 1:n_f1
    rind0 = [snum_f1,1,k1];
    % residues = model.structures{snum_f2}(1).residues{k1};
    r0 = 57;
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kr = 58:154 % RRM1
        rind = rind0;
        rind(4) = kr - r0;
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    f1_RRM1{k1} = ch_coor(1:cpoi,:);
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kr = 182:283 % RRM2
        rind = rind0;
        rind(4) = kr - r0;
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    f1_RRM2{k1} = ch_coor(1:cpoi,:);
end

% extract backbone coordinates of RRM2 and RRM34 for link(2,3) models
f2_RRM2 = cell(1,n_f2);
f2_RRM34 = cell(1,n_f2);
for k2 = 1:n_f2
    rind0 = [snum_f2,2,k2];
    % residues = model.structures{snum_f2}(1).residues{k1};
    r0 = 181;
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kr = 182:283 % RRM2
        rind = rind0;
        rind(4) = kr - r0;
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    f2_RRM2{k2} = ch_coor(1:cpoi,:);
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kr = 337:531 % RRM34
        rind = rind0;
        rind(4) = kr - r0;
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    f2_RRM34{k2} = ch_coor(1:cpoi,:);
end

% Determine best-matching peptide-linked RBAs for all RNA-linked RBAs 
matches = zeros(n_RBA,2);
all_matches_f1 = zeros(n_RBA,n_f1);
all_matches_f2 = zeros(n_RBA,n_f2);
rmsds = matches;
stitch_indices = zeros(5,3);
stitch_indices(2,:) = [snum_f1,4,0];
stitch_indices(4,:) = [snum_f2,4,0];
nindices = [snum_new,length(model.structures{snum_new})+1,0];
chain_id = 'P';
init = true;
for kr = 1:n_RBA
    mtag = sprintf('%i',kr);
    nindices(3) = kr;
    cmind0 = [snum_RBA,0,kr];
    % RRM1
    cmind = cmind0;
    cmind(2) = 1;
    % [~,rrm1] = get_chain_model(cmind,'xyz_heavy');
    stitch_indices(1,:) = cmind;
    r0 = 57;
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kres = 58:154 % RRM1
        rind = [cmind kres-r0];
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    rrm1_bb = ch_coor(1:cpoi,:);
    % RRM2
    cmind = cmind0;
    cmind(2) = 2;
    stitch_indices(3,:) = cmind;
    % [~,rrm2] = get_chain_model(cmind,'xyz_heavy');
    r0 = 181;
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kres = 182:283 % RRM2
        rind = [cmind kres-r0];
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    rrm2_bb = ch_coor(1:cpoi,:);
    % RRM34
    cmind = cmind0;
    cmind(2) = 3;
    stitch_indices(5,:) = cmind;
    % [~,rrm34] = get_chain_model(cmind,'xyz_heavy');
    r0 = 336;
    ch_coor = zeros(10000,3);
    cpoi = 0;
    for kres = 337:531 % RRM2
        rind = [cmind kres-r0];
        [~,coor] = get_residue(rind,'xyz_backbone');
        [m,~] = size(coor);
        ch_coor(cpoi+1:cpoi+m,:) = coor;
        cpoi = cpoi + m;
    end
    rrm34_bb = ch_coor(1:cpoi,:);
    min_rmsd = 1e6;
    best_fit =0;
    for k1 = 1:n_f1
        rrm1_f1 = f1_RRM1{k1};
        rrm2_f1 = f1_RRM2{k1};
        coor1 = [rrm1_bb;rrm2_bb];
        coor2 = [rrm1_f1;rrm2_f1];
        rmsd = rmsd_superimpose(coor1,coor2);
        all_matches_f1(kr,k1) = rmsd;
        if rmsd < min_rmsd
            min_rmsd = rmsd;
            best_fit = k1;
        end
    end
    matches(kr,1) = best_fit;
    stitch_indices(2,3) = best_fit;
    rmsds(kr,1) = min_rmsd;
    fprintf(fid,'RBA %i: Best RRM1-RRM2 match %i at %4.2f Å\n',kr,best_fit,min_rmsd);
    mf1tag = sprintf('%i',best_fit);
    min_rmsd = 1e6;
    best_fit =0;
    for k2 = 1:n_f2
        rrm2_f2 = f2_RRM2{k2};
        rrm34_f2 = f2_RRM34{k2};
        coor1 = [rrm2_bb;rrm34_bb];
        coor2 = [rrm2_f2;rrm34_f2];
        rmsd = rmsd_superimpose(coor1,coor2);
        all_matches_f2(kr,k2) = rmsd;
        if rmsd < min_rmsd
            min_rmsd = rmsd;
            best_fit = k2;
        end
    end
    matches(kr,2) = best_fit;
    stitch_indices(4,3) = best_fit;
    rmsds(kr,2) = min_rmsd;
    mf2tag = sprintf('%i',best_fit);
    fprintf(fid,'RBA %i: Best RRM2-RRM34 match %i at %4.2f Å\n',kr,best_fit,min_rmsd);
    stitch_chain(stitch_indices,nindices,chain_id,init);
    replace_section(['[PTBP](P){' mtag '}154'],['[PTBP](P){' mtag '}182'],['[PTF1](A){' mf1tag '}154'],['[PTF1](A){' mf1tag '}182']);
    replace_section(['[PTBP](P){' mtag '}283'],['[PTBP](P){' mtag '}337'],['[PTF2](B){' mf2tag '}283'],['[PTF2](B){' mf2tag '}337']);
    init = false;
    protein = resolve_address(['[PTBP](P){' mtag '}']);
    RNA = resolve_address(['[PTBP](D){' mtag '}']);
    [chi2_SANS,chi2_SAXS] = fit_SAS(protein,RNA,fid);
    chi2_SAS(kr,1:2) = chi2_SANS;
    chi2_SAS(kr,3) = chi2_SAXS;
end

figure(1); clf; hold on;
plot(1:n_RBA,chi2_SAS(:,1),'.','Color',[0,0,0.75]);
plot(1:n_RBA,chi2_SAS(:,1),'-','Color',[0,0,0.75]);
plot(1:n_RBA,chi2_SAS(:,2),'.','Color',[0,0.75,0]);
plot(1:n_RBA,chi2_SAS(:,2),'-','Color',[0,0.75,0]);
plot(1:n_RBA,chi2_SAS(:,3),'.','Color',[0.75,0,0]);
plot(1:n_RBA,chi2_SAS(:,3),'-','Color',[0.75,0,0]);

full_chi2 = sum(chi2_SAS,2);
[sorted_chi2,order] = sort(full_chi2);

figure(2); clf; hold on;
plot(1:n_RBA,sorted_chi2,'k.');
plot(1:n_RBA,sorted_chi2,'k');

fprintf(fid,'\n\n--- Models sorted by SAS chi 2 ---\n\n');

for kr = 1:n_RBA
    fprintf(fid,'Model %i is R%i.%i (%i) at %5.2f: SANS(1.2 m) = %5.2f, SANS(4 m) = %5.2f, SAXS = %5.2f\n',...
        kr,solutions(order(kr),:),order(kr),full_chi2(order(kr)),chi2_SAS(order(kr),:));
end

fclose(fid);

for kr = 1:n_RBA
    poi = order(kr);
    if isfield(model,'selected')
        model = rmfield(model,'selected');
    end
    mtag = sprintf('%i',poi);
    protein = resolve_address(['[PTBP](P){' mtag '}']);
    RNA = resolve_address(['[PTBP](D){' mtag '}']);
    model.selected{1} = protein;
    model.selected{2} = RNA;
    model_name = sprintf('PTB1_f%i_R%i_%i_m%i',poi,solutions(poi,:),kr);
    score = full_chi2(poi);
    wr_pdb_selected(model_name,'PTB1',score);
    PTBP1_quality_check(model_name,restraint_file,list_file);
end

save(match_file,'all_matches_f1','all_matches_f2');

toc,

function excise_chain(snum,cnum,ncnum,resi,rese,ctag)

global model

model.structures{snum}(ncnum).name = ctag;
model.structures{snum}(ncnum).seqtype = model.structures{snum}(cnum).seqtype;
model.structures{snum}(ncnum).modified = 0;
model.structures{snum}(ncnum).nonstandard = 0;
model.structures{snum}(ncnum).helices = 0;
model.structures{snum}(ncnum).strands = 0;
model.structures{snum}(ncnum).helix = [];
model.structures{snum}(ncnum).strand = [];
model.structures{snum}(ncnum).loop_defs = {};
model.structures{snum}(ncnum).helix_defs = {};
model.structures{snum}(ncnum).sheet_defs = {};
model.structures{snum}(ncnum).maxconn = model.structures{snum}(cnum).maxconn;

resnum = model.structures{snum}(cnum).resnum;
resmask = zeros(1,resnum);
xyz = model.structures{snum}(cnum).xyz{1};
[m,~] = size(xyz);
xyzmask = zeros(1,m);

restags = ':';
residue_tags = ':';
for kr = 1:resnum
    info = model.structures{snum}(cnum).residues{1}.info(kr);
    if info.number >= resi && info.number <= rese
        resmask(kr) = 1;
        restags = [restags model.structures{1}(1).residues{1}.info(kr).name ':'];
        residue_tags = [residue_tags sprintf('%i',info.number) ':'];
        for ka = 1:length(model.structures{1}(1).residues{1}.info(kr).atom_numbers)
            anum = model.structures{1}(1).residues{1}.info(kr).atom_numbers{ka};
            [ml,~] = size(anum);
            for kl = 1:ml
                xyzmask(anum(kl,1)) = 1;
            end
        end
    end
end

seq = model.structures{snum}(cnum).sequence;
for k = 1:resnum
    if ~resmask(k)
        seq(k) = '?';
    end
end
model.structures{snum}(ncnum).sequence = seq;
model.structures{snum}(ncnum).restags = restags;
model.structures{snum}(ncnum).resnum = sum(resmask);
isotopes = model.structures{snum}(cnum).isotopes;
isotopes = isotopes(xyzmask>0,:);
model.structures{snum}(ncnum).isotopes = isotopes;
conn = model.structures{snum}(cnum).conn;
conn = conn(xyzmask>0,:);
model.structures{snum}(ncnum).conn = conn;

for km = 1:length(model.structures{snum}(cnum).xyz)
    xyz = model.structures{snum}(cnum).xyz{km};
    xyz = xyz(xyzmask>0,:);
    model.structures{snum}(ncnum).xyz{km} = xyz;
    Bfactor = model.structures{snum}(cnum).Bfactor{km};
    Bfactor = Bfactor(xyzmask>0);
    model.structures{snum}(ncnum).Bfactor{km} = Bfactor;
    Btensor = model.structures{snum}(cnum).Btensor{km};
    Btensor = Btensor(xyzmask>0,:);
    model.structures{snum}(ncnum).Btensor{km} = Btensor;
    model.structures{snum}(ncnum).atoms{km} = sum(xyzmask);
    model.structures{snum}(ncnum).residues{km}.residue_tags = residue_tags;
    model.structures{snum}(ncnum).residues{km}.info = model.structures{snum}(cnum).residues{km}.info(resmask>0);
    for kr = 1:sum(resmask)
        model.structures{snum}(ncnum).residues{km}.info(kr).dssp = ' ';
    end    
end

function [chi2_SANS,chi2_SAXS] = fit_SAS(protein,RNA,fid)

global model

options.err = true;
options.lm = 50;
options.fb = 18;
options.D2O = 0.66;

restraints.SANS(1).data = '1-4sD2O_66_1p2m_atsas.dat';
restraints.SANS(1).illres = 'ill_1p2m.res';
restraints.SANS(2).data = '1-4sD2O_66_4m_atsas.dat';
restraints.SANS(2).illres = 'ill_4m.res';

chi2_SANS = zeros(1,2);
to_be_deleted = '';
for ks = 1:length(restraints.SANS)
    if isfield(model,'selected')
        model = rmfield(model,'selected');
    end
    model.selected{1} = protein;
    model.selected{2} = RNA;
    pdbfile = sprintf('t_%i',ks);
    to_be_deleted = 't*.*';
    wr_pdb_selected(pdbfile,'SANS');
    [chi2,~,~,result] = fit_SANS_by_cryson(restraints.SANS(ks).data,pdbfile,restraints.SANS(ks).illres,options);
    if isempty(chi2) || isnan(chi2)
        chi2_SANS(ks) = 1e6;
            fprintf(2,'Warning: SANS fitting failed\n');
            fprintf(fid,'Warning: SANS fitting of curve %i failed\n',ks);
            fprintf(2,'%s',result);
    else
        fprintf(fid,'SANS curve %i fitted with chi^2 of %6.3f\n',ks,chi2);
        chi2_SANS(ks) = chi2;
    end
end
fprintf(fid,'Sum of SANS curve chi^2: %6.3f\n',sum(chi2_SANS));
delete(to_be_deleted);

if isfield(model,'selected')
    model = rmfield(model,'selected');
end
model.selected{1} = protein;
model.selected{2} = RNA;
pdbfile = sprintf('t_%i',ks);
to_be_deleted = 't*.*';
wr_pdb_selected(pdbfile,'SAXS');
SAXS_curve = load_SAXS_curve('1-4s-Buffer.dat');
options.sm = 10*max(SAXS_curve(:,1));
[chi2_SAXS,~,~,result] = fit_SAXS_by_crysol('1-4s-Buffer.dat',pdbfile,options);
if isempty(chi2_SAXS) || isnan(chi2_SAXS)
    fprintf(2,'Warning: SAXS fitting failed\n');
    fprintf(fid,'Warning: SAXS fitting of curve %i failed\n',ks);
    fprintf(2,'%s',result);
else
    fprintf(fid,'SAXS curve fitted with chi^2 of %6.3f\n',chi2_SAXS);
end
delete(to_be_deleted);
    