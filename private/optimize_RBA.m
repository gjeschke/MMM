function [snum,distances] = optimize_RBA(distances,rigid,randomize)
% [snum,distance_restraints] = optimize_RBA(distance_restraints,rigid)
%
% optimizes a rigid-body arrangement with respect to distance restraints
% flexible linkers are adapted by distributing relative translation and
% rotation between the anchor residues evenly
% restraints are culled to involve only those ones that involve at least
% one site in a flexible linker
%
% distances     distance restraints
% rigid         rigid-body definitions by residue ranges
% randomize     flag indicating that randomized distances are to be used
%               defualts to false
%
% snum          structure number of the optimized model
% distances     culled distance restraints
%
% G. Jeschke, 1.11.2018

global model

if ~exist('randomize','var') || isempty(randomize)
    randomize = false;
end

snum = model.current_structure;

if isempty(rigid) || isempty(rigid(1).indices)
    return
end

% assign distance restraints to rigid bodies and cull them
% make restraint set for RBA optimization
dpoi = 0;
spoi = 0;
maxrba = 0;
for kr = 1:length(distances)
    rb1 = assign_rigid_body(distances(kr).indices1,rigid);
    rb2 = assign_rigid_body(distances(kr).indices2,rigid);
    if rb1 > maxrba
        maxrba = rb1;
    end
    if rb2 > maxrba
        maxrba = rb2;
    end
    if rb1 == 0 || rb2 == 0 % restraint involves at least one flexible linker site
        dpoi = dpoi + 1;
        distances(dpoi) = distances(kr);
    else
        spoi = spoi + 1;
        sites(spoi).rb1 = rb1;
        sites(spoi).rb2 = rb2;
        sites(spoi).coor1 = distances(kr).ecoor1(1:3);
        sites(spoi).coor2 = distances(kr).ecoor2(1:3);
        if randomize
            sites(spoi).r = distances(kr).r_rand;
%            sites(spoi).sigr = distances(kr).sigr_rand;        
%             sites(spoi).r = distances(kr).r;
            sites(spoi).sigr = distances(kr).sigr;
        else
            sites(spoi).r = distances(kr).r;
            sites(spoi).sigr = distances(kr).sigr;
        end
    end
end

distances = distances(1:dpoi);

% get the affine transformation matrices for the optimal RBA
if maxrba < 2
    return
end

v = get_optimal_RBA(sites,maxrba);

rb = round(length(v)/6); % number of rigid bodies
transmats = cell(1,rb);

for kr = 1:rb
    baspoi6 = 6*(kr-1);
    trans = v(baspoi6+1:baspoi6+3); % reconvert to Angstroem
    euler0 = v(baspoi6+4:baspoi6+6);
    if max(abs(euler0)) > 0.05
        euler_corr = 0.05*euler0/max(abs(euler0));
    else
        euler_corr = euler0;
    end
    transmats{kr} = transrot2affine(trans,euler_corr);
end

snum0 = snum;

snum = copy_structure(snum0,'+000');
for kr = 2:length(rigid)
    transmat = transmats{kr-1};
    for kp = 1:length(rigid(kr).indices)
        indices = rigid(kr).indices{kp};
        [m,~] = size(indices);
        for kres = 1:m
            cind = indices(kres,:);
            cind(1) = snum;
            set_residue(cind,'transform',transmat);
        end
    end
end

% identify flexible sections and adapt them
chains = length(model.structures{snum0});
flex_ranges = zeros(100,8);
fpoi = 0;
apoi = 0;
for kc = 1:chains
    nres = length(model.structures{snum0}(kc).residues{1}.info);
    flex = false;
    cind0 = [snum0 kc 1 1];
    rb0 = 0;
    for kr = 1:nres
        cind = [snum0 kc 1 kr];
        rb = assign_rigid_body(cind,rigid);
        if ~flex && rb == 0
            flex = true;
            fpoi = fpoi + 1;
            flex_ranges(fpoi,1:4) = cind;
            flex_ranges(fpoi,1) = snum;
            if rb0 ~= 0
                apoi = apoi + 1;
                anchors(apoi).rb = rb0;
                anchors(apoi).indices = cind0;
                [~,xyz] = get_residue(cind0,'xyz_paradigm');
                anchors(apoi).coor = mean(xyz);
            end
        end
        if flex && rb ~= 0
            flex = false;
            flex_ranges(fpoi,5:8) = cind;
            flex_ranges(fpoi,5) = snum;
            flex_ranges(fpoi,8) = flex_ranges(fpoi,8) - 1;
            apoi = apoi + 1;
            anchors(apoi).rb = rb;
            anchors(apoi).indices = cind;
            [~,xyz] = get_residue(cind,'xyz_paradigm');
            anchors(apoi).coor = mean(xyz);
        end
        if rb ~= 0
            rb0 = rb;
        end
        cind0 = cind;
    end        
end
anchors = anchors(1:apoi);

for kf = 1:fpoi
    if isempty(find(flex_ranges(kf,:) == 0,1)) % valid index range
        adjust_section(flex_ranges(kf,:),anchors,v);
    end
end



function rb = assign_rigid_body(residue,rigid)

rb = 0;
for kr = 1:length(rigid)
    for kp = 1:length(rigid(kr).indices)
        indices = rigid(kr).indices{kp};
        [m,~] = size(indices);
        det = sum(abs(indices - repmat(residue,m,1)),2);
        if ~isempty(find(det == 0,1))
            rb = kr;
            break
        end
    end
    if rb ~= 0
        break
    end
end

function v = get_optimal_RBA(sites,maxrba)

v0 = zeros(1,6*(maxrba-1)); % first rigid body remains in place
cost = rba_cost_fct(v0,sites);
fprintf(1,'Initial cost  : %8.4f\n',cost);

options = optimset('Display','none','TolFun',0.001,'TolX',0.03,'MaxFunEvals',50000);
[v,cost] = fminsearch(@rba_cost_fct,v0,options,sites);
fprintf(1,'Optimized cost: %8.4f\n',cost);
fprintf(1,'Translation of RBA2: (%4.1f, %4.1f, %4.1f) Å\n',v(1:3));
fprintf(1,'Translation of RBA3: (%4.1f, %4.1f, %4.1f) Å\n',v(7:9));
fprintf(1,'Rotation of RBA2   : (%4.1f, %4.1f, %4.1f)°\n',180*v(4:6)/pi);
fprintf(1,'Rotation of RBA3   : (%4.1f, %4.1f, %4.1f)°\n',180*v(10:12)/pi);



function cost = rba_cost_fct(v,sites)

rb = round(length(v)/6); % number of rigid bodies
transmats = cell(1,rb);

for kr = 1:rb
    baspoi6 = 6*(kr-1);
    trans = v(baspoi6+1:baspoi6+3); % reconvert to Angstroem
    euler = v(baspoi6+4:baspoi6+6);
    transmats{kr} = transrot2affine(trans,euler);
end

cost = 0;
for ks = 1:length(sites)
    coor1 = sites(ks).coor1;
    if sites(ks).rb1 ~= 1
        coor1 = affine_trafo_point(coor1,transmats{sites(ks).rb1 - 1});
    end
    coor2 = sites(ks).coor2;
    if sites(ks).rb2 ~= 1
        coor2 = affine_trafo_point(coor2,transmats{sites(ks).rb2 - 1});
    end
    r0 = norm(coor1-coor2);
    % cost = cost + sqrt((r0-sites(ks).r)^2);
    cost = cost + sqrt((r0-sites(ks).r)^2/sites(ks).sigr^2);
end