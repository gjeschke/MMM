function [popcoor,rotamers,part_fct,mean_Cu_coor,rmsd] = two_pronged_label(adr1,adr2,label,ecoor1,torsionpot)

global model
global ligand_libraries
global label_defs

popcoor = [];
rotamers = [];
part_fct = 0;
mean_Cu_coor = [];
rmsd = [];

gas_un = 8.314472;    % universal gas constant in CI (J/(mol*K)       
T = 298;
run_opt.forgive = 0.5;

if ~exist('torsionpot','var')
    torsionpot = [];
end

for k = 1:length(ligand_libraries)
    if strcmpi(label,ligand_libraries(k).label) || strcmpi(label,ligand_libraries(k).tc)
        ligand_def = ligand_libraries(k);
    end
end

if ~exist('ligand_def','var')
    add_msg_board(sprintf('ERROR: Two-pronged label %s is not registered with MMM',label));
    return
end

if ~isempty(torsionpot)
    ligand_def.V_ksi = torsionpot;
end

res1 = id2tag(1,ligand_def.res);
res2 = id2tag(2,ligand_def.res);

ind1 = resolve_address(adr1);
msg1 = compute_rotamers(ind1,res1,'ambient');
if msg1.error
    add_msg_board(sprintf('ERROR in rotamer computation for attachment residue %s: %s',res1,msg1.text));
    return
end

ind2 = resolve_address(adr2);
msg2 = compute_rotamers(ind2,res2,'ambient');
if msg2.error
    add_msg_board(sprintf('ERROR in rotamer computation for attachment residue %s: %s',res2,msg2.text));
    return
end

for k0 = 1:length(model.sites)
    for k1 = 1:length(model.sites{k0})
        for k = 1:length(model.sites{k0}(k1).residue)
            indices = model.sites{k0}(k1).residue(k).indices;
            id = tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
            tag = label_defs.residues(id).short_name;
            tc = label_defs.residues(id).tc;
            if length(indices) == length(ind1) && sum(abs(indices-ind1)) == 0
                if strcmpi(tag,res1) || strcmpi(tc,res1)
                    pos1 = model.sites{k0}(k1).residue(k).NOpos;
                end
            end
            if length(indices) == length(ind2) && sum(abs(indices-ind2)) == 0
                if strcmpi(tag,res2) || strcmpi(tc,res2)
                    pos2 = model.sites{k0}(k1).residue(k).NOpos;
                end
            end
        end
    end
end

if ~exist('pos1','var')
    add_msg_board(sprintf('ERROR: Residue %s could not be labeled with %s',adr1,res1));
    return
end

if ~exist('pos2','var')
    add_msg_board(sprintf('ERROR: Residue %s could not be labeled with %s',adr2,res2));
    return
end

[m1,~] = size(pos1);
[m2,~] = size(pos2);

dist = zeros(m1*m2,4);
poi = 0;
for k1 = 1:m1
    for k2 = 1:m2
        cdist = norm(pos1(k1,1:3)-pos2(k2,1:3));
        % fprintf(1,'His.N3-His.N3'' distance: %5.2f Å\n',cdist);
        lib = 0;
        for k3 = 1:length(ligand_def.dist)
            if cdist > ligand_def.dist(k3)
                lib = lib + 1;
            end
        end
        if lib > 0 && lib < length(ligand_def.dist)
            poi = poi + 1;
            dist(poi,1) = cdist;
            dist(poi,2) = k1;
            dist(poi,3) = k2;
            dist(poi,4) = lib;
        end
    end
end
dist = dist(1:poi,:);

% colors = [0,0,0.6;0.75,0,0];

[msg,coor]=get_object(ind1(1:2),'xyz');

if msg.error
    add_msg_board(sprintf('ERROR: Center of mass of the chain could not be determined: %s',mag.text));
    return
end

if iscell(coor)
    coor=coor{1};
end

cofm = mean(coor,1);

liblen = 72;
poten_ext = zeros(poi*liblen,1);
poten_full = zeros(poi*liblen,1);
pos_Cu = zeros(poi*liblen,3);
rotamers = zeros(poi*liblen,3);

for k = 1:poi
    basind = (k-1)*liblen;
    N3 = pos1(dist(k,2),1:3);
    N3p = pos2(dist(k,3),1:3);
    orig = (N3+N3p)/2;
    x = N3p-N3;
    x = x/norm(x);
    yp = orig - cofm;
    yp = yp/norm(yp);
    z = cross_rowvec(x,yp);
    z = z/norm(z);
    y = cross_rowvec(z,x);
    Rp = [x;y;z];
    libname = id2tag(dist(poi,4),ligand_def.files);
    load(libname);
    [ml,~] = size(ligand_lib.library(1).ecoor);
    for r = 1:liblen
        rotamers(basind+r,:) = [dist(k,2:3) r];
        ksi = pi*ligand_lib.library(r).dihedrals/180;
        int_en = ligand_def.V_ksi*(1-cos(ksi))/2;
        ecoor2 = ligand_lib.library(r).ecoor;
        coor2 = ecoor2(:,2:4)*Rp;
        coor2 = coor2 + repmat(orig,ml,1);
        ecoor2(:,2:4) = coor2;
        pos_Cu(basind+r,:) = coor2(1,:);
        out = get_energyLJ_UFF(ecoor1,ecoor2,run_opt);
        ext_en = out.energy;
        poten_ext(basind+r)=exp(-ext_en/(gas_un*T));  % Boltzmann factor for LJ potential for k-th rotamer
        poten_full(basind+r)=exp(-(int_en+ext_en)/(gas_un*T));
    end
%     poten_ext = poten_ext/sum(poten_ext);
%     poten_full = poten_full/sum(poten_full);
end
[mean_Cu_coor,rmsd,popcoor,part_fct,rotamers] = Cupos_rmsd(poten_full,pos_Cu,rotamers);
if isfield(model,'bisites')
    exsites = length(model.bisites);
else
    exsites = 0;
end
[nrot,~] = size(popcoor);
if nrot > 0
    add_msg_board(sprintf('Site pair %s,%s has %i rotamers of %s with normalized partition function %5.3f',adr1,adr2,nrot,ligand_def.tc,part_fct));
    model.bisites{exsites+1}.mutants{1} = res1;
    model.bisites{exsites+1}.mutants{2} = res2;
    model.bisites{exsites+1}.library = libname;
    model.bisites{exsites+1}.sites{1}.label = ligand_def.tc;
    model.bisites{exsites+1}.sites{1}.indices = [ind1;ind2];
    model.bisites{exsites+1}.sites{1}.coor = [popcoor rotamers];
    model.bisites{exsites+1}.sites{1}.mean_coor = mean_Cu_coor;
    model.bisites{exsites+1}.sites{1}.rmsd = rmsd;
    model.bisites{exsites+1}.sites{1}.partition_function = part_fct;
else
    add_msg_board(sprintf('Site pair %s,%s could not be labelled with %s',adr1,adr2,ligand_def.tc));
end

function [mean_Cu_coor,rmsd,popcoor,part_fct,rotamers] = Cupos_rmsd(pop,coor,rotamers)
% in nm(!)

part_fct = sum(pop)/length(pop);
popcoor = coor(pop >= 0.01,:);
rotamers = rotamers(pop >= 0.01,:);
pop = pop(pop > 0.01);
pop=pop/sum(pop);
xmean=sum(popcoor(:,1).*pop);
ymean=sum(popcoor(:,2).*pop);
zmean=sum(popcoor(:,3).*pop);
mean_Cu_coor = [xmean,ymean,zmean];
dx=(popcoor(:,1)-xmean);
dy=(popcoor(:,2)-ymean);
dz=(popcoor(:,3)-zmean);
nCu=length(dx);
rmsd=sqrt(0.005+nCu*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nCu-1)); 
popcoor = [popcoor pop];