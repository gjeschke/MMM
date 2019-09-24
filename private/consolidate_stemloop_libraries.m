function consolidate_stemloop_libraries(defs,links,sites)

for lib = 1:length(defs)
    clear linksites
    clear labelsites
    chaintag = defs(lib).motif;
    name = defs(lib).name;
    lpoi = 0;
    for link = 1:length(links)
        poi1 = strfind(links(link).adr1,chaintag);
        if ~isempty(poi1)
            lpoi = lpoi + 1;
            linksites(lpoi).adr = links(link).adr1;
        end
        poi2 = strfind(links(link).adr2,chaintag);
        if ~isempty(poi2)
            lpoi = lpoi + 1;
            linksites(lpoi).adr = links(link).adr2;
        end
    end
    spoi = 0;
    labelsites = [];
    for site = 1:length(sites)
        poi = strfind(sites(site).adr,chaintag);
        if ~isempty(poi)
            spoi = spoi + 1;
            labelsites(spoi).adr = sites(site).adr;
            labelsites(spoi).label = sites(site).label;
        end        
    end
    library = get_stemloop_library(name,chaintag,linksites,labelsites);
    save(name,'library');
end

function library = get_stemloop_library(name,chaintag,linksites,labelsites)

global model
global hMain

mydir = pwd;
cd(name);

RNA_list = dir('*.pdb');

ndecoys = length(RNA_list);

library.chaintag = chaintag;
library.name = name;
library.linksites = linksites;
library.labelsites = labelsites;
library.chains = cell(1,ndecoys);
library.assignment = cell(1,ndecoys);
for k = 1:length(linksites)
    library.linksites(k).coor = zeros(ndecoys,3);
end
for k = 1:length(labelsites)
    library.labelsites(k).coor = zeros(ndecoys,3);
    library.labelsites(k).rmsd = zeros(1,ndecoys);
end

success = 0;
tic,
for kd = 1:ndecoys
    add_msg_board(sprintf('Processing decoy %s.',RNA_list(kd).name));
    [msg,snum] = add_pdb(RNA_list(kd).name);
    if msg.error
        add_msg_board(sprintf('Warning: PDB file %s could not be read (%s). Aborting.',RNA_list(kd).name,msg.text));
        continue
    end
    ctag = chaintag(2);
    model.current_chain = ctag;
    model.chain_tags{snum} = [':' ctag ':'];
    model.structures{snum}.name = ctag;
    success = success + 1;
    library.chains{success} = model.structures{snum}(1);
    library.assignment{success} = RNA_list(kd).name;
    for link = 1:length(linksites)
        [msg,coor] = get_object(sprintf('%s.C5''',linksites(link).adr),'coor');
        if msg.error
            add_msg_board(sprintf('Warning: RNA anchor %s.C5'' coordinate could not be obtained. Aborting.\n',linksites(link).adr));
            success = success - 1;
            continue
        end
        library.linksites(link).coor(success,:) = coor;
    end
    for site = 1:length(labelsites)
        [indices,message] = resolve_address(labelsites(site).adr);
        if message.error
            add_msg_board(sprintf('Warning: Label site %s indices could not be obtained. Aborting.\n',labelsites(site).adr));
            success = success - 1;
            continue
        else
            command=sprintf('rotamers %s %s %i',labelsites(site).adr,labelsites(site).label,298);
            hMain.store_undo=false;
            hMain.dynamic_rotamers=false;
            cmd(hMain,command);
            labels = label_information(model.sites);
            found = false;
            for k = 1:length(labels)
                if abs(sum(indices - labels(k).indices)) == 0
                    found = true;
                    library.labelsites(site).coor(success,:) = labels(k).xyz;
                    library.labelsites(site).rmsd(success) = labels(k).rmsd;
                end
            end
            if ~found
                add_msg_board(sprintf('Warning: Label site %s indices could not be labeled. Aborting.\n',labelsites(site).adr));
                success = success - 1;
                continue
            end
        end
    end
end
toc,

library.chains = library.chains(1:success);
library.assignment = library.assignment(1:success);
for k = 1:length(linksites)
    library.linksites(k).coor = library.linksites(k).coor(1:success,:);
end
for k = 1:length(labelsites)
    library.labelsites(k).coor = library.labelsites(k).coor(1:success,:);
    library.labelsites(k).rmsd = library.labelsites(k).rmsd(1:success);
end

for k = 1:length(library.chains)
    heavy_coor = eliminate_H(library.chains{k}.xyz{1},library.chains{k}.isotopes);
    faces = convhulln(heavy_coor);
    [na,~] = size(faces);
    [kpa,ar] = reducepatch(faces,heavy_coor,na);
    library.hulls(k).vertices = double(ar);
    library.hulls(k).faces = kpa;
end

cd(mydir);

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
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function coor = eliminate_H(coor0,isotopes)

[m,~] = size(coor0);
coor = zeros(m,3);
poi = 0;
for k = 1:m
    if isotopes(k,1) ~= 1
        poi = poi + 1;
        coor(poi,:) = coor0(k,:);
    end
end
coor = coor(1:poi,:);