function make_test_restraints_ENM

fname = '1G6S_20_distances_TM_35.dat';
ostub = '1G6S';

[~,diagnostics0] = bilabel_site_scan([1 1],'IDA');
[~,diagnostics,sites,Cu_coor] = bilabel_site_scan([2 1],'IDA');

fid=fopen(fname);
if fid==-1,
    add_msg_board('ERROR: Restraint file does not exist');
    return;
end;

ofid_dir = fopen([ostub '_direct.dat'],'wt');
ofid_MTSL = fopen([ostub '_MTSL.dat'],'wt');
ofid_IDA = fopen([ostub '_IDA.dat'],'wt');

nl = 0;
mode = 0;
while 1
    tline = fgetl(fid);
    nl = nl + 1;
    if ~ischar(tline) || mode<0, break, end
    if ~isempty(tline),
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k),
            if k(1)>1,
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end;
        end;
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#'),
            if strcmpi(char(args(2)),'DEER')
                mode = 1;
                fprintf(ofid_dir,'# DIRECT\n');
                fprintf(ofid_MTSL,'# DEER MTSL 298\n');
                fprintf(ofid_IDA,'# ODEER IDA IDA 298\n');
                continue
            else
                mode = 0;
            end
        end
        if mode == 1
            adr1 = char(args(1));
            adr2 = char(args(2));
            [~,xyz1] = get_object(strcat(adr1,'.CA'),'xyz');
            [~,xyz2] = get_object(strcat(adr2,'.CA'),'xyz');
            r = norm(xyz1-xyz2)/10;
            fprintf(ofid_dir,'%8s%8s%8.2f%8.2f\n',adr1,adr2,r,0.03);
            [rmean,rmsd] = get_MTSL_distance(adr1,adr2);
            fprintf(ofid_MTSL,'%8s%8s%8.2f%8.2f\n',adr1,adr2,rmean,rmsd);
            [rmean,rmsd,adr1n,adr2n,adr3n,adr4n] = get_IDA_distance(adr1,adr2,diagnostics,sites,Cu_coor,diagnostics0);
            if isempty(rmean)
                fprintf(ofid_IDA,'ERROR: no Cu label for site pair %s (%s,%s), %s (%s,%s)\n',adr1,adr1n,adr2n,adr2,adr3n,adr4n);
            else
                fprintf(ofid_IDA,'%8s%8s%8s%8s%8.2f%8.2f\n',adr1n,adr2n,adr3n,adr4n,rmean,rmsd);
            end
        end
    end
end

fclose(fid);
fclose(ofid_dir);
fclose(ofid_MTSL);
fclose(ofid_IDA);

function [rmean,rmsd] = get_MTSL_distance(adr1,adr2)

global rotamer_libraries

label='MTSL';
rmean = 0;
rmsd = 0;

rindices1 = resolve_address(adr1);
[m1,n1] = size(rindices1);
if m1 ~= 1 || n1 ~=4,
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr1));
    return;
end;
rindices2 = resolve_address(adr2);
[m2,n2] = size(rindices2);
if m2 ~= 1 || n2 ~=4,
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr2));
    return
end;

rot_lib_name = '';

for k = 1:length(rotamer_libraries),
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label),
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T),
            if rotamer_libraries(k).T(kk) == 298,
                id = kk;
            end;
        end;
        if id >0,
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
        end;
    end;
end;

load(rot_lib_name);
label1 = rotamer_populations(rindices1,rot_lib);
label2 = rotamer_populations(rindices2,rot_lib);
[rax,distr]=get_distribution(label1(1).NOpos,label2(1).NOpos,0.1);
distr = distr/sum(distr);
rmean = sum(rax.*distr);
add_msg_board(sprintf('Mean distance : %5.1f Å',10*rmean));
rdev = rax - rmean;
mom2 = sum(rdev.^2.*distr);
rmsd = sqrt(mom2);

function [rmean,rmsd,adr1,adr2,adr3,adr4] = get_IDA_distance(adr1,adr2,diagnostics,sites,Cu_coor,diagnostics0)

rmean = [];
rmsd = [];
ind1 = resolve_address(adr1);
ind2 = resolve_address(adr2);
[adr1,adr2,popcoor1] = dHis_site(ind1,diagnostics,sites,Cu_coor,diagnostics0);
[adr3,adr4,popcoor2] = dHis_site(ind2,diagnostics,sites,Cu_coor,diagnostics0);
if isempty(popcoor1) || isempty(popcoor2)
    return
end
[rax,distr]=get_distribution(popcoor1,popcoor2,0.1);
distr = distr/sum(distr);
rmean = sum(rax.*distr);
add_msg_board(sprintf('Mean distance : %5.1f Å',10*rmean));
rdev = rax - rmean;
mom2 = sum(rdev.^2.*distr);
rmsd = sqrt(mom2);
indices = resolve_address(adr1);
[~,ctag,~,resnum] = mk_address_parts(indices);
adr1 = sprintf('(%s)%i',ctag,resnum);
indices = resolve_address(adr2);
[~,ctag,~,resnum] = mk_address_parts(indices);
adr2 = sprintf('(%s)%i',ctag,resnum);
indices = resolve_address(adr3);
[~,ctag,~,resnum] = mk_address_parts(indices);
adr3 = sprintf('(%s)%i',ctag,resnum);
indices = resolve_address(adr4);
[~,ctag,~,resnum] = mk_address_parts(indices);
adr4 = sprintf('(%s)%i',ctag,resnum);

function [adr1,adr2,popcoor] = dHis_site(ind,diagnostics,sites,Cu_coor,diagnostics0)

adr1 = '';
adr2 = '';

[md0,~] = size(diagnostics0);
[md,~] = size(diagnostics);

success = false;
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 2 % most preferred labelling
        if ind(4) == diagnostics(k,7)+1
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 4 % second best labelling
        if ind(4) == diagnostics(k,7)+2
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 1 % third best labelling
        if ind(4) == diagnostics(k,7) || ind(4) == diagnostics(k,8)
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 3 % fourth best labelling
        if ind(4) == diagnostics(k,7)+1 || ind(4) == diagnostics(k,7) + 2 
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 2 
        if ind(4) == diagnostics(k,7) || ind(4) == diagnostics(k,8)
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 4 
        if ind(4) == diagnostics(k,7)+1 || ind(4) == diagnostics(k,7)+3
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 5 
        if ind(4) == diagnostics(k,7)+2 || ind(4) == diagnostics(k,7)+3
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 3 
        if ind(4) == diagnostics(k,7) || ind(4) == diagnostics(k,8)
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 1 
        if ind(4) == diagnostics(k,7)-1 || ind(4) == diagnostics(k,8)+1
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 4 
        if ind(4) == diagnostics(k,7) || ind(4) == diagnostics(k,8)
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 5 
        if ind(4) == diagnostics(k,7)+1 || ind(4) == diagnostics(k,7)+4
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 5 
        if ind(4) == diagnostics(k,7) || ind(4) == diagnostics(k,8)
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 4 
        if ind(4) == diagnostics(k,7)-1 || ind(4) == diagnostics(k,8)+1
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 2 
        if ind(4) == diagnostics(k,7)-1 || ind(4) == diagnostics(k,8)+1
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 5 
        if ind(4) == diagnostics(k,7)-1 || ind(4) == diagnostics(k,8)+1
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 1 
        if ind(4) == diagnostics(k,7)-2 || ind(4) == diagnostics(k,8)+2
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 2 
        if ind(4) == diagnostics(k,7)-2 || ind(4) == diagnostics(k,8)+2
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if success
    return
end
for k = 1:md
    % check whether this dHis can also be labelled for the initial
    % structure
    matcher = repmat(diagnostics(k,7:8),md0,1);
    diff = sum(abs(diagnostics0(:,7:8)-matcher),2);
    if min(diff) > 0
        continue
    end
    if diagnostics(k,6) == 4 
        if ind(4) == diagnostics(k,7)-2 || ind(4) == diagnostics(k,8)+2
            popcoor = Cu_coor{k};
            dadr = sites{k};
            spoi = strfind(dadr,'|');
            adr1 = dadr(1:spoi-1);
            adr2 = dadr(spoi+1:end);
            success = true;
            break
        end
    end
end
if ~success
    popcoor = [];
end
