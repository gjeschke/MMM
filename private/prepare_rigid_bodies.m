function prepare_rigid_bodies(fname)
% Given a draft structure with chains belonging to a set of rigid bodies
% and a RigiFlex restraint file, a PDB input file for RigiFlex is written
% where individual rigid bodies are well separated
%
% fname     file name of restraint file (with extension)
%
% output file name is determined by the # PDB command in the restraint
% file
%
% G. Jeschke, 10.9.2019

global model

clear restraints

restraints.PDB = 'RBD0';
restraints.rb(1).chains = [];
restraints.rb(1).chainid = cell(1,10);

rb_poi=0;

% Read restraint file to determine chain assignment and output filename

fid=fopen(fname);
if fid==-1
    add_msg_board('ERROR: Restraint file does not exist');
    return;
end

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            if k(1)>1
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end
        end
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#')
            switch upper(char(args(2)))
                case 'PDB'
                    pdbid = strtok(char(args(3)),':');
                    restraints.PDB = pdbid; % to conform to old format
                case 'RIGID'
                    if length(args)<= 2
                        add_msg_board(sprintf('ERROR: Rigid-body definition in line %i misses chain identifiers. Aborting.',nl));
                        fclose(fid);
                        return
                    end
                    rb_poi = rb_poi + 1;
                    chains = zeros(1,length(args)-2);
                    for kc = 1:length(args)-2
                        indices = resolve_address(char(args(kc+2)));
                        restraints.rb(rb_poi).chainid{kc} = char(args(kc+2));
                        if length(indices) ~= 2
                            add_msg_board(sprintf('ERROR: Rigid-body definition contains a wrong chain address. Aborting.',nl));
                            fclose(fid);
                            return
                        end
                        snum = indices(1);
                        chains(kc) = indices(2);
                    end
                    restraints.rb(rb_poi).chains  = chains;
            end
        end
    end
end

fclose(fid);

maxr = 0;
centers = zeros(rb_poi,3);  % center coordinates of each rigid body

for krb = 1:rb_poi
    atnum = 0;
    bcenter = zeros(1,3);
    for kc = 1:length(restraints.rb(krb).chains)
        [msg,xyz] = get_chain([snum restraints.rb(krb).chains(kc)],'xyz');
        if msg.error
            add_msg_board(sprintf('ERROR: Chain %s does not exist in current model',restraints.rb(krb).chainid{kc}));
            return
        else
            [m,~] = size(xyz);
            atnum = atnum + m;
            bcenter = bcenter + m*mean(xyz);
        end
    end
    centers(krb,:) = bcenter/atnum;
end

for krb = 1:rb_poi
    for kc = 1:length(restraints.rb(krb).chains)
        [msg,xyz] = get_chain([snum restraints.rb(krb).chains(kc)],'xyz');
        if msg.error
            add_msg_board(sprintf('ERROR: Chain %s does not exist in current model',restraints.rb(krb).chainid{kc}));
            return
        else
            [m,~] = size(xyz);
            xyzs = xyz - repmat(centers(krb,:),m,1);
            rmax = max(sqrt(sum(xyzs.^2,2)));
            if rmax > maxr
                maxr = rmax;
            end
        end
    end
end

[midpoints,mindist] = fibonacci_sphere(rb_poi);
midpoints = midpoints*2*maxr/mindist;

for krb = 1:rb_poi
    shift = midpoints(krb,:) - centers(krb,:);
    for kc = 1:length(restraints.rb(krb).chains)
        for km = 1: length(model.structures{snum}(restraints.rb(krb).chains(kc)).xyz)
            xyz = model.structures{snum}(restraints.rb(krb).chains(kc)).xyz{km};
            [m,~] = size(xyz);
            xyz = xyz + repmat(shift,m,1);
            model.structures{snum}(restraints.rb(krb).chains(kc)).xyz{km} = xyz;
        end
    end
end

if isfield(model,'selected')
    model = rmfield(model,'selected');
end
model.selected{1} = snum;
message = wr_pdb_selected(restraints.PDB,restraints.PDB);
if message.error
    add_msg_board(sprintf('ERROR: PDB file %s.pdb could not be written.',restraints.PDB));
end