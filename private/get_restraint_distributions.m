function [rax,distributions] = get_restraint_distributions(restraints)

for k = 1:length(restraints.DEER)
    ind1 = resolve_address(restraints.DEER(k).adr1);
    NO_pos1 = get_NO_pos(ind1,restraints.DEER(k).label,298);
    ind2 = resolve_address(restraints.DEER(k).adr2);
    NO_pos2 = get_NO_pos(ind2,restraints.DEER(k).label,298);
    [rax,sim_distr] = get_distribution(NO_pos1,NO_pos2,0.05);
    if k == 1
        distributions = zeros(length(restraints.DEER),length(rax));
        distributions(1,:) = sim_distr;
    else
        distributions(k,:) = sim_distr;
    end
end


function NO_pos = get_NO_pos(indices,label,T)

global model
global label_defs
global hMain

if strcmpi(label,'CA')
    adr = sprintf('%s.CA',mk_address(indices));
    [~,xyz] = get_object(adr,'coor');
    NO_pos = [xyz 1];
    return
end

NO_pos = [];
if isfield(model,'sites')
    for k0=1:length(model.sites)
        for k1=1:length(model.sites{k0})
            for k=1:length(model.sites{k0}(k1).residue)
                if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0
                    id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                    if strcmpi(label,label_defs.residues(id).short_name)  || strcmpi(label,label_defs.residues(id).tc)
                        if T == model.sites{k0}(k1).residue(k).T
                            NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                        end
                    end
                end
            end
        end
    end
end

if isempty(NO_pos)
    adr = mk_address(indices);
    command=sprintf('rotamers %s %s %i',adr,label,T);
    hMain.store_undo=false;
    hMain.dynamic_rotamers=false;
    cmd(hMain,command);
end

for k0=1:length(model.sites)
    for k1=1:length(model.sites{k0})
        for k=1:length(model.sites{k0}(k1).residue)
            if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0
                id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                if strcmpi(label,label_defs.residues(id).short_name) || strcmpi(label,label_defs.residues(id).tc)
                    if T == model.sites{k0}(k1).residue(k).T
                        NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                    end
                end
            end
        end
    end
end

