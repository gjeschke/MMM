function adjust_section(indices,anchors,v)
% adjust_section(indices,anchors,v)
%
% - adjusts coordinates of flexible residue ranges to rigid-body
%   transformations
% - behavior is defined only if the sections are continuous
% - it is recommended to optimize the chain or segment afterwards
%
% indices   initial residue address, N-terminal anchor, where the sections 
%           are superimposed
% anchors   anchor points, for which rigid-body transformations are defined
% v         vector of translation and Euler angle parameters
%
% G. Jeschke, 1.11.2018

global model

indi = indices(1:4);
indf = indices(5:8);

inconsistent = sum(abs(indi(1:3)-indf(1:3)));

if inconsistent
    add_msg_board(sprintf('ERROR: Target initial (%s) and final (%s) residue address do not belong to the same chain model in section adjustment. Aborting.',res1i,res1f));
    return
end

rb = round(length(v)/6); % number of rigid bodies
euler = zeros(rb,3);
trans = zeros(rb,3);

for kr = 1:rb
    baspoi6 = 6*(kr-1);
    trans(kr,:) = v(baspoi6+1:baspoi6+3); % reconvert to Angstroem
    euler0 = v(baspoi6+4:baspoi6+6);
    if max(abs(euler0)) > 0.05
        euler_corr = 0.05*euler0/max(abs(euler0));
    else
        euler_corr = euler0;
    end
    euler(kr,:) = euler_corr;
end

target_ind = indi;
for k = indi(4):indf(4) % loop over all residues that have to be adjusted
    target_ind(4) = k; 
    atom_numbers = model.structures{target_ind(1)}(target_ind(2)).residues{target_ind(3)}.info(target_ind(4)).atom_numbers;
    nat = length(atom_numbers);
    xyz = model.structures{target_ind(1)}(target_ind(2)).xyz{target_ind(3)};
    for kat = 1:nat
        targ_anums = atom_numbers{kat};
        [mtarg,~] = size(targ_anums);
        for kl = 1:mtarg
            poi_targ = targ_anums(kl,1);
            coor0 = xyz(poi_targ,:);
            coor = adjust_xyz(coor0,anchors,trans,euler);
            xyz(poi_targ,:) = coor;
        end
    end
    model.structures{target_ind(1)}(target_ind(2)).xyz{target_ind(3)} = xyz;
end


function coor = adjust_xyz(coor0,anchors,trans,euler)

weights = zeros(1,length(anchors));
for ka = 1:length(anchors)
    d = norm(coor0 - anchors(ka).coor) + 1e-6; % avoid accidental zero coordinate
    weights(ka) = 1/d;
end
weights = weights/sum(weights); % normalization
mean_trans = zeros(1,3);
mean_euler = zeros(1,3);
for ka = 1:length(anchors)
    rb = anchors(ka).rb;
    if rb > 1 % the first rigid body has no translation and no rotation
        mean_trans = mean_trans + weights(ka)*trans(rb-1,:);
        mean_euler = mean_euler + weights(ka)*euler(rb-1,:);
    end
end

transmat =  transrot2affine(mean_trans,mean_euler);
coor = affine_trafo_coor(coor0,transmat);