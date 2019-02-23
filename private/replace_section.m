function [trans,euler,anchors] = replace_section(res1i,res1f,res2i,res2f,test)
% [trans,euler,anchor] = replace_section(adr1,res1a,res1e,adr2,res2a,res2e,test)
%
% - replaces a section of amino acid residues in a chain model by
%   corresponding residues from another chain model
% - behavior is defined only if the sections are continuous
% - anchor residues are not replaced but used only for superposition
% - translation and rotation between the two optimal superpositions are
%   uniformly distributed over the segment
% - it is recommended to optimize the chain or segment afterwards
%
% res1i     initial residue address, N-terminal anchor, where the sections 
%           are superimposed
% res1f     final residue adddress, C-terminal anchor for superposition
% res2i     initial residue address of the new conformation
% res2f     final residue address of the new conformation
% test      optional flag for test mode, defaults to false
%
% trans     translation vectors for initial and final anchor residues
% euler     euler rotations for initial and final anchor residues
% anchors   mean backbone coordinates of initial and final anchor residues
%
% G. Jeschke, 18.10.2018

global model

if ~exist('test','var')
    test = false;
end

[indr1i,msg] = resolve_address(res1i);
if msg.error
    add_msg_board(sprintf('ERROR: Unknown target initial residue %s in section replacement (%s). Aborting.',res1i,msg.text));
    return
end
    
if length(indr1i) ~= 4
    add_msg_board(sprintf('ERROR: Target initial residue addresa %s does not specify a residue in section replacement. Aborting.',res1i));
    return
end

[indr1f,msg] = resolve_address(res1f);
if msg.error
    add_msg_board(sprintf('ERROR: Unknown target final residue %s in section replacement (%s). Aborting.',res1f,msg.text));
    return
end
if length(indr1f) ~= 4
    add_msg_board(sprintf('ERROR: Target final residue addresa %s does not specify a residue in section replacement. Aborting.',res1i));
    return
end

[indr2i,msg] = resolve_address(res2i);
if msg.error
    add_msg_board(sprintf('ERROR: Unknown template initial residue %s in section replacement (%s). Aborting.',res2f,msg.text));
    return
end
if length(indr2i) ~= 4
    add_msg_board(sprintf('ERROR: Template initial residue addresa %s does not specify a residue in section replacement. Aborting.',res2i));
    return
end

[indr2f,msg] = resolve_address(res2f);
if msg.error
    add_msg_board(sprintf('ERROR: Unknown template final residue %s in section replacement (%s). Aborting.',res2f,msg.text));
    return
end
if length(indr2f) ~= 4
    add_msg_board(sprintf('ERROR: Template final residue addresa %s does not specify a residue in section replacement. Aborting.',res2f));
    return
end

[msg,coor1i] = get_residue(indr1i,'xyz_backbone'); 
if msg.error
    add_msg_board(sprintf('ERROR: Coordinates of target initial residue %s could not be retrieved in section replacement (%s). Aborting.',res1i,msg.text));
    return
end

[msg,coor1f] = get_residue(indr1f,'xyz_backbone'); 
if msg.error
    add_msg_board(sprintf('ERROR: Coordinates of target final residue %s could not be retrieved in section replacement (%s). Aborting.',res1f,msg.text));
    return
end

[msg,coor2i] = get_residue(indr2i,'xyz_backbone'); 
if msg.error
    add_msg_board(sprintf('ERROR: Coordinates of template initial residue %s could not be retrieved in section replacement (%s). Aborting.',res2i,msg.text));
    return
end

[msg,coor2f] = get_residue(indr2f,'xyz_backbone'); 
if msg.error
    add_msg_board(sprintf('ERROR: Coordinates of template final residue %s could not be retrieved in section replacement (%s). Aborting.',res2f,msg.text));
    return
end

[m1i,~] = size(coor1i);
[m2i,~] = size(coor2i);
if m1i ~= m2i
    add_msg_board('ERROR: Inconsistent initial residues in section replacement. Aborting.');
    return
end

[m1f,~] = size(coor1f);
[m2f,~] = size(coor2f);
if m1f ~= m2f
    add_msg_board('ERROR: Inconsistent final residues in section replacement. Aborting.');
    return
end
    
[~,~,transmati] = rmsd_superimpose(coor1i,coor2i); % affine transformation for initial residue
[transi,euleri] = affine2transrot(transmati); % additional translation vector and euler angles for final residue

[~,~,transmata] = rmsd_superimpose(coor1f,coor2f); % additional affine transformation required at final residue
[transa,eulera] = affine2transrot(transmata); % additional translation vector and euler angles for final residue

trans = zeros(2,3);
trans(1,:) = transi;
trans(2,:) = transa;

euler = zeros(2,3);
euler(1,:) = euleri;
euler(2,:) = eulera;

anchors = zeros(2,3);
anchors(1,:) = mean(coor2i);
anchors(2,:) = mean(coor2f);

if test
    return
end

numlinks0 = indr1f(4) - indr1i(4);

numlinks = indr2f(4) - indr2i(4);

if numlinks ~= numlinks0
    add_msg_board(sprintf('ERROR: Template (%i) and target (%i) section have different length in section replacement. Aborting.',numlinks,numlinks0));
    return
end

inconsistent = sum(abs(indr1i(1:3)-indr1f(1:3)));

if inconsistent
    add_msg_board(sprintf('ERROR: Target initial (%s) and final (%s) residue address do not belong to the same chain model in section replacement. Aborting.',res1i,res1f));
    return
end

inconsistent = sum(abs(indr2i(1:3)-indr2f(1:3)));

if inconsistent
    add_msg_board(sprintf('ERROR: Template initial (%s) and final (%s) residue address do not belong to the same chain model in section replacement. Aborting.',res2i,res2f));
    return
end


target_ind = indr1i;
template_ind = indr2i;

xyz = model.structures{target_ind(1)}(target_ind(2)).xyz{target_ind(3)};
xyz0 = model.structures{template_ind(1)}(template_ind(2)).xyz{template_ind(3)};

for k = 1:numlinks-1 % loop over all residues that have to be replaced
    target_ind(4) = indr1i(4)+k;
    template_ind(4) = indr2i(4)+k;
%     target_adr = mk_address(target_ind,true);
%     template_adr = mk_address(template_ind,true);
%     fprintf(1,'Replacing residue %s by residue %s\n',target_adr,template_adr);
    tags = model.structures{target_ind(1)}(target_ind(2)).residues{target_ind(3)}.info(target_ind(4)).atom_tags;
    atom_numbers = model.structures{target_ind(1)}(target_ind(2)).residues{target_ind(3)}.info(target_ind(4)).atom_numbers;
    nat = length(atom_numbers);
    temp_tags = model.structures{template_ind(1)}(template_ind(2)).residues{template_ind(3)}.info(template_ind(4)).atom_tags;
    temp_atom_numbers = model.structures{template_ind(1)}(template_ind(2)).residues{template_ind(3)}.info(template_ind(4)).atom_numbers;
    new_tags = ':';
    new_anums = cell(1,nat);
    atpoi = 0;
    for kat = 1:nat
        atag = id2tag(kat,tags);
        targ_anums = atom_numbers{kat};
        temp_kat = tag2id(atag,temp_tags);
        if ~isempty(temp_kat)
            new_tags = [new_tags atag ':'];
            atpoi = atpoi + 1;
            new_anums{atpoi} = targ_anums;
            temp_anums = temp_atom_numbers{temp_kat};
            [mtarg,~] = size(targ_anums);
            [mtemp,~] = size(temp_anums);
            if mtarg ~= mtemp
                add_msg_board(sprintf('ERROR: Mismatch of the number of locations for atom %s in template chain at residue %i. Aborting.',atag,template_ind(4)));
                return
            end
            for kl = 1:mtemp
                poi_temp = temp_anums(kl,1);
                xyz_temp = xyz0(poi_temp,:);
                xyz_targ = adjust_xyz(xyz_temp,anchors,trans,euler);
                poi_targ = targ_anums(kl,1);
                xyz(poi_targ,:) = xyz_targ;
            end       
        end
    end
    model.structures{target_ind(1)}(target_ind(2)).residues{target_ind(3)}.info(target_ind(4)).atom_tags = new_tags;
    model.structures{target_ind(1)}(target_ind(2)).residues{target_ind(3)}.info(target_ind(4)).atom_numbers = new_anums(1:atpoi);
end
model.structures{target_ind(1)}(target_ind(2)).xyz{target_ind(3)} = xyz;

function coor = adjust_xyz(coor0,anchors,trans,euler)

weights = zeros(1,length(anchors));
for ka = 1:2
    d = norm(coor0 - anchors(ka,:)) + 1e-6; % avoid accidental zero coordinate
    weights(ka) = 1/d;
end
weights = weights/sum(weights); % normalization
mean_trans = zeros(1,3);
mean_euler = zeros(1,3);
for ka = 1:2
    mean_trans = mean_trans + weights(ka)*trans(ka,:);
    mean_euler = mean_euler + weights(ka)*euler(ka,:);
end
transmat =  transrot2affine(mean_trans,mean_euler);
coor = affine_trafo_coor(coor0,transmat);

