function [rmsd,transmat,msg] = overlap_nucleic_acid(adr1,adr2)
% [rmsd,transmat,msg] = overlap_nucleic_acid(adr1,adr2)
%
% Determine overlap by backbone rmsd of two nucleotide chain sections
%
% adr1      residue addresses of target model
% adr2      residue addresses of other model (must correspond to same 
%           number of residues_
%
% rmsd      backbone rmsd in Å
% transmat  affine transformation matrix that superimposes second onto
%           first section
% msg       error code (msg.error) and text (msg.text)

rmsd = [];
transmat = [];

indices1 = resolve_address(adr1);
indices2 = resolve_address(adr2);
[m1,n1] = size(indices1);
if m1 == 0 || n1 ~= 4
    msg.error = 1;
    msg.text = sprintf('Nucleotide indices of section %s could not be retrieved',adr1);
    return
end
[m2,n2] = size(indices2);
if m2 == 0 || n2 ~= 4
    msg.error = 2;
    msg.text = sprintf('Nucleotide indices of section %s could not be retrieved',adr2);
    return
end
if m1 ~= m2
    msg.error = 3;
    msg.text = sprintf('Different numbers of nucleotides in sections %s and %s',adr1,adr2);
    return
end

coor1 = zeros(100*m1,3);
coor2 = zeros(100*m1,3);
poi1 = 0;
poi2 = 0;

for k = 1:m1
    coor1a = get_nucleotide_coor(indices1(k,:));
    % [msg,coor1a] = get_residue(indices1(k,:),'xyz_base');
    if isempty(coor1a)
        msg.error = 4;
        msg.text = sprintf('Nucleotide base coordinates %s (%i) could not be retrieved',adr1,k);
        return
    end
    [mc,~] = size(coor1a);
    coor1(poi1+1:poi1+mc,:) = coor1a;
    poi1 = poi1 + mc;
    coor2a = get_nucleotide_coor(indices2(k,:));
    %[msg,coor2a] = get_residue(indices2(k,:),'xyz_base');
    if isempty(coor2a)
        msg.error = 4;
        msg.text = sprintf('Nucleotide base coordinates %s (%i) could not be retrieved',adr2,k);
        return
    end
    [mc,~] = size(coor2a);
    coor2(poi2+1:poi2+mc,:) = coor2a;
    poi2 = poi2 + mc;
end
coor1 = coor1(1:poi1,:);
coor2 = coor2(1:poi2,:);

[rmsd,~,transmat] = rmsd_superimpose(coor1,coor2);

msg.error = 0;
msg.text = 'No error.';

