function [rmsd,transmat,msg] = backbone_overlap_peptide(adr1,adr2)
% [rmsd,transmat,msg] = backbone_overlap_peptide(adr1,adr2)
%
% Determine overlap by backbone rmsd of two peptide sections
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

madr1 = sprintf('%s.CA,N,C,O',adr1);
madr2 = sprintf('%s.CA,N,C,O',adr2);
indices1 = resolve_address(madr1);
indices2 = resolve_address(madr2);

[msg,coor1c] = get_object(indices1,'coor');
if msg.error
    return
end

m1 = length(coor1c);
coor1 = zeros(m1,3);
for k = 1:length(coor1c)
    coor1(k,:) = coor1c{k};
end

[msg,coor2c] = get_object(indices2,'coor');
if msg.error
    return
end

m2 = length(coor2c);
coor2 = zeros(m2,3);
for k = 1:length(coor2c)
    coor2(k,:) = coor2c{k};
end

if m1 ~= m2
    msg.error = 7;
    msg.text = 'Mismatch in the number of backbone atoms';
    return
end

[rmsd,~,transmat] = rmsd_superimpose(coor1,coor2);
