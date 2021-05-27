load definitions/residues
for k = 1:20
    resname = id2tag(k,residue_defs.restags);
    residues.(upper(resname)).tlc = resname;
    residues.(upper(resname)).slc = residue_defs.single_letter_code(k);
    residues.(upper(resname)).id = k;
    residue_index(k).tlc = upper(resname);
end
    
for k = 1:8
    resname = strtrim(id2tag(k,residue_defs.nucleotide_tags));
    residues.(upper(resname)).tlc = resname;
    residues.(upper(resname)).slc = residue_defs.nucleotide_slc(k);
    residues.(upper(resname)).id = k + 20;
    residue_index(k+20).tlc = upper(resname);
end

for k = 1:8
    resname = strtrim(id2tag(k,residue_defs.cyana_nucleotide_tags));
    residues.(upper(resname)).tlc = resname;
    residues.(upper(resname)).slc = residue_defs.nucleotide_slc(k);
    residues.(upper(resname)).id = k + 20;
end

residue_slc = strcat(residue_defs.single_letter_code(1:20),residue_defs.nucleotide_slc);

save monomers residues residue_index residue_slc
