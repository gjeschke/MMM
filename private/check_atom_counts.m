function check_atom_counts(ename)

pdb_list = get_file_list(ename);
nhexpected = [1753,3624];

for k = 1:length(pdb_list)
    structure = rd_pdb(pdb_list{k});
    fprintf(1,'\n--- Conformer %s ---\n',pdb_list{k}); 
    for kc = 1:length(structure)
        num_heavy = sum(structure(kc).isotopes(:,1) ~= 1);
        num_prot = sum(structure(kc).isotopes(:,1) == 1);
        if num_heavy == nhexpected(kc)
            fprintf(1,'Chain %s has %i heavy atoms and %i protons\n',structure(kc).name,num_heavy,num_prot);
        else
            fprintf(2,'*** Chain %s has %i heavy atoms and %i protons ***\n',structure(kc).name,num_heavy,num_prot);
        end
    end
end