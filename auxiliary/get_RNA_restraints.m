for rrm = 1:4,
    
    clear RNA_restraints

    switch rrm
        % translation table
        case 1 % RRM1/SL-E
            % RRM1
            trans(1).range = [57,154];
            trans(1).chain = '[RIB1](A)';
            trans(1).offset = 0;
            % SL-E, CLIR-MS
            trans(2).range = [228,232];
            trans(2).chain = '[1](A)';
            trans(2).offset = -224;
            upl = read_Cyana_restraints('unambigous_RRM1.upl',trans);
            H_lol = read_Cyana_restraints('RRM1_SLE_hbonds.lol',trans);
            H_upl = read_Cyana_restraints('RRM1_SLE_hbonds.upl',trans);
        case 2
            % RRM2
            trans(1).range = [182,283];
            trans(1).chain = '[RIB2](A)';
            trans(1).offset = 0;
            % SL-F, NMR
            trans(2).range = [348,370];
            trans(2).chain = '[1](A)';
            trans(2).offset = -347;
            upl = read_Cyana_restraints('RRM2_NMR_for_Gunnar.upl',trans);
            H_lol = read_Cyana_restraints('RRM2_inter_hbonds_for_Gunnar.lol',trans);
            H_upl = read_Cyana_restraints('RRM2_inter_hbonds_for_Gunnar.upl',trans);
        case 3
            % RRM34
            trans(1).range = [337,531];
            trans(1).chain = '[RIB3](A)';
            trans(1).offset = 0;
            % SL-D
            trans(2).range = [594,598];
            trans(2).chain = '[1](A)';
            trans(2).offset = -580;
            upl = read_Cyana_restraints('unambigous_RRM3.upl',trans);
            H_lol = read_Cyana_restraints('RRM3_hbonds_for_Gunnar.lol',trans);
            H_upl = read_Cyana_restraints('RRM3_hbonds_for_Gunnar.upl',trans);
        case 4
            % RRM34
            trans(1).range = [337,531];
            trans(1).chain = '[RIB3](A)';
            trans(1).offset = 0;
            % Link E-F
            trans(2).range = [592,601];
            trans(2).chain = '[1](A)';
            trans(2).offset = -591;
            upl = read_Cyana_restraints('unambigous_RRM4.upl',trans);
            H_lol = read_Cyana_restraints('RRM4_hbonds_for_Gunnar.lol',trans);
            H_upl = read_Cyana_restraints('RRM4_hbonds_for_Gunnar.upl',trans);
    end
    
    poi = 0;
    for k = 1:length(upl)
        poi = poi + 1;
        RNA_restraints(poi).r = upl(k).r;
        RNA_restraints(poi).type1 = upl(k).type1;
        RNA_restraints(poi).type2 = upl(k).type2;
        RNA_restraints(poi).adr1 = upl(k).adr1;
        RNA_restraints(poi).adr2 = upl(k).adr2;
    end
    for k = 1:length(H_lol)
        poi = poi + 1;
        RNA_restraints(poi).r = -H_lol(k).r;
        RNA_restraints(poi).type1 = H_lol(k).type1;
        RNA_restraints(poi).type2 = H_lol(k).type2;
        RNA_restraints(poi).adr1 = H_lol(k).adr1;
        RNA_restraints(poi).adr2 = H_lol(k).adr2;
    end
    for k = 1:length(H_upl)
        poi = poi + 1;
        RNA_restraints(poi).r = H_upl(k).r;
        RNA_restraints(poi).type1 = H_upl(k).type1;
        RNA_restraints(poi).type2 = H_upl(k).type2;
        RNA_restraints(poi).adr1 = H_upl(k).adr1;
        RNA_restraints(poi).adr2 = H_upl(k).adr2;
    end
    
    save(sprintf('RRM%i_RNA_restraints.mat',rrm),'RNA_restraints');
end

