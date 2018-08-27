function site_scan_cofactors()
%SITE_SCAN_COFACTORS Sets global variables for a cofactor site scan 

global hMain

hMain.site_scan=1;
hMain.z_analysis=false;
hMain.statistics=false;
hMain.rotamer_PDB=false;
hMain.no_rot_pop=false;
hMain.dynamic_rotamers=false;
hMain.residue_pattern = '*';
hMain.site_scan_type = 'cofactor';
hMain.site_scan_inter = 0;
hMain.site_scan_intra = 0;
hMain.site_scan_homooligomer = 0;

end

