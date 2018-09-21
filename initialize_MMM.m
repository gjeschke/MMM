function initialize_MMM
% function initialize_MMM
%
% Loads definition files and declares global variables
%
% G. Jeschke, 2009-2014

global MMM_info
global graph_settings
global geometry_settings
global chemistry
global residue_defs
global label_defs
global model
global glass_transitions
global rotamer_libraries
global ligand_libraries
global web_adr
global help_files
global rootdir
global MMM_icon
global PDBwiki_icon
global proteopedia_icon
global general
global queries
global formats
global third_party
global annotation_keys
global membrane_profiles
global eav
global DPPC
global DOPC
global ENM_param
global browser
global Modeller_info
global LJ
global uff

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
warning('off','MATLAB:hg:uicontrol:SingleLineEditCannotHaveMultiLineText');
model=[];

browser='system';

gowhere=which('MMM.m');

% addpath('C:\FCCC\Scwrl4\'); % adds SCWRL4 path

addpath(genpath(gowhere(1:end-6)));

% information on current release
MMM_info.title='MMM version 2018.1';
MMM_info.authors='G. Jeschke, Ye. Polyhach, S. Stoll';
MMM_info.date='20.07.2018';
MMM_info.institution='ETH Zürich, D-CHAB & Univ. Washington, Dept. Chem.';
MMM_info.contact{1}='gjeschke@ethz.ch';
MMM_info.contact{2}='yevhen.polyhach@phys.chem.ethz.ch';
MMM_info.contact{3}='stst@uw.edu';

% general settings

general.tmp_check=30; % time interval (in days) at which MMM asks about deleting old temporary files
general.old_file=30; % time interval (in days) after which a temporary file is considered to be old
general.cpu=6; % number of CPU cores to be used in parallel computing, is overwritten by user preferences
general.timers={};
general.large_windows = false; % allow for windows larger than 800*600 screen resolution

%output formats
formats.rtf_preamble='{\rtf1\ansi\ansicpg1252\deff0\deflang2055{\fonttbl{\f0\froman\fcharset0 Times New Roman;}{\f1\fswiss\fcharset0 Arial;}}';
formats.rtf_origin='{\*\generator MMM 2015.1;}\viewkind4\uc1\pard\f0\fs24';
formats.rtf_colors='{\colortbl ;\red255\green0\blue0;\red0\green255\blue0;\red0\green0\blue255;}'; % a basic primary color table for RTF files
% Web addresses of data bases, prediction servers etc.

web_adr.SFX='sfx.ethz.ch:9003/sfx_locater';
% web_adr.SFX='http://sfx.kobv.de/sfx_tub';
web_adr.PDB='ftp.ebi.ac.uk';
web_adr.PDB_EU='ftp.ebi.ac.uk';
web_adr.PDB_US='ftp.wwpdb.org';
web_adr.PDB_J='ftp.protein.osaka-u.ac.jp';
web_adr.PDB_explore='http://www.pdb.org/pdb/explore/';
web_adr.PDB_queries='http://www.pdb.org/pdb/files/';
web_adr.proteopedia='http://www.proteopedia.org/wiki/index.php/';
web_adr.PDBwiki='http://www.pdbwiki.org/index.php/';
web_adr.OCA='http://oca.weizmann.ac.il';
web_adr.EDS='http://eds.bmc.uu.se';
web_adr.PubMed='https://www.ncbi.nlm.nih.gov/pubmed/';
web_adr.PubMed_interactive='https://www.ncbi.nlm.nih.gov/pubmed/';
web_adr.UniProt='http://www.uniprot.org/uniprot/';
web_adr.DOI_resolver='http://dx.doi.org/';
web_adr.StoneHinge='http://stonehinge.molmovdb.org/submit.html';
web_adr.HingeMaster='http://molmovdb.org/cgi-bin/submit-flexoracle.cgi';
web_adr.Spritz='http://distill.ucd.ie/spritz/';
web_adr.MEMSAT3='http://bioinf.cs.ucl.ac.uk/psipred/';
web_adr.TMHMM='http://www.cbs.dtu.dk/services/TMHMM/';
web_adr.Philius='http://www.yeastrc.org/philius/pages/philius/runPhilius.jsp';
web_adr.Hex='http://www.csd.abdn.ac.uk/hex_server/';
web_adr.RosettaDock='http://rosettadock.graylab.jhu.edu/';
web_adr.Haddock='http://haddock.chem.uu.nl/Haddock/haddock.php';
web_adr.HHpred='http://toolkit.tuebingen.mpg.de/hhpred';
web_adr.Phyre='http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index';
web_adr.SAM_sec='http://compbio.soe.ucsc.edu/SAM_T08/T08-query.html';
web_adr.SAM_tert='http://compbio.soe.ucsc.edu/SAM_T08/T08-query.html';
web_adr.PSIpred='http://bioinf.cs.ucl.ac.uk/psipred/';
web_adr.Porter='http://distill.ucd.ie/porter/';
web_adr.PredictProtein='http://www.predictprotein.org/';
web_adr.Sable='http://sable.cchmc.org/';
web_adr.DichroCalc='http://comp.chem.nottingham.ac.uk/cgi-bin/dichrocalc/bin/getparams.cgi';
web_adr.CryCo='http://ligin.weizmann.ac.il/~lpgerzon/cryco5.0/cryco';
web_adr.PRONOX='http://rockscluster.hsc.usc.edu/research/software/pronox/pronox.html';
web_adr.NCBI_domains='http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi';
web_adr.ExPASy='http://www.expasy.org/tools/';
%web_adr.OCA='http://oca.weizmann.ac.il/oca-bin';

% Web queries
queries.PMID='?report=medline&format=text';
queries.search_PubMed='&report=medline&format=text';
queries.tool='&tool=MMM';
queries.email='&email=gunnar.jeschke@phys.chem.ethz.ch';
queries.PubMed_file='PubMed.txt';
queries.update=7; % minimum time between reference updates in days
queries.PDB_monomers='pub/pdb/data/monomers/'; % monomer directoiry of PDB
queries.PDB_structures='/pub/databases/rcsb/pdb-remediated/data/structures/all/pdb/'; % structure directory of PDB Europe
queries.PDB_structures_EU='/pub/databases/rcsb/pdb-remediated/data/structures/all/pdb/'; % structure directory of PDB Europe
queries.PDB_structures_US='/pub/pdb/data/structures/all/pdb/';
queries.PDB_structures_J='/pub/pdb/data/structures/all/pdb/';
queries.PDB_explore_id='explore.do?structureId=';
queries.PDB_fasta='fasta.txt?structureIdList='; % use str=urlread([web_adr.PDB_queries queries.PDB_fasta pdbid]); to retrieve FASTA
queries.OCA='/oca-bin/ocashort?id=';
queries.EDS='/cgi-bin/eds/uusfs?pdbCode=';
queries.SFX_preamble='?&url_ver=Z39.88-2004&url_ctx_fmt=info:ofi/fmt:kev:mtx:ctx&rft_val_fmt=info:ofi/fmt:kev:mtx:journal';
queries.SFX_title='&rft.atitle=';
queries.SFX_initial='&rft.auinit=';
queries.SFX_surname='&rft.aulast=';
queries.SFX_journal='&rft.btitle=';
queries.SFX_year='&rft.date=';
queries.SFX_first_page='&rft.spage=';
queries.SFX_last_page='&rft.epage=';
queries.SFX_issue='&rft.issue';
queries.SFX_volume='&rft.volume=';
queries.SFX_author='&rft.au=';
queries.SFX_sep='%2C%20';
queries.SFX_space='%20';
queries.SFX_DOI='&rft_id=info:doi';
queries.SFX_genre='&rft.genre=article';
queries.UniProt_fasta='.fasta';
queries.CryCo='/cryst2.cgi?id=';
queries.CSU='/oca-bin/csu?PDB_ID=';

load definitions/PSE
chemistry.pse=pse;
chemistry.element_tags=element_tags;

load definitions/colors

load definitions/LJ_charmm27

graph_settings.scroll_speed=0.02; % speed of zooming with the scroll wheel
graph_settings.max_full_objects=10000;    % maximum atom graphics objects
                                          % before MMM suppresses atoms
                                          % in rotate, zoom, and pan
graph_settings.max_reduced_objects=25000; % maximum reduced atom graphics 
                                          % before MMM suppresses reduced
                                          % atom display in rotate, zoom,
                                          % and pan

graph_settings.az=49; % standard 3D view on loading protein
graph_settings.el=26;
graph_settings.colors=colors;
graph_settings.color_tags=color_tags;
graph_settings.color_names=names;
graph_settings.selected_color=[1 1 0]; % color for selected elements
graph_settings.maxbonds=12;
graph_settings.wire_width=1; % line width for bond wires
graph_settings.CaWire_width=2; % line width for C_alpha trace wires 
graph_settings.stick_radius=0.2;
graph_settings.CaStick_radius=0.35;
graph_settings.ball_radius=0.4;
graph_settings.sprh=3; %segments per half residue in coil plots
graph_settings.spr=5; %segments per residue in ribbon plots
graph_settings.NA_spr=3; %segments per residue in DNA/RNA ribbon plots
graph_settings.loop_color=[250 250 210]/255; %[0,1,0];
graph_settings.loop_radius=0.4;
graph_settings.coil_radius=0.35;
graph_settings.helix_color=[255 160 122]/255; % lightsalmon
graph_settings.helix_width=2.5; % 3.0
graph_settings.helix_height=0.4; % 0.4
graph_settings.TM_helix_color=[255  69   0]/255; % orangered
graph_settings.helix_cartoon_radius = 1.5;
graph_settings.helix_cartoon_points = 15;
% graph_settings.sheet_color=[176 196 222]/255; %[0,0,1];
graph_settings.sheet_color=[70 130 180]/255; % steelblue
graph_settings.sheet_width=2.0; % 2.0
graph_settings.sheet_height=0.6; % 0.6
graph_settings.sheet_arrow=1.75;
graph_settings.TM_sheet_color=[0 128 128]/255; % teal
graph_settings.NA_color=[255, 228, 181]/255;
graph_settings.NA_width=1.5; % 1.5
graph_settings.NA_height=0.6; % 0.6
graph_settings.NA_tube_radius=0.35;
graph_settings.CaWire_color=[128,0,0]/255; % maroon
graph_settings.CaStick_color=[128,64,0]/255; % slightly lighter than maroon
graph_settings.coil_color=[50,205,50]/255; % lime green
graph_settings.DSSP_color_H=graph_settings.helix_color;
graph_settings.DSSP_color_B=[135 206 250]/255;
graph_settings.DSSP_color_E=graph_settings.sheet_color;
graph_settings.DSSP_color_G=[244 164 96]/255;
graph_settings.DSSP_color_I=[255 69 0]/255;
graph_settings.DSSP_color_T=[0 100 0]/255;
graph_settings.DSSP_color_S=[34 139 34]/255;
graph_settings.DSSP_color_C=graph_settings.coil_color;
graph_settings.bilayer_midplane_color=[152 251 152]/255; % palegreen
graph_settings.bilayer_outer_color=[255 105 180]/255; % hotpink
graph_settings.decent_string=[0.94,1.00,0.94]; % honeydew
graph_settings.motion_arrow=[0.86,0.08,0.24]; % crimson

graph_settings.min_hydropathy=-4.5;
graph_settings.max_hydropathy=4.5;
graph_settings.min_helix_propensity=0;
graph_settings.max_helix_propensity=3.16;

graph_settings.Bmin=round(sqrt(8*pi^2*0.5^2)); % sqrt(B factor) corresponding to "minimum color" on B factor scale
graph_settings.Bmax=round(sqrt(8*pi^2*1.5^2)); % sqrt(B factor) corresponding to "maximum color" on B factor scale

graph_settings.label_radius=(1.55+1.262+1.52)/2; % label sphere radius for 100% of population, 
                                                 % from van-der-Waals radii of N, O and N-O bond length
                                                 
graph_settings.SAS_transparency=0.3; % default transparency for solvent accessible surface
graph_settings.SAS_color=[127,255,212]/255; % default color for solvent accessible surface (aquamarine)

geometry_settings.atom_context=3.5;
geometry_settings.residue_context=8;
geometry_settings.chain_context=40;

glass_transitions.protein=175; % see W. Doster, Eur Biophys J (2008) 37:591–602

% Membrane profiles 
membrane_profiles.O2.mindel=9.0; % [usec^{-1}]
membrane_profiles.O2.maxdel=16.9; % [usec^{-1}]
membrane_profiles.O2.d0=9.8; % [Å]
membrane_profiles.O2.lambda=1.4; % [Å]

% Computation of effective accessible volume of spin labels
eav.NiEDDA.name='NiEDDA';
eav.NiEDDA.decay_length=5.5;
eav.NiEDDA.blocking_radius=5.5;
eav.NiEDDA.kex_factor=0.632;
eav.NiEDDA.Pi_factor=3*eav.NiEDDA.kex_factor; % Hubbell, definition for 3 mM NiEDDA
eav.NiEDDA.sig=5.55;
eav.NiEDDA.bulk=0;
eav.NiEDDA.amp=10;
eav.NiEDDA.buffer=10;
eav.NiEDDA.offset=-30; % offset of layer with maximum relaxivity from bilayer centre in Å, must be negative when non-zero
eav.O2.name='O2';
eav.O2.decay_length=1.0;
eav.O2.blocking_radius=2.5;
eav.O2.kex_factor=9.835271; % 2.72374; % 1.640848; % 1.27267;
eav.O2.Pi_factor=0.27*eav.O2.kex_factor; % Hubbell, definition for 0.27 mM O2, equil. 298 K, 1 atm, air/water
eav.O2.sig=4.5;
eav.O2.bulk=2.0;
eav.O2.buffer=2;
eav.O2.amp=47;
eav.O2.offset=0; % offset of layer with maximum relaxivity from bilayer centre in Å
eav.water.name='D2O';
eav.water.decay_length=3;
eav.water.blocking_radius=2.5;
eav.water.kex_factor=3.395;
eav.water.Pi_factor=1; 
eav.water.sig=5.55;
eav.water.bulk=0.1;
eav.water.amp=57.36;
eav.water.buffer=57.46;
eav.water.weight_grid_name='ESEEM_weight_grid';

% Definition of existing rotamer libraries, all libraries of a given label
% type are listed in one structure, further information is in the library
% files themselves

rotamer_libraries(1).label='MTSL';
rotamer_libraries(1).tc='R1A';
rotamer_libraries(1).T=[175,298];
rotamer_libraries(1).files=':R1A_175K_UFF_216_CASD:R1A_298K_UFF_216_r1_CASD:';
rotamer_libraries(1).MC='R1A_298K_UFF_216_r#_CASD';
% rotamer_libraries(1).files=':R1A_175K_090619:R1A_298K_090619:'; 
rotamer_libraries(1).type='peptide';
rotamer_libraries(1).exclude=':H2:OXT:HXT:H:'; % for peptides, these atoms are defined in the rotamer, but not attached

rotamer_libraries(2).label='IA-PROXYL';
rotamer_libraries(2).tc='IA1';
rotamer_libraries(2).T=[175,298];
rotamer_libraries(2).files=':IA1_175K_090814:IA1_298K_090814:';
rotamer_libraries(2).type='peptide';
rotamer_libraries(2).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(3).label = 'IAP-4TU';
rotamer_libraries(3).tc = 'TUP';
rotamer_libraries(3).T = 298;
rotamer_libraries(3).files = ':TUP_298K_UFF_72:';
rotamer_libraries(3).type = 'nucleotide';
rotamer_libraries(3).exclude = ':N1:C2:O2:N3:H3:C4:O4:C5:H5:C6:H6:N4:H42:H41:'; % for uracil anc cytosine, these atoms come from the rotamer or are not present

rotamer_libraries(4).label='MA-PROXYL';
rotamer_libraries(4).tc='MA1';
rotamer_libraries(4).T=[175,298];
rotamer_libraries(4).files=':MalP_175K_UFF_108:MA1_298K_UFF_sticky_108:';
rotamer_libraries(4).type='peptide';
rotamer_libraries(4).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(5).label='Br-MTSL';
rotamer_libraries(5).tc='R7A';
% rotamer_libraries(5).T=[175,298];
% rotamer_libraries(5).files=':R7A_175K_UFF_216_CASD:R7A_298K_UFF_216_CASD:';
rotamer_libraries(5).T= 298;
rotamer_libraries(5).files=':R7A_298K_UFF_216_r1_CASD:';
rotamer_libraries(5).type='peptide';
rotamer_libraries(5).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(6).label='M8DOTA-Gd';
rotamer_libraries(6).tc='M8D';
rotamer_libraries(6).T=298;
rotamer_libraries(6).files=':M8D_298K_UFF_1944_r2:';
rotamer_libraries(6).type='peptide';
rotamer_libraries(6).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(7).label = 'MTS-4TU';
rotamer_libraries(7).tc = 'TUM';
rotamer_libraries(7).T = 298;
rotamer_libraries(7).files = ':TUM_298K_UFF_192:';
rotamer_libraries(7).type = 'nucleotide';
rotamer_libraries(7).exclude = ':N1:C2:O2:N3:H3:C4:O4:C5:H5:C6:H6:N4:H42:H41:'; % for uracil and cytosine, these atoms come from the rotamer or are not present

rotamer_libraries(8).label = 'R5-TP';
rotamer_libraries(8).tc = 'R5T';
rotamer_libraries(8).T = 298;
rotamer_libraries(8).files = ':R5T_298K_UFF_360:';
rotamer_libraries(8).type = 'nucleotide';
rotamer_libraries(8).exclude = ':OP1:'; % this atom comes from the rotamer (as an S atom)

rotamer_libraries(9).label = 'CNC-NO';
rotamer_libraries(9).tc = 'CNR';
rotamer_libraries(9).T = 298;
rotamer_libraries(9).files = ':CNR_298K_UFF_144:';
rotamer_libraries(9).type = 'cofactor';
rotamer_libraries(9).exclude = ':O8R:'; % this atom comes from the rotamer

rotamer_libraries(10).label='DOTA-Gd';
rotamer_libraries(10).tc='GMO';
rotamer_libraries(10).T=298;
rotamer_libraries(10).files=':GMO_298K_UFF_648:';
rotamer_libraries(10).type='peptide';
rotamer_libraries(10).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(11).label='DTPA-Gd';
rotamer_libraries(11).tc='GTO';
rotamer_libraries(11).T=298;
rotamer_libraries(11).files=':GTO_298K_UFF_2430:'; %':GTO_298K_UFF_noH_T6_648:';
rotamer_libraries(11).type='peptide';
rotamer_libraries(11).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(12).label='HF-K1';
rotamer_libraries(12).tc='K1H';
rotamer_libraries(12).T=298;
rotamer_libraries(12).files=':RUA_298K_UFF_288:';
rotamer_libraries(12).type='peptide';
rotamer_libraries(12).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(13).label='R5-TPT';
rotamer_libraries(13).tc='RTT';
rotamer_libraries(13).T=298;
rotamer_libraries(13).files=':RTT_298K_UFF_576:';
rotamer_libraries(13).type='nucleotide';
rotamer_libraries(13).exclude=':O5'':OP1:OP2:P:H05'':H5'':H5''''';

rotamer_libraries(14).label='Tormyshev-Trityl';
rotamer_libraries(14).tc='TMT';
rotamer_libraries(14).T=298;
rotamer_libraries(14).files=':TMT_298K_UFF_3888_r2:';
rotamer_libraries(14).type='peptide';
rotamer_libraries(14).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(15).label='V1';
rotamer_libraries(15).tc='V1A';
rotamer_libraries(15).T=298;
rotamer_libraries(15).files=':V1A_298K_UFF_72_r3:';
rotamer_libraries(15).type='peptide';
rotamer_libraries(15).exclude=':H2:OXT:HXT:H:'; % for peptides, these atoms are defined in the rotamer, but not attached

rotamer_libraries(16).label='cNox@Tyr';
rotamer_libraries(16).tc='Nc1';
rotamer_libraries(16).T=298;
rotamer_libraries(16).files=':Nc1_298K_UFF_128_r1:';
rotamer_libraries(16).type='peptide';
rotamer_libraries(16).exclude=':H2:OXT:HO2:H:'; % for peptides, these atoms are defined in the rotamer, but not attached

rotamer_libraries(17).label='lNox@Tyr';
rotamer_libraries(17).tc='Nx1';
rotamer_libraries(17).T=298;
rotamer_libraries(17).files=':Nx1_298K_UFF_256_r1:';
rotamer_libraries(17).type='peptide';
rotamer_libraries(17).exclude=':H2:OXT:HO2:H:'; % for peptides, these atoms are defined in the rotamer, but not attached

rotamer_libraries(18).label='GPymiMTA';
rotamer_libraries(18).tc='GPM';
rotamer_libraries(18).T=298;
rotamer_libraries(18).files=':GPM_298K_UFF_432.mat:';
rotamer_libraries(18).type='peptide';
rotamer_libraries(18).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(19).label='dHis-Cu';
rotamer_libraries(19).tc='HCU';
rotamer_libraries(19).T= 298;
rotamer_libraries(19).files=':HCU_298K_UFF_12_r5:';
rotamer_libraries(19).type='peptide';
rotamer_libraries(19).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(20).label='DZD';
rotamer_libraries(20).tc='DZD';
rotamer_libraries(20).T= 298;
rotamer_libraries(20).files=':DZD_298K_UFF_216:';
rotamer_libraries(20).type='peptide';
rotamer_libraries(20).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(21).label='IAG';
rotamer_libraries(21).tc='GDI';
rotamer_libraries(21).T= 298;
rotamer_libraries(21).files=':GDI_298K_UFF_f09_2461:';
rotamer_libraries(21).type='peptide';
rotamer_libraries(21).exclude=':H2:OXT:HXT:H:';

rotamer_libraries(22).label='MAG';
rotamer_libraries(22).tc='GDM';
rotamer_libraries(22).T= 298;
rotamer_libraries(22).files=':GDM_298K_UFF_f09_1367:';
rotamer_libraries(22).type='peptide';
rotamer_libraries(22).exclude=':H2:OXT:HXT:H:';

ligand_libraries(1).label='dHisCuIDA';
ligand_libraries(1).tc='IDA';
ligand_libraries(1).files=':dHis_Cu_IDA_asym:';
ligand_libraries(1).dist = [2.6,3.2];
ligand_libraries(1).V_ksi = 10000; % J/mol
ligand_libraries(1).res=':HCU:HCU:';
ligand_libraries(1).type='peptide';

ligand_libraries(2).label='dHisCuNTA';
ligand_libraries(2).tc='NTA';
ligand_libraries(2).files=':dHis_Cu_NTA_sym:';
ligand_libraries(2).V_ksi = 10000; % J/mol
ligand_libraries(2).dist = [2.6,3.2];
ligand_libraries(2).res=':HCU:HCU:';
ligand_libraries(2).type='peptide';

ligand_libraries(3).label='dHisCuIDAsym';
ligand_libraries(3).tc='IDS';
ligand_libraries(3).files=':dHis_Cu_IDA_sym:';
ligand_libraries(3).dist = [2.6,3.2];
ligand_libraries(3).V_ksi = 10000; % J/mol
ligand_libraries(3).res=':HCU:HCU:';
ligand_libraries(3).type='peptide';

% ligand_libraries(3).label='fake';
% ligand_libraries(3).tc='NTA';
% ligand_libraries(3).files=':dHis_Cu_NTA_sym:dHis_Cu_NTA_asym:';
% ligand_libraries(3).dist = [2.6,3.0,3.3];
% ligand_libraries(3).res=':R1A:R1A:';
% ligand_libraries(3).type='peptide';

Modeller_info.label_codes='!;';
Modeller_info.het_codes='i$z3obmh';

load definitions/residues
load definitions/labels

rootdir=which('MMM.m');
rootdir=rootdir(1:end-5);
help_files=['file:///' rootdir 'help/'];

general.tmp_files=[rootdir 'tmp' filesep];
general.rootdir=rootdir;
general.help_files=[rootdir 'help' filesep];
general.help_home='overview.html';
general.help_changes='changes.html';
general.scripts=strcat(rootdir,'scripts/');
general.Ramachandran=strcat(rootdir,'Ramachandran/');
general.virgin_path=true;
general.exp_files=pwd;
general.pdb_files=pwd;
general.DEER_files=pwd;
general.reports=pwd;
general.restraint_files=pwd;
general.biblio_files=pwd;
general.reports=pwd;

load third_party_references

third_party.tags=auto_ref_tags;
third_party.references=references;
third_party.modeller_version='mod9.9';

load bilayer_definitions

load uff_parameters % universal force field as implemented in Towhee

% elastic network model (ENM) parameters
% note that part of these parameters are overwritten by user preferences
% that are processeed by set_ANM

ENM_param.gamma=2.06; % generic force constant divided by k_B T in Å^{-2}
ENM_param.rc=7.3; % cutoff distance in Å, should be 7.3
ENM_param.rc_ANM=sqrt(3)*ENM_param.rc; % cutoff distance for anisotropic network model
ENM_param.p_GNM=2; % exponent for distance-dependent force constants GNM
ENM_param.p_GNM_gamma=2.06; % exponent for distance-dependent force constants GNM
ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
ENM_param.p_ANM_gamma=2.06; %  exponent for distance-dependent force constants ANM
ENM_param.imANM_scaling=16; % see T. R. Lezon, I. Bahar, Biophys. J. 2012, 102, 1331-1340

ENM_param.fit_basis=20; % number of slow modes used in structure fitting
ENM_param.fix_local=2; % number of neighbors to which distances are fixed during fitting,
                       % 1 stabilizes only C_alpha-C_alpha distances
                       % 2 stabilizes C_alpha-C_alpha distances as well as
                       % C_alpha-C_alpha-C_alpha angles
		               % -2 stabilizes only C_alpha-C_alpha-C_alpha angles,
		               % assuming that the ENM stabilizes C_alpha-C_alpha
		               % distances
                       % 3 stabilizes local secondary structure (not yet
                       % implemented)
ENM_param.tol_local=0.05; % tolerance (r.m.s.d. in Å) for local distortions
ENM_param.fix_force=10; % see Zheng/Brooks (2006), below Eq. (4)
ENM_param.sec_force=0.02; % see Zheng/Brooks (2006), below Eq. (4), different 
                        % restoring force for second neighbor
ENM_param.tol_dist=3.0; % tolerance of distance restraints in structure fitting
ENM_param.tol_model=0.5; % expected uncertainty of the model
ENM_param.track=false;   % initial mode set is tracked during fitting, if true, 
                        % otherwise the lowest-energy modes of new diagonalization are used
ENM_param.cycles=100; % maximum number of iterations in ENM-based fitting
ENM_param.diagonalize=false; % if true, always recompute normal modes by Hessian rediagonalization
ENM_param.reorientate=false; % if true, always use reorientation of normal modes rather than recomputation
ENM_param.update_frequency=5; % number of cycles, after which normal modes are updated by Hessian diagonalization
ENM_param.tif=0.075; % threshold for active space extension, can be overwritten in constraint file
ENM_param.mmax=30; % iteration at which active space achieves maximum dimension, can be overwritten in constraint file
ENM_param.relax=false;

ENM_param.imANM=false; % is overwritten by user preference
ENM_param.parametrization='ed-ENM'; % is overwritten by user preference
ENM_param.mass_weighting=false; % is overwritten by user preference

pfname=which('preferences.mat');
load preferences
if isfield(user_preferences,'ANM'), 
    ENM_param=set_ANM(ENM_param,user_preferences.ANM);
    ENM_param=set_ANM(ENM_param,user_preferences.ANM);
end;
if isfield(user_preferences,'cpu'),
    general.cpu=user_preferences.cpu;
end;
if isfield(user_preferences,'large'),
    general.large_windows= user_preferences.large;
end;

if isfield(user_preferences,'Modeller_call'),
    third_party.modeller_version=user_preferences.Modeller_call;
end;

% import preferences
if exist([general.rootdir 'preferences.mat'],'file'),
    load([general.rootdir 'preferences.mat']);
    if exist('user_preferences','var'),
        if isfield(user_preferences,'SFX'),
            web_adr.SFX=user_preferences.SFX;
        end;
        if isfield(user_preferences,'PDB'),
            web_adr.PDB=user_preferences.PDB;
        end;
        if isfield(user_preferences,'SFX'),
            queries.PDB_structures=user_preferences.PDB_structures;
        end;
        if isfield(user_preferences,'browser'),
            browser=user_preferences.browser;
        end;
    end;
end;


if virgin,
    home=[general.help_files general.help_changes];
    webcall(home,'-helpbrowser');
    virgin=false;
    save(pfname,'virgin','user_preferences');
end;

annotation_keys(1).type='Alternate';
annotation_keys(1).keywords={'alternate locations'};
annotation_keys(2).type='Binding';
annotation_keys(2).keywords={'binding sites'};
annotation_keys(3).type='Coordinating';
annotation_keys(3).keywords={'metal centers'};
annotation_keys(4).type='General';
annotation_keys(4).keywords={};
annotation_keys(5).type='Inserted';
annotation_keys(5).keywords={'inserted residues'};
annotation_keys(6).type='Metal';
annotation_keys(6).keywords={'metal centers'};
annotation_keys(7).type='Missing';
annotation_keys(7).keywords={'missing atoms'};
annotation_keys(8).type='Mutation';
annotation_keys(8).keywords={'mutations'};
annotation_keys(9).type='Spin';
annotation_keys(9).keywords={'spin label attached' 'rotamers computed'};

MMM_icon(:,:,1) = [ ...
	191 191 191 191 191 191 254 225 174 144 135 131 128 125 122 120 117 114 111 109 106 103 102 107 105 102 100  97  94  92  89  86 ; ...
	191 191 191 191 191 239 162 137 137 137 134 131 128 126 123 120 117 115 112 109 106  97 167 240 240 239 239 239 239 238 240 142 ; ...
	191 191 191 191 242 168 143 143 140 137 134 131 129 126 123 120 117 115 112 109 105 107 234 255 255 255 255 255 255 255 211  77 ; ...
	191 191 191 244 205 164 144 143 140 137 134 131 129 126 123 120 117 115 112 109 104 117 182 182 180 174 226 255 255 242 103  71 ; ...
	191 191 243 171 166 149 144 142 139 137 134 131 128 125 123 120 117 114 112 109 106 102  94  91  82 127 243 255 252 146  71  75 ; ...
	191 191 166 143 145 144 143 140 138 136 133 130 128 125 122 119 117 114 111 108 106 103 100  94 120 239 255 255 236 101  72  75 ; ...
	191 231 137 143 143 142 140 139 136 134 132 129 127 124 121 119 116 113 111 108 105 102  99  95 184 254 255 255 255 228  96  71 ; ...
	191 179 137 140 140 139 138 136 135 132 130 128 125 123 120 118 115 112 110 107 104 102  99  96  94 127 229 255 255 255 190  70 ; ...
	238 146 137 137 137 137 136 134 132 130 128 126 124 122 119 117 114 111 109 106 104  97  92  89  91  82 123 250 255 255 236  83 ; ...
	189 135 134 135 135 134 133 132 130 128 127 124 122 120 118 115 113 110 108 105  99 139 183 165 101  86  90 230 255 255 234  81 ; ...
	136 131 131 131 130 130 129 128 127 125 123 121 120 118 116 114 111 109 106  99 147 252 255 255 209  86  99 242 255 255 192  68 ; ...
	127 128 128 131 134 134 133 132 131 130 128 126 117 116 115 112 110 108 105  96 195 255 255 255 248 107 179 255 255 236 101  68 ; ...
	125 125 129 220 239 239 239 239 238 238 238 237 201 116 112 111 108 106 102  95 152 255 255 255 255 235 255 253 206 105  68  71 ; ...
	123 118 177 255 255 255 255 255 255 255 255 255 247 123 109 109 105 100 109 106  92 148 200 205 204 192 160 115  77  71  73  69 ; ...
	116 132 240 255 255 255 255 255 255 255 255 255 218 109 109 103 117 197 235 233 142  84  86  85  82  78  74  75  76  74  71  68 ; ...
	115 126 148 147 147 147 144 163 252 255 255 255 174 103 104 111 220 255 255 255 202  86  89  87  85  79  76  76  75  73  70  67 ; ...
	114 114 111 111 111 111 105 167 255 255 255 251 132 101 105 206 255 255 255 255 171  84  88  85  80 104 116  85  71  72  69  66 ; ...
	111 112 112 112 112 112 108 212 255 255 255 224 105 100 192 255 255 255 255 252 120  85  85  81 173 247 255 218  87  68  67  65 ; ...
	109 109 109 109 109 107 123 245 255 255 255 180  91 177 255 255 255 255 255 224  90  84  82 179 255 255 255 255 129  62  66  63 ; ...
	106 106 107 107 107 101 156 255 255 255 253 128 158 254 253 246 255 255 255 177  77  88 190 255 255 255 255 250 101  63  64  62 ; ...
	103 103 104 104 104  98 201 255 255 255 224 150 250 255 188 233 255 255 255 123  89 205 255 255 255 255 255 208  66  65  63  60 ; ...
	100 101 101 101  99 110 239 255 255 255 221 241 255 195 136 254 255 255 227 109 216 255 255 255 255 255 255 140  59  64  61  59 ; ...
	 97  98  98  98  93 141 255 255 255 255 255 255 219  92 182 255 255 255 215 224 255 249 201 253 255 255 240  83  62  62  60  57 ; ...
	 95  95  95  95  89 187 255 255 255 255 255 235 109  89 230 255 255 255 255 255 249 134 193 255 255 255 184  58  63  60  58  56 ; ...
	 92  92  92  91  97 231 255 255 255 255 246 126  77 123 253 255 255 255 255 247 131  91 246 255 255 255 114  57  61  59  56  54 ; ...
	 89  89  90  85 126 254 255 255 255 254 146  79  78 171 255 255 255 255 245 126  59 156 255 255 255 223  67  60  59  57  55  52 ; ...
	 86  87  87  80 173 255 255 255 255 169  78  82  83 219 255 255 255 241 122  65  74 220 255 255 255 149  44  51  49  46  49  51 ; ...
	 84  84  83  85 221 255 255 255 191  79  80  76 114 251 255 255 239 115  65  66 115 254 255 255 254 173 148 149 148 150  99  45 ; ...
	 81  81  77 115 255 255 255 218  86  77  79  72 131 255 255 251 119  63  70  62 133 255 255 255 255 255 255 255 255 239  86  42 ; ...
	 78  78  74 117 179 178 163  94  74  77  76  74  82 175 207 133  65  69  68  65  70 140 170 170 170 169 168 166 138  71  36  98 ; ...
	 75  75  76  73  69  69  68  73  74  74  73  73  70  65  67  64  67  66  65  64  61  54  52  50  48  47  45  43  41  43  92 241 ; ...
	 73  73  73  73  73  73  73  73  72  72  71  70  70  68  67  67  66  64  63  62  60  59  57  56  54  52  51  49  44 106 244 191 ; ...
	]/255;
MMM_icon(:,:,2) = [ ...
	191 191 191 191 191 191 254 223 167 134 124 120 117 114 111 108 105 102  99  96  93  90  89  95  92  89  86  84  81  78  76  72 ; ...
	191 191 191 191 191 237 154 127 128 126 124 121 118 115 112 109 106 103 100  97  94  84 159 239 236 235 235 235 235 235 238 131 ; ...
	191 191 191 191 241 160 134 133 130 127 124 121 118 115 112 109 106 103 100  97  92  95 229 255 255 255 255 255 255 255 205  63 ; ...
	191 191 191 243 201 157 135 133 130 127 124 121 118 115 112 109 106 103 100  97  92 106 175 174 172 166 221 253 254 239  90  56 ; ...
	191 191 241 164 158 140 135 132 129 126 124 121 118 115 111 109 106 103  99  96  93  90  81  77  67 115 240 254 251 135  56  60 ; ...
	191 191 159 134 136 135 133 131 128 125 123 120 117 114 111 108 105 102  99  96  93  90  87  81 109 235 255 254 233  88  58  60 ; ...
	191 229 128 134 133 132 131 129 126 124 122 119 116 113 110 107 104 101  98  95  93  90  86  82 177 253 255 252 255 224  83  56 ; ...
	191 173 127 131 131 130 128 126 124 122 120 117 114 112 109 106 103 100  97  95  92  89  86  83  80 116 225 254 252 255 182  54 ; ...
	236 136 127 127 127 127 125 124 122 120 118 115 113 110 108 105 102  99  96  94  91  83  78  75  77  68 111 247 253 255 232  68 ; ...
	184 125 124 124 124 124 123 122 120 118 116 113 111 109 106 104 101  98  95  93  86 128 177 156  89  72  76 226 255 255 230  66 ; ...
	126 121 121 121 120 120 119 118 116 114 112 110 109 107 105 102  99  97  94  87 137 250 255 255 202  72  86 239 254 255 185  52 ; ...
	117 118 117 121 124 124 123 122 120 119 117 115 106 105 103 100  98  95  93  83 188 255 252 254 244  94 171 255 255 233  88  52 ; ...
	114 114 119 212 235 232 232 232 232 231 231 233 193 105 100  98  96  93  90  82 142 255 255 255 253 230 255 251 201  92  52  55 ; ...
	111 106 168 255 252 252 252 252 251 251 252 254 240 112  97  96  92  88  97  93  79 138 194 199 198 186 151 103  62  55  57  54 ; ...
	105 121 234 252 249 250 250 249 249 249 249 254 211  97  96  91 104 190 230 230 130  70  72  70  67  63  59  60  61  59  56  53 ; ...
	104 115 139 137 137 137 134 153 246 250 249 255 165  90  91  99 213 255 252 255 194  73  76  73  70  65  60  61  60  57  54  52 ; ...
	103 102  99 100 100  99  93 157 254 249 250 246 121  89  92 198 255 249 249 255 162  70  74  70  65  91 104  71  56  56  53  50 ; ...
	 99 100 100 100 100 100  96 204 254 249 253 217  92  87 183 255 249 249 250 247 108  71  70  66 164 244 253 211  72  52  52  49 ; ...
	 96  97  97  97  97  94 111 239 251 249 255 171  78 167 255 250 249 249 253 216  77  69  68 170 255 252 251 255 117  46  50  48 ; ...
	 93  94  94  94  94  88 146 254 249 250 249 116 147 252 248 240 250 249 255 168  62  73 182 255 251 250 251 246  87  47  49  46 ; ...
	 91  91  91  91  91  86 192 255 249 253 216 140 247 255 179 226 251 249 251 110  75 198 255 250 250 250 255 201  50  49  47  44 ; ...
	 87  88  88  88  86  97 232 252 249 252 214 234 255 186 125 251 249 253 220  96 210 255 251 249 250 250 255 129  43  48  45  43 ; ...
	 84  85  85  85  80 130 252 249 249 249 250 255 212  79 173 255 249 253 207 217 255 244 193 250 250 253 235  68  46  46  44  41 ; ...
	 81  82  82  82  75 178 255 249 249 249 253 230  96  75 222 253 249 249 249 255 246 122 184 255 250 255 175  42  47  44  42  39 ; ...
	 79  79  79  78  84 223 253 249 249 252 242 114  63 111 249 250 249 249 253 244 119  77 242 252 251 252 101  41  45  42  40  37 ; ...
	 76  76  76  71 115 250 250 249 250 251 135  64  63 161 255 249 249 253 241 113  43 146 255 250 254 216  51  44  43  40  38  36 ; ...
	 72  73  73  65 163 255 249 249 255 158  63  68  68 211 254 249 253 237 109  49  59 213 254 250 255 138  27  34  31  29  32  34 ; ...
	 69  69  69  70 213 254 249 255 182  65  66  62 101 247 250 252 235 102  49  50 102 251 251 250 250 164 136 138 137 142  85  27 ; ...
	 66  66  62 102 253 255 255 211  72  62  64  57 119 255 255 248 106  48  55  46 121 255 255 255 255 255 255 255 255 235  71  25 ; ...
	 63  64  60 104 172 170 154  80  59  62  61  59  67 165 202 121  49  53  52  49  55 129 162 161 160 159 158 157 127  55  18  85 ; ...
	 60  60  61  58  53  53  53  57  59  59  58  57  55  49  51  48  52  51  49  48  45  37  35  33  31  29  27  25  23  26  79 240 ; ...
	 58  58  58  58  58  58  58  57  57  56  56  55  54  53  51  51  50  49  47  46  44  43  41  39  37  36  34  32  27  94 243 191 ; ...
	]/255;
MMM_icon(:,:,3) = [ ...
	191 191 191 191 191 191 254 236 200 180 175 172 170 168 166 164 162 160 159 157 155 153 152 152 151 149 147 146 144 142 140 138 ; ...
	191 191 191 191 191 244 192 176 176 176 174 172 170 168 166 164 163 161 159 157 155 151 172 197 195 195 195 194 194 194 196 158 ; ...
	191 191 191 191 246 197 180 180 178 176 174 172 170 168 166 165 163 161 159 157 155 155 193 204 203 203 202 201 201 203 184 135 ; ...
	191 191 191 248 221 194 181 180 178 176 174 172 170 168 166 165 163 161 159 157 155 157 177 176 176 173 190 201 201 196 145 133 ; ...
	191 191 247 199 195 184 181 179 177 175 174 172 170 168 166 164 163 161 159 157 155 153 149 147 143 157 197 201 200 161 134 134 ; ...
	191 191 196 180 181 181 180 178 177 175 173 171 169 168 166 164 162 160 158 157 155 153 151 148 156 195 201 201 194 145 134 134 ; ...
	191 239 176 180 180 179 178 177 175 174 172 171 169 167 165 164 162 160 158 156 154 153 151 148 177 201 201 200 201 191 143 133 ; ...
	191 204 176 178 178 177 177 175 174 173 171 170 168 166 165 163 161 159 158 156 154 152 150 148 147 157 192 201 200 203 176 132 ; ...
	243 181 176 176 176 176 175 174 173 172 170 168 167 165 164 162 160 159 157 155 153 150 148 146 145 141 155 199 200 201 193 137 ; ...
	211 175 174 174 174 174 173 172 171 170 169 167 166 164 163 161 159 158 156 154 152 163 178 171 149 143 143 191 201 201 193 136 ; ...
	175 172 172 173 173 172 172 171 170 169 168 166 165 163 162 160 159 157 155 152 166 200 202 203 184 142 146 196 201 203 177 131 ; ...
	169 170 171 170 168 168 168 167 166 166 164 163 164 162 161 159 158 156 155 151 181 203 200 200 198 149 174 202 203 194 144 131 ; ...
	168 168 167 132 125 127 129 130 131 133 134 134 144 161 159 158 157 155 154 150 167 202 202 202 200 193 203 200 183 145 131 131 ; ...
	166 168 146 116 119 120 122 123 125 127 128 129 133 157 158 157 156 154 152 150 149 165 182 183 183 179 167 150 136 133 133 131 ; ...
	165 160 122 117 120 121 123 124 126 128 129 130 139 158 157 156 153 145 143 144 147 145 143 142 140 138 136 135 135 134 132 130 ; ...
	163 160 152 153 153 154 154 149 127 127 129 130 146 158 156 153 142 139 141 143 145 146 144 142 141 139 137 136 134 133 131 129 ; ...
	160 161 162 162 162 162 163 147 125 128 129 131 152 156 154 142 138 140 141 143 145 145 143 141 140 140 140 136 133 132 130 128 ; ...
	159 159 159 159 159 159 159 135 125 128 128 136 155 154 143 136 138 140 141 143 145 143 142 140 145 152 154 151 134 131 129 127 ; ...
	157 157 157 157 157 158 154 127 126 128 128 142 155 143 134 137 138 140 141 143 144 142 140 145 151 152 154 156 139 129 128 127 ; ...
	155 155 155 155 155 156 144 123 126 127 129 148 145 133 135 138 138 140 141 143 142 141 145 149 151 152 154 155 135 128 127 126 ; ...
	153 153 154 154 154 154 134 123 126 127 133 144 132 133 140 138 138 140 141 142 141 144 148 149 151 152 154 149 129 128 126 125 ; ...
	151 151 152 152 152 150 126 124 126 127 133 132 131 138 144 136 138 140 141 141 144 146 148 149 151 152 155 139 127 127 125 124 ; ...
	149 150 150 150 151 142 122 124 126 128 129 130 135 145 140 136 138 140 141 143 145 146 144 149 151 152 152 131 127 126 124 123 ; ...
	148 148 148 148 149 133 122 124 126 128 129 132 143 144 136 137 138 140 141 143 144 138 143 149 151 153 144 127 126 125 123 121 ; ...
	146 146 146 146 145 125 122 124 126 127 130 140 143 140 135 137 138 140 141 143 138 135 147 149 151 153 134 125 125 123 122 120 ; ...
	144 144 144 145 139 121 123 124 126 127 137 142 141 137 135 137 138 140 141 137 133 139 148 149 151 148 127 125 124 122 121 119 ; ...
	142 142 142 143 131 120 123 124 126 134 141 140 139 135 135 137 138 139 136 133 132 143 148 149 152 137 123 122 121 118 118 118 ; ...
	140 140 140 140 124 121 123 124 131 139 139 138 137 134 135 137 138 135 132 131 134 146 148 149 151 141 137 138 137 147 136 115 ; ...
	138 138 139 135 119 120 122 127 137 137 137 136 135 133 135 137 134 132 131 130 135 147 149 150 152 154 157 158 160 168 130 114 ; ...
	136 136 137 133 126 127 129 135 136 135 135 134 134 134 134 133 131 130 129 128 128 134 138 138 139 139 139 140 135 122 110 149 ; ...
	134 134 134 135 135 135 135 134 134 134 133 133 132 131 131 130 129 128 128 127 126 124 123 122 120 119 118 117 115 115 146 246 ; ...
	133 133 133 133 133 133 133 133 133 132 132 131 131 130 129 129 128 127 126 125 125 124 122 122 120 119 118 117 114 156 248 191 ; ...
	]/255;

PDBwiki_icon(:,:,1) = [ ...
	  0   0   0   0 191 NaN   0   0   0   0 NaN NaN   0   0 NaN NaN ; ...
	  0 191 191 191   0 191   0 191 NaN   0   0 NaN   0 NaN   0 NaN ; ...
	  0 191 NaN NaN   0 191   0 191 NaN NaN   0 191   0 NaN NaN   0 ; ...
	  0 191 NaN NaN   0 191   0 191 NaN NaN   0 191   0 NaN NaN   0 ; ...
	  0   0   0   0 191 191   0 191 NaN NaN   0 191   0   0   0 NaN ; ...
	  0   0 191 191 191 NaN   0 191 NaN NaN   0 191   0 NaN NaN   0 ; ...
	  0 191 NaN NaN NaN NaN   0 191 NaN NaN   0 NaN   0 NaN NaN   0 ; ...
	  0 191 NaN NaN NaN NaN   0 191 NaN   0   0 NaN   0 NaN   0 NaN ; ...
	  0 191 NaN NaN NaN NaN   0   0   0   0 NaN NaN   0   0 NaN NaN ; ...
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
	NaN 191 191 191 191 191 NaN 255 NaN 255 NaN NaN NaN 255 NaN NaN ; ...
	NaN 191 NaN 191 NaN 191 NaN NaN NaN 255 NaN NaN NaN NaN NaN NaN ; ...
	NaN 255 191 255 191 255 NaN 255 NaN 255 NaN 255 NaN 255 NaN NaN ; ...
	NaN 255 191 255 191 255 NaN 255 NaN 255 255 NaN NaN 255 NaN NaN ; ...
	NaN NaN 255 NaN 255 NaN NaN 255 NaN 255 NaN 255 NaN 255 NaN NaN ; ...
	NaN NaN 255 NaN 255 NaN NaN 255 NaN 255 NaN 255 NaN 255 NaN NaN ; ...
	]/255;
PDBwiki_icon(:,:,2) = [ ...
	  0   0   0   0 191 NaN   0   0   0   0 NaN NaN   0   0 NaN NaN ; ...
	  0 191 191 191   0 191   0 191 NaN   0   0 NaN   0 NaN   0 NaN ; ...
	  0 191 NaN NaN   0 191   0 191 NaN NaN   0 191   0 NaN NaN   0 ; ...
	  0 191 NaN NaN   0 191   0 191 NaN NaN   0 191   0 NaN NaN   0 ; ...
	  0   0   0   0 191 191   0 191 NaN NaN   0 191   0   0   0 NaN ; ...
	  0   0 191 191 191 NaN   0 191 NaN NaN   0 191   0 NaN NaN   0 ; ...
	  0 191 NaN NaN NaN NaN   0 191 NaN NaN   0 NaN   0 NaN NaN   0 ; ...
	  0 191 NaN NaN NaN NaN   0 191 NaN   0   0 NaN   0 NaN   0 NaN ; ...
	  0 191 NaN NaN NaN NaN   0   0   0   0 NaN NaN   0   0 NaN NaN ; ...
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
	NaN 191 191 191 191 191 NaN   0 NaN   0 NaN NaN NaN   0 NaN NaN ; ...
	NaN 191 NaN 191 NaN 191 NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN   0 191   0 191   0 NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	NaN   0 191   0 191   0 NaN   0 NaN   0   0 NaN NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN   0 NaN NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN   0 NaN NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	]/255;
PDBwiki_icon(:,:,3) = [ ...
	128 128 128 128 191 NaN 128 128 128 128 NaN NaN 128 128 NaN NaN ; ...
	128 191 191 191 128 191 128 191 NaN 128 128 NaN 128 NaN 128 NaN ; ...
	128 191 NaN NaN 128 191 128 191 NaN NaN 128 191 128 NaN NaN 128 ; ...
	128 191 NaN NaN 128 191 128 191 NaN NaN 128 191 128 NaN NaN 128 ; ...
	128 128 128 128 191 191 128 191 NaN NaN 128 191 128 128 128 NaN ; ...
	128 128 191 191 191 NaN 128 191 NaN NaN 128 191 128 NaN NaN 128 ; ...
	128 191 NaN NaN NaN NaN 128 191 NaN NaN 128 NaN 128 NaN NaN 128 ; ...
	128 191 NaN NaN NaN NaN 128 191 NaN 128 128 NaN 128 NaN 128 NaN ; ...
	128 191 NaN NaN NaN NaN 128 128 128 128 NaN NaN 128 128 NaN NaN ; ...
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
	NaN 191 191 191 191 191 NaN   0 NaN   0 NaN NaN NaN   0 NaN NaN ; ...
	NaN 191 NaN 191 NaN 191 NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN   0 191   0 191   0 NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	NaN   0 191   0 191   0 NaN   0 NaN   0   0 NaN NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN   0 NaN NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN   0 NaN NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	]/255;

proteopedia_icon(:,:,1) = [ ...
	NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN   0 NaN   0   0   0 NaN   0   0 NaN NaN NaN NaN ; ...
	NaN NaN NaN   0 NaN NaN   0 NaN   0   0 NaN NaN   0 NaN NaN NaN ; ...
	NaN NaN   0 NaN NaN NaN   0 NaN 255   0   0 NaN NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN NaN   0   0 255 255 255 NaN   0 NaN NaN   0 NaN ; ...
	NaN   0 NaN NaN   0 255 255 128 NaN   0 NaN   0 NaN NaN NaN   0 ; ...
	NaN   0 NaN   0 255   0 128 NaN NaN   0 NaN   0 255 255 NaN   0 ; ...
	  0   0 NaN   0 255   0 NaN 255 128   0 255 NaN   0 NaN   0 NaN ; ...
	NaN   0 NaN   0 NaN   0 255 NaN 255 255 NaN 255 NaN NaN   0 NaN ; ...
	NaN   0 NaN   0 NaN NaN   0 NaN NaN   0 NaN   0 NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN   0 NaN   0 NaN   0 NaN   0 NaN   0 NaN NaN NaN ; ...
	NaN NaN NaN   0   0 NaN NaN   0 NaN   0   0   0 NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN   0   0 NaN   0   0 NaN   0 NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
	]/255;
proteopedia_icon(:,:,2) = [ ...
	NaN NaN NaN NaN NaN NaN   0   0 128 NaN NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN   0 128   0 128   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN   0 NaN 128   0   0 NaN   0   0 NaN NaN NaN NaN ; ...
	NaN NaN NaN   0 NaN NaN   0 NaN 128   0 NaN NaN   0 NaN NaN NaN ; ...
	NaN NaN   0 NaN NaN NaN   0 NaN 128 128   0 NaN NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN NaN   0 128 128 128 128 NaN   0 NaN NaN   0 NaN ; ...
	NaN   0 NaN NaN   0 128 128   0 NaN 128 NaN   0 NaN NaN NaN   0 ; ...
	NaN   0 NaN   0 128 128   0 NaN NaN 128 NaN   0 128 128 NaN   0 ; ...
	  0   0 NaN   0 128 128 NaN 128   0 128 128 NaN   0 NaN   0 NaN ; ...
	NaN   0 NaN   0 NaN 128 128 NaN 128 128 NaN 128 NaN NaN   0 NaN ; ...
	NaN   0 NaN   0 NaN NaN 128 NaN NaN 128 NaN   0 NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN   0 NaN 128 NaN 128 NaN   0 NaN   0 NaN NaN NaN ; ...
	NaN NaN NaN   0   0 NaN NaN 128 NaN   0   0   0 NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN   0   0 NaN 128   0 NaN   0 NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
	]/255;
proteopedia_icon(:,:,3) = [ ...
	NaN NaN NaN NaN NaN NaN   0   0  64 NaN NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN   0  64 128  64   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN   0 NaN  64 128 128 NaN   0   0 NaN NaN NaN NaN ; ...
	NaN NaN NaN   0 NaN NaN 128 NaN  64 128 NaN NaN   0 NaN NaN NaN ; ...
	NaN NaN   0 NaN NaN NaN 128 NaN   0  64 128 NaN NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN NaN 128  64   0  64   0 NaN 128 NaN NaN   0 NaN ; ...
	NaN   0 NaN NaN 128   0  64   0 NaN  64 NaN 128 NaN NaN NaN   0 ; ...
	NaN   0 NaN 128   0  64   0 NaN NaN  64 NaN 128  64  64 NaN   0 ; ...
	  0   0 NaN 128   0  64 NaN   0   0  64   0 NaN 128 NaN   0 NaN ; ...
	NaN   0 NaN 128 NaN  64   0 NaN  64   0 NaN   0 NaN NaN   0 NaN ; ...
	NaN   0 NaN 128 NaN NaN  64 NaN NaN  64 NaN 128 NaN   0 NaN NaN ; ...
	NaN NaN   0 NaN 128 NaN  64 NaN  64 NaN 128 NaN   0 NaN NaN NaN ; ...
	NaN NaN NaN   0 128 NaN NaN  64 NaN 128 128   0 NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN   0 128 NaN  64 128 NaN   0 NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN   0   0 128 128   0 NaN NaN NaN NaN NaN NaN ; ...
	NaN NaN NaN NaN NaN NaN 128   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
	]/255;



% determine which third-party programs are installed and accessible
dospath=which('scwrl4.exe');
if isempty(dospath),
    third_party.scwrl4=false;
else
    third_party.scwrl4=true;
end;

dospath=which('msms.exe');
if isempty(dospath),
    third_party.msms=false;
else
    third_party.msms=true;
end;

dospath=which('dssp.exe');
if isempty(dospath),
    dospath=which('dsspcmbi.exe');
end;
if isempty(dospath),
    third_party.dssp=false;
else
    third_party.dssp=true;
end;

dospath=which([third_party.modeller_version  '.exe']);
if isempty(dospath),
    third_party.modeller=false;
else
    third_party.modeller=true;
end;
