fname='IA1_298K_090814';
load(fname);

if exist('library')~=2,
    library=rot_lib.library;
    usefull_atoms=rot_lib.usefull_atoms;
    connect=rot_lib.connect;
    stdframe=rot_lib.stdframe;
    calibration=rot_lib.calibration;
    chi_select=rot_lib.chi_select;
    forcefield=rot_lib.forcefield;
    types=rot_lib.types;
    chi_def=rot_lib.chi_def;
end;

label='IA1';
MDtags=':N:CA:C:O:CB:SG:CD:CE:OE:NZ:C3:C4:C5:N1:C2:C6:C7:O1:C8:C9:H:H2:HA:OXT:';
MDtags=[MDtags 'HB2:HB3:HD2:HD3:HZ:H31:H41:H42:H61:H62:H63:H71:H72:H73:H81:H82:H83:H91:H92:H93:'];
load labels
lid=tag2id(label,label_defs.restags);
labeldef=label_defs.residues(lid);
atoms=length(labeldef.elements);
[atoms2,n]=size(library(1).ecoor);
PDB2MD=zeros(1,atoms2);
for k=1:atoms2,
    atag=id2tag(k,MDtags);
    aid=tag2id(atag,labeldef.atoms);
    if isempty(aid), aid=0; end;
    PDB2MD(k)=aid;
end;
MD2PDB=zeros(1,atoms);
for k=1:atoms,
    atag=id2tag(k,labeldef.atoms);
    aid=tag2id(atag,MDtags);
    if isempty(aid), aid=0; end;
    MD2PDB(k)=aid;
end;

Ca0=2;
for il=1:length(library)
    ecoor0=library(il).ecoor;
    ecoor=ecoor0;
    rot_zero=ecoor0(Ca0,2:4); % get Ca coordinates for each rotamer
    for iil=1:length(ecoor0(:,1))
        ecoor(iil,2:4)=ecoor0(iil,2:4)-rot_zero;
    end
    library(il).ecoor=ecoor;
end


rot_lib.library=library;
rot_lib.usefull_atoms=usefull_atoms;
rot_lib.connect=connect;
rot_lib.stdframe=stdframe;
rot_lib.calibration=calibration;
rot_lib.chi_select=chi_select;
rot_lib.forcefield=forcefield;
rot_lib.types=types;
rot_lib.chi_def=chi_def;
rot_lib.PDB2MD=PDB2MD;
rot_lib.MD2PDB=MD2PDB;
rot_lib.label=label;

% save(fname,'chi_select','library','stdframe','calibration','connect','forcefield','types','chi_def','usefull_atoms','PDB2MD','MD2PDB','label');

save(fname,'rot_lib');
