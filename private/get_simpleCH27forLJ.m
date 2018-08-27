function [r_vdw,q,Rmin2,eps]=get_simpleCH27forLJ(index)

% (for the same purpose previously (before 0109) funtion poten_param.m was used)

% returns SIMPLIFIED charmm27 parameters for the Lennard-Jones potential
% (see homepage for CHARMM or Alex MacKernell homepage;
% SIMPLE: because there is no differentiation between different charmm
% types of the same atom. For instance, for all different carbons from
% charmm  - the same set of the force field parameetrs is assigned;
% to eliminate hydrogens - epsilon=0 is ssigned for all H.
% (installation packege fro CHARMM toppar_c35b2_c36a2 was downloaded from the CHARMM homepage at 150609)

% CHARMM van der Waals term (see toppar_c35b2_c36a2):
%   V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
%   epsilon: kcal/mol, Eps,i,j = sqrt(eps,i * eps,j)
%   Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j

% index - number of the atom in the structure. (corresponds to the atom number in the Periodic Table)

% Van der Waals radii (A)   (http://www.webelements.com/ & Gunnar's full_vdw.mat - they match for the specified elements)
% (this additional value for the VdW radius is returned for convinience here!)

% @160609NOTE: zero is assigned for all atoms just because it is not in use
% presently. This means that undefined atom types do not contribute to
% clash energy, which often makes sense, as they were probably added for phasing only  
%--------------------------------------------------------------------------

vdw=zeros(1,92);
p_q=zeros(1,92);
p_Rmin2=zeros(1,92);
p_eps=zeros(1,92);

% Hydrogen, 1
vdw(1)=1.20; % VdW radius
p_q(1)=0; p_Rmin2(1)=0.0; p_eps(1)=0.0;

% Carbon, 6
vdw(6)=1.70; % VdW radius
p_q(6)=0; p_Rmin2(6)=2.275; p_eps(6)=-0.020; % CT1 "peptide Alpha carbon"
% p_q(6)=0; p_Rmin2(6)=1.850; p_eps(6)=-0.200; % CT2 "GLY Alpha carbon"

% Nitrogen, 7
vdw(7)=1.55; % VdW radius
p_q(7)=0; p_Rmin2(7)=1.850; p_eps(7)=-0.200; % NH1 "Peptide Nitrogen"

% Oxygen, 8
vdw(8)=1.52; % VdW radius
p_q(8)=0; p_Rmin2(8)=1.700; p_eps(8)=-0.120; % O "Peptide Oxygen"

% Sodium, 11
vdw(11)=1.66;
p_q(11)=0; p_Rmin2(12)=1.36475; p_eps(12)=-0.0469; % SOD "Sodium Ion"

% Magnesium, 12
vdw(12)=1.73; % VdW radius
p_q(12)=0; p_Rmin2(12)=1.185; p_eps(12)=-0.015; % MG "Magnesium Ion"

% Sulphur, 16
vdw(16)=1.80; % VdW radius
p_q(16)=0; p_Rmin2(16)=2.000; p_eps(16)=-0.450; % S "Sulfide Sulfur"
% p_q(16)=0; p_Rmin2(16)=1.975; p_eps(16)=-0.380; % SM "Disulfide Sulfur"

% Chloride, 17
vdw(17)=1.75;
p_q(17)=0; p_Rmin2(17)=2.27; p_eps(19)=-0.150;

% Phosphorus, 15 (taken EQUAL to Sulfur!)
vdw(15)=1.80; % VdW radius
p_q(15)=0; p_Rmin2(15)=p_Rmin2(16); p_eps(15)=p_eps(16);

% Potassium, 19
vdw(19)=2.75;
p_q(19)=0; p_Rmin2(19)=1.76375; p_eps(19)=-0.0870;

% Calcium, 20
vdw(20)=1.80;
p_q(20)=0; p_Rmin2(20)=1.367; p_eps(20)=-0.120;

% Zinc, 30
vdw(30)=1.39; % VdW radius
p_q(30)=0; p_Rmin2(26)=1.090; p_eps(26)=-0.25; % ZN

% Cesium, 55
vdw(20)=1.80;
p_q(20)=0; p_Rmin2(20)=2.100; p_eps(20)=-0.190;

% all first- and second-row transition metal ions are treated like
% zinc
for k=21:30,
    p_q(k)=0; p_Rmin2(k)=p_Rmin2(30); p_eps(k)=p_eps(30);
end;

for k=39:48,
    p_q(k)=0; p_Rmin2(k)=p_Rmin2(30); p_eps(k)=p_eps(30);
end;

% Selen, 34
% we use sulfur values here, assuming that the protein for spin-label
% experiments is made with native amino acids
vdw(34)=1.80; % VdW radius
p_q(34)=0; p_Rmin2(34)=2.000; p_eps(34)=-0.450; % S "Sulfide Sulfur"
% p_q(16)=0; p_Rmin2(16)=1.975; p_eps(16)=-0.380; % SM "Disulfide Sulfur"

conv_factor=(4.185*1000); % conversion to CI, - if working per mol

if index<=length(vdw),
    r_vdw=vdw(index); %VdW radius of index-th atom in a given structure
    q=p_q(index); 
    Rmin2=p_Rmin2(index); 
    eps=p_eps(index)*conv_factor;
else % all unknown atoms are ignored
    r_vdw=0;
    q=0;
    Rmin2=0;
    eps=0;
end;