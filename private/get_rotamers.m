function [rotamers_stats,msg]=get_rotamers(coor,site,calc_opt,rot_lib,Bfac)

% function [rotamers_stats,msg]=get_rotamers(coor,site,calc_opt,rot_lib,Bfac);

% put rotamers from the selected library onto selected mutation posisiotn
% of the protein and calculate the whole interaction energy between protein
% and rotamer atoms.

% coor          Nx4 array of xyz protein coordinates preceeded by the atom number (carbon=6, nitrogen=7 etc)
% site          structure, contains information about the requested mutation
%               position (indices of the protein atoms belonging to the
%               mutated residue etc); see line 77ff
% calc_opt      various calculation options (temperature for Boltzman factors,
%               statistics flag etc); see line 44ff for definitions and defaults
% rot_lib       structure with choosen rotamer library (atom coordinate for all
%               rotamers, useful atomic numbers, internal probabilities for each rotamer etc);
% Bfac          B factors corresponding to coor
%
% rotamer_stats
%
% msg   structure with fields .error =0: no error
%       .text clear text error or warning message
%
% Ye. Polyhach, 2009
%
%--------------------------------------------------------------------------
% bolt=1.3806503e-23;	% Boltzmann constant in CI (J/K), OLD STUFF

global label_defs

gas_un=8.314472;    % universal gas constant in CI (J/(mol*K)       

temperature_warning=0.1; % relative temperature difference from calibration at which error 1 (a warning) is set

msg.error=0;
msg.text='Rotamer analysis successful.';
rotamers_stats.all_potentials=[]; % an empty structure should be returned when an error is encountered

    
%extract information of rotamer library
library=rot_lib.library;
calibration=rot_lib.calibration;
usefull_atoms=rot_lib.usefull_atoms;
% midNO = usefull_atoms.midNO;  % indexes for N and O label atoms
sidech=usefull_atoms.side_chain; % index for the first side-chain atoms
maxdist=calibration.maxdist; % max distance rotamer atom from the origin (for the entire library)
n_rot=length(library);  % get number of rotamers in the library

% code for dealing with libraries generated before version 2014.1
labnum=tag2id(rot_lib.label,label_defs.restags);
if isfield(label_defs.residues(labnum),'class') && ~isempty(label_defs.residues(labnum).class),
    class = label_defs.residues(labnum).class; 
else
    class = 'nitroxide';
end;
if isfield(rot_lib,'spin_density') && ~ isempty(rot_lib.spin_density),
    spin_density = rot_lib.spin_density;
elseif strcmpi(class,'nitroxide'),
    spin_density = zeros(2,2);
    spin_density(:,1) = usefull_atoms.midNO;
    spin_density(:,2) = [0.5;0.5];
else
    msg.error=5;
    msg.text=sprintf('Spin density center not defined. Setting populations to zero.');
end;

% determine attachment frame
clabel = rot_lib.label;
cli = tag2id(clabel,label_defs.restags);
attach_frame = label_defs.residues(cli).res_frame;
if isempty(attach_frame),
    attach_frame = ':C:N:CA:';
end;

% extract information on calculation options

if isfield(calc_opt,'T'),
    T=calc_opt.T;	% absolute temperature (K);
else
    T=calibration.T;
end;

if isfield(calc_opt,'ext_potential'),
    switch calc_opt.ext_potential
        case 'OPLS'
            type_poten=1;
        case 'charmm27'
            type_poten=2;
        case 'avgLJ'
            type_poten=3;
            load mean_potentials
            [~,ppoi]=min(abs(forces-calc_opt.force));
            interEn_opt.p=polynomials(ppoi,:);
            interEn_opt.mu=mu_vecs(:,ppoi);
        case 'repulsion'
            type_poten=4;
        case 'RA'
            type_poten=5;
        case 'clash'
            type_poten=6;
        case 'UFF'
            type_poten=7;
        case 'UFF_Bavg'
            type_poten=8;
        case 'LJ_UFF_B'
            type_poten=9;
        case 'UFF_attract'
            type_poten=10;
    end;
else
    type_poten=2; % charmm27 as default
end;

if isfield(calc_opt,'forgive'), % forgive factor for minor clashes (reduction factor for effective van-der-Waals radius)
    forgive=calc_opt.forgive;
else
    forgive=0.5; % default forgive factor
end;

if isfield(calc_opt,'repulsion'), % repulsion exponent
    repulsion=calc_opt.repulsion;
else
    repulsion=7; % default repulsion factor
end;

if isfield(calc_opt,'pair_stats'),
    pair_stats=calc_opt.pair_stats;	% flag: 0 - no statistics made (default)
    %                                       1 - statistics for every pair of atoms is saved (results in a huge file, time consuming)
else
    pair_stats=0;
end;

% extract information on mutation site

index_array=site.res_atoms;

% ypax_name = id2tag(1,attach_frame);
% xax_name = id2tag(2,attach_frame);
% orig_name = id2tag(3,attach_frame);

xax = site.xax;
ypax = site.ypax;
orig = site.orig;

%--------------------------------------------------------------------------
n_rot=length(library);  % get number of rotamers in the library

NOcoor=zeros(n_rot,5);
mut_coor1=cell(1,n_rot);
labels_own=cell(1,n_rot);

int_pop0=calibration.pop;    % get non-normalized internal populations from the rotamers
int_pop0=int_pop0/sum(calibration.pop);   % normalize populations
int_pop=int_pop0(1:n_rot); % debugging trick: to be able to reduce n_rot easily

relative_Delta_T=T/calibration.T-1; % relative temperature difference from calibration

if abs(relative_Delta_T)>temperature_warning,
    msg.error=1;
    msg.text=sprintf('Target temperature deviates by %2.0f%% from calibration temperature. Reliable limit is %2.0f%%.',100*abs(relative_Delta_T),100*temperature_warning);
end

[m,n]=size(coor);

poten_ext=zeros(1,n_rot);
pop_rot=poten_ext;
rot_clash=cell(1,n_rot);

offset=coor(orig,2:4);	% origin of attachment frame
for kp=1:m,
    coor(kp,2:4)=coor(kp,2:4)-offset;
end;
x=coor(xax,2:4)-coor(orig,2:4); % x axis is along C_alpha-N bond for peptides
x=x/norm(x);    % unit vector along x
yp=coor(ypax,2:4)-coor(orig,2:4); % y axis is in the plane spanned by x axis and C-Ca bond for peptide
yp=yp/norm(yp);
z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z=z/norm(z);
y=cross_rowvec(z,x); % real (corrected) y axis 
dircos=[x;y;z];
Rp=dircos; % rotation matrix for conversion to standard frame
% ### Alternative code, this is problematic (mirror image)
% coor1=[coor(N,2:4);coor(Ca,2:4);coor(C,2:4)];
% coor2=[1.462,0,0;0,0,0;-0.562,1.417,0];
% [rmsd,coor2b,transmat]=rmsd_superimpose(coor1,coor2);
% Rp=transmat(1:3,1:3);
% ### end alternative code
%-----
% prepare protein atoms for the analysis:
% atoms of the labelled residue must be excluded
all_indices = 1:m;
if ~isempty(index_array),
    for kk = 1:length(index_array),
        all_indices(index_array(kk)) = 0;
    end;
end;
coor0ind = all_indices(all_indices~=0);
coor0 = coor(coor0ind,:);
Bfac_protein = Bfac(coor0ind);
%-----
% from all protein atoms select those who are within a certain sphere
% around the attachment point (Calpha)
coor10=zeros(length(coor0ind),4);
Bfac10=zeros(length(coor0ind),1);
coor0ind_orig=zeros(1,length(coor0ind));
kp1=0;
for kp=1:length(coor0ind)
    inter_atom=(sum((coor0(kp,2:4)-coor(orig,2:4)).^2))^(1/2);
    if inter_atom<(maxdist+4)  % 5 Angstroem - is absolutely safe.
% 4 Angstroem: 10% faster, and practically, no difference in the rot populations. (part. function changes at most 0.5% at some checked positions)
% 3 Angstroem: 20% faster, saw change of about 1.8% in the partition funciton at some positions
        kp1=kp1+1;
        coor10(kp1,:)=coor0(kp,:);
        Bfac10(kp1)=Bfac_protein(kp);
        coor0ind_orig(kp1)=coor0ind(kp); % store original indices
    end
end
clear coor0;
coor0=coor10(1:kp1,:);
coor0ind_orig=coor0ind_orig(1:kp1);
Bfac_protein = Bfac10(1:kp1);
%-----
for k=1:n_rot % number of rotamers in the library
    
    lcoor=library(k).ecoor;
    Bfac_label = calibration.Bfactors(:,k);
    new_lcoor=lcoor;
    rot_lcoor=lcoor;
    own_lcoor=lcoor;
    [ml,nl]=size(lcoor);
    for kr=1:ml % conversion to standard frame
        vec=lcoor(kr,2:4);
        newvec=vec*Rp; 
        new_lcoor(kr,2:4)=newvec+coor(orig,2:4);
        rot_lcoor(kr,2:4)=newvec+offset;
        own_lcoor(kr,2:4)=newvec;        % labels are stored in their own residue frame (Ca is always [0,0,0]);
    end;
    
    mut_coor1{k}=rot_lcoor;	% stores coordinates for running rotamer (k)
    labels_own{k}=own_lcoor;
    
    % prepare rotamer coordinates for analysis:
    % only side-chain has to be considered during LJenergy calculation
    lcoor0=new_lcoor(sidech:ml,:);
    lcoor0ind=(sidech:ml);
    Bfac_label=Bfac_label(sidech:ml);

    interEn_opt.energy_stats=pair_stats;
    interEn_opt.forgive=forgive;
    interEn_opt.repulsion=repulsion;
    if pair_stats==1 % !!! must be consistant: or_ind1 must correspond to coor1!!!
        interEn_opt.ind_coor1=lcoor0ind;   % (see get_energyLJ_OPLS for that)
        interEn_opt.ind_coor2=coor0ind_orig;
%         rot_clash=cell(1,n_rot);
    end

    switch type_poten
        case 1
            interEn=get_energyLJ_OPLS(lcoor0,coor0,interEn_opt);
            type_poten_string='OPLS parametrized Lennard-Jones';
            ext_en=interEn.energy;
        case 2
            interEn=get_energyLJ_charmm27(lcoor0,coor0,interEn_opt);
            type_poten_string='charmm27 parametrized Lennard-Jones';
            ext_en=interEn.energy;
        case 3
            interEn=get_energy_avgLJ(lcoor0,coor0,interEn_opt);
            type_poten_string='harmonic average Lennard-Jones';
            ext_en=interEn.energy;
        case 4
            interEn=get_energy_repulsion(lcoor0,coor0,interEn_opt);
            type_poten_string='pure repulsion';
            ext_en=interEn.energy;
        case 5
            interEn=get_energyRA_charmm27(lcoor0,coor0,interEn_opt);
            type_poten_string='charmm27 parametrized repulsion-attraction';
            ext_en=interEn.energy;
        case 6
            interEn=get_energy_clash(lcoor0,coor0,interEn_opt);
            type_poten_string='clash rejection';
            ext_en=interEn.energy;
        case 7
            interEn=get_energyLJ_UFF(lcoor0,coor0,interEn_opt);
            type_poten_string='UFF parametrized Lennard-Jones';
            ext_en=interEn.energy;
        case 8
            interEn=get_energyLJ_UFF_Bavg(lcoor0,coor0,interEn_opt,Bfac_label,Bfac_protein);
            type_poten_string='UFF parametrized Lennard-Jones B factor averaged';
            ext_en=interEn.energy;
        case 9
            interEn=get_energyLJ_UFF_B(lcoor0,coor0,interEn_opt,Bfac_label,Bfac_protein);
            type_poten_string='UFF parametrized Lennard-Jones B factor scaled';
            ext_en=interEn.energy;
        case 10
            interEn=get_energyLJ_UFF_attract(lcoor0,coor0,interEn_opt);
            type_poten_string='UFF parametrized Lennard-Jones with enhanced attaction';
            ext_en=interEn.energy;
        otherwise 
            msg.error=2;
            msg.text='No external potential corresponds to the ext_potential flag used! Reverting to charmm27.';
            interEn=get_energyLJ_charmm27(lcoor0,coor0,interEn_opt);
            type_poten_string='charmm27 parametrized Lennard-Jones';
            ext_en=interEn.energy;
    end
    poten_ext(k)=exp(-ext_en/(gas_un*T));  % Boltzmann factor for LJ potential for k-th rotamer
    if pair_stats==1
        rot_clash{1,k}=LJenergy.clash_stats;
    end
end



poten_ext_sum=sum(poten_ext);   % net Boltzmann factor for the mutation position (partition sum)
partition_function=sum(poten_ext.*int_pop)/sum(int_pop); % partition function relative to the one of the free label
if poten_ext_sum==0 && pair_stats==0
    msg.error=3;
    msg.text=sprintf('External energies are too high. External populations set to zero.');
%     msg.text=sprintf('%s Try to rerun with pair_stats=1 to analyze problem',msg.text); 
    ext_pop=zeros(1,n_rot);
elseif poten_ext_sum==0 && pair_stats==1
    msg.error=4;
    msg.text=sprintf('External energies are too high. External populations set to zero.');
%     msg.text=sprintf('%s Clash statistics was saved.',msg.text); 
    ext_pop=zeros(n_rot);
end

if poten_ext_sum~=0 && pair_stats==0
    ext_pop=poten_ext/poten_ext_sum;
    pop_rot=ext_pop.*int_pop; % get populations for all rotamers
    pop_rot=pop_rot/sum(pop_rot); % normalize full rotamer populations
elseif poten_ext_sum~=0 && pair_stats==1
    ext_pop=poten_ext/poten_ext_sum;
    pop_rot=ext_pop.*int_pop; % get populations for all rotamers
    pop_rot=pop_rot/sum(pop_rot); % normalize full rotamer populations
end

% at that point variables mut_coor1 and rot_coef contain all rotamer
% coordinates and potentials (respectively) in the protein coordinate
% frame. The following creates arrays of the N-O centers coordinates vs
% rotamer potentials for chosen residue.

%keyboard
% old-style spin center coordinate
% lN=midNO(1);
% lO=midNO(2);
% for k=1:n_rot
%     mut1=mut_coor1{k}; % k-th rotamer in the mutation pos. 1
%     NOmut1=(mut1(lN,2:4)+mut1(lO,2:4))./2; % "NO" coordinates for the k-th rotamer
%     NOcoor(k,1:3)=NOmut1; 
%     NOcoor(k,4)=pop_rot(k);
%     NOcoor(k,5)=k;
% end
spinindices = spin_density(:,1);
spindens = spin_density(:,2);
spindens = spindens/sum(spindens);
for k=1:n_rot
    mut1=mut_coor1{k}; % k-th rotamer in the mutation pos. 1
    spincent = mut1(spinindices,2:4);
    NOmut1 = spindens'*spincent; % mean spin center coordinates for the k-th rotamer
    NOcoor(k,1:3)=NOmut1; 
    NOcoor(k,4)=pop_rot(k);
    NOcoor(k,5)=k;
end
NOstats_all=NOcoor;

% all potentials in one structure:
all_pops.pop_rot=pop_rot; % net population: ext*int for each rotamer
all_pops.ext_net=poten_ext_sum; % net Boltzmann factor for external potential (partition sum)
all_pops.ext_pop=ext_pop; % populations based on the ext potentials only
all_pops.partition_function=partition_function; % full partition function
% statistics for the chosen mutation position is saved as a structure
% rotamer_stats with the following fields:
rotamers_stats.NOall=NOstats_all; % NO-centers and weights for all rotamers in a single matrix;
rotamers_stats.all_rotamers=mut_coor1; % contains coordinates for all rotated rotamers for current position
rotamers_stats.all_potentials=all_pops; % total potential for the mutation position;
rotamers_stats.labels_own=labels_own; % rotamer coordinates in the local residue frame (Ca is always [0,0,0]);
rotamers_stats.loc_frame_Ca=offset; % global (protein) coordinates of the local frame origin (Ca alpha of the mutated residue for peptides);
rotamers_stats.ext_poten_type=type_poten_string; % type of the interatomic potential used
if pair_stats==1
    rotamers_stats.rot_clash=rot_clash;  % save clash statistics if needed
end
