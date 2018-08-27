function calc_positions=rotamer_populations(allindices,rot_lib,T,no_compute,stat_file,PDB_path,library,no_context)
% Computes populations of rotamers for a set of residues defined by index
% array allindices, using the rotamer library rot_lib at temperature T,
% T defaults to the calibration temperature of the rotamer library (as
% coded in get_rotamers)
% if allindices contains references to objects other than residues, these
% are ignored
% populations are determined independently, i.e. mutual clashes between
% labels are not considered
%
% no_compute        optional flag that requests to return only pre-computed
%                   rotamers, defaults to false
% stat_file         optional file handle for statistics file
% PDB_path          optional path for individual rotamer PDB files
% library           optional name of rotamer library
% no_context        optional flag, if true, the mean coordinates of
%                   non-interacting labels are returned with associated
%                   weights 1, rotamer structures are empty, potential is 1
%                   in this case, defaults to false
%
% calc_positions    vector of structures (length corresponds to number of
%                   residues in allindices minus number of residues that 
%                   could not be labelled) with fields:
%                   .label  three-letter label code
%                   .rotamers() vector of structures for all significant
%                   rotamers with fields
%                   .indices    residue indices into MMM object structure
%                   .coor       coordinate array
%                   .pop        population vector
%                   .NOpos      NO group coordinates (columns 1-3) with 
%                               populations (column 4)
%                   .T          labeling temperature
%                   .potentials full potential vector
%

global hMain
global model
global third_party
global label_defs

% Determine attachment frame

clabel = rot_lib.label;
cli = tag2id(clabel,label_defs.restags);
attach_frame = label_defs.residues(cli).res_frame;
if isempty(attach_frame)
    attach_frame = ':C:N:CA:';
end


% The following is for testing averaged Lennard-Jones potentials
% initial tests did not show any advantage, GUJE, 3.1.2014
% global avg_LJ_pot
% 
% load avg_LJ_potentials_175K_100cal
% [rm,sigrm] = meshgrid(avg_LJ_pot.r,avg_LJ_pot.sigr);
% avg_LJ_pot.rm = rm;
% avg_LJ_pot.sigrm = sigrm;

if ~exist('no_context','var')
    no_context = false;
end;

if nargin<4 || isempty(no_compute), no_compute=false; end;

if nargin<5 || isempty(stat_file) || stat_file==-1, hMain.statistics=0; end;

if nargin<6 || isempty(PDB_path), PDB_path=''; end;
    
if nargin<7 || isempty(library), library='UNKNOWN'; end;

threshold=0.005; % 0.5% of population is neglected

consider_water=false; % flag for also considering water collisions, should 
                      % be true only for comparison with old computations
                      % or for testing
                      
packing_radius=20; % radius in which sidechains are repacked for dynamic computations

forgive= 0.5; % 0.5 for LJ, charmm27, 0.35 for pure repulsion, 0.25 for clash

potential_type='LJ';

% for UFF-derived libraries, set optimal forgive of 0.50 
if strcmpi(rot_lib.calibration.force_field,'UFF_Towhee')
    forgive = 0.50; % 0.50
    potential_type='LJ_UFF';
    % potential_type='repulsion';
end

if isfield(rot_lib,'forgive')
    forgive=rot_lib.forgive;
    if length(forgive) > 1
        potential_type = 'LJ_UFF_attract';
    end
end

repulsion=6.25;

msg.error=0;
msg.text='';

% hMain.dynamic_rotamers=true;

switch potential_type
    case 'LJ'
        calc_opt.ext_potential='charmm27'; % default are charmm27 Lennard-Jones energies
        calc_opt.forgive=forgive;
    case 'LJ_UFF'
        calc_opt.ext_potential='UFF'; % UFF Lennard-Jones energies
        calc_opt.forgive=forgive;
    case 'RA'
        calc_opt.ext_potential='RA'; % repulsion/attraction
        calc_opt.forgive=forgive;
        calc_opt.repulsion=repulsion;
    case 'LJ_UFF_B'
        calc_opt.ext_potential='LJ_UFF_B'; 
        calc_opt.forgive=forgive;
    case 'avgLJ'
        calc_opt.ext_potential='avgLJ'; 
        calc_opt.forgive=forgive;
        calc_opt.force=1000;
    case 'Bavg_LJ_UFF'
        calc_opt.ext_potential='UFF_Bavg'; % UFF Lennard-Jones energies B factor averaged
        calc_opt.forgive=forgive;
    case 'LJ_UFF_attract'
        calc_opt.ext_potential = 'UFF_attract';
        calc_opt.forgive = forgive;
        calc_opt.enhance_attraction = forgive(2);
    case 'repulsion'
        calc_opt.ext_potential='repulsion'; 
        calc_opt.forgive=forgive;
    case 'clash'
        calc_opt.ext_potential='clash'; 
        calc_opt.forgive=forgive;
end;    
if nargin>2 && T > 0,
    calc_opt.T=T; % set target temperature, if given
else
    T=rot_lib.calibration.T;
    calc_opt.T=T;
end;

calc_positions=[]; % initialize empty output array

stags=model.structure_tags;
stags=stags(2:end); % remove leading delimiter
nonsense=textscan(stags,'%s','Delimiter',':');
stag_list=nonsense{1};

respoi=0; % pointer for successfully labeled residues
[mres,nres]=size(allindices);
add_msg_board(sprintf('Rotamer computation for %i selected object(s).',mres));
msgpoi=1;

if mres>5,
    comp_status=status_figure('Rotamers: Close to stop.');
end;

for kres=1:mres, % loop over all objects
    if mres>5 && mod(kres,round(mres/10))==0    % just a counter
        comp_status=status_figure((kres-1)/mres,comp_status);
        if ~comp_status, 
            calc_positions=[];
            return 
        end;
    end;
    cindices=allindices(kres,:);
    cindices=cindices(cindices>0); % cut zero indices
    if length(cindices)==4, % consider only residue indices
        computed=0;
        if isfield(model,'sites'),
            for ks=1:length(model.sites),
                for kc=1:length(model.sites{ks}),
                    for kr=1:length(model.sites{ks}(kc).residue),
                        eindices=model.sites{ks}(kc).residue(kr).indices;
                        if cindices==eindices,
                            % this residue was already computed in a site scan
                            if abs(T-model.sites{ks}(kc).residue(kr).T)<1 && strcmp(model.sites{ks}(kc).residue(kr).label,rot_lib.label),
                                computed=1;
                                my_rotamers=model.sites{ks}(kc).residue(kr);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        if ~computed, % compute only residues that were not previously computed in a site scan
            name=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).name;
            lid = tag2id(name,label_defs.restags);
            if ~isempty(lid),
                resadr=mk_address(cindices);
                add_msg_board(sprintf('WARNING: Residue %s was already labeled with %s',resadr,name));
            end;
            if no_compute,
                resadr=mk_address(cindices,1);
                potentials=[];
                msg.error=13;
                msg.text=sprintf('Warning: Rotamers of %s were not previously computed. Please use EPR/Site scan menu item for rotamer computation.',resadr);
            else
                resadr=mk_address(cindices);
                % get extended coordinates of the whole structure where the labeled
                % residue is
                
                if isfield(rot_lib,'noclash'),
                    noclash = rot_lib.noclash;
                else
                    noclash = '*';
                end;
                if hMain.dynamic_rotamers,
                    [bcoor,site,Bfac]=ecoor_structure(cindices,consider_water,false,attach_frame,noclash); % backbone only
                    if isempty(bcoor), % site preparation failed (no amino acid or crucial backbone atoms missing)
                        continue;
                    end;
                else
                    [coor,site,Bfac]=ecoor_structure(cindices,consider_water,true,attach_frame,noclash);
                    if isempty(coor), % site preparation failed (no amino acid or crucial backbone atoms missing)
                        continue;
                    end;
                end;

                add_msg_board(sprintf('Computing rotamer populations for residue %s',mk_address(cindices,1)));
                add_msg_board(sprintf('with label %s at a temperature of %4.0f K',rot_lib.label,T));
                set(hMain.figure,'Pointer','watch');
                drawnow;

                if no_context,
                    backbone = [coor(site.xax,2:4);coor(site.orig,2:4);coor(site.ypax,2:4)];
                    xyz = get_non_interacting_label_position(backbone,rot_lib);
                    all_pops.pop_rot=1; % net population: ext*int for each rotamer
                    all_pops.ext_net=1; % net Boltzmann factor for external potential (partition sum)
                    all_pops.ext_pop=1; % populations based on the ext potentials only
                    all_pops.partition_function=1; % full partition function
                    rotamers_stats.NOall=[xyz 1]; % NO-centers and weights for all rotamers in a single matrix;
                    rotamers_stats.all_rotamers=[]; % contains coordinates for all rotated rotamers for current position
                    rotamers_stats.all_potentials=all_pops; % total potential for the mutation position;
                    rotamers_stats.labels_own=xyz; % rotamer coordinates in the local residue frame (Ca is always [0,0,0]);
                    rotamers_stats.loc_frame_Ca=site.orig; % global (protein) coordinates of the local frame origin (Ca alpha of the mutated residue for peptides);
                    rotamers_stats.ext_poten_type = 'no_context'; % type of the interatomic potential used
                    rotamers_stats.rot_clash=[];  % save clash statistics if needed
                    msg.error=0;
                    msg.text='Context-free label computation successful.';
                else
                    if hMain.dynamic_rotamers,
                        add_msg_board('Computing static backbone interactions');
                        [backbone_stats,msg]=get_rotamers(bcoor,site,calc_opt,rot_lib,Bfac); % this is required since SCWRL4 ignores backbone interactions of frame atoms
                        add_msg_board(msg.text);
                        add_msg_board('Computing dynamic sidechain interactions');
                        [rotamers_stats,msg]=dynamic_rotamer_populations(resadr,backbone_stats,calc_opt.T,rot_lib.label,rot_lib.MD2PDB,rot_lib.calibration,packing_radius);
                    else
                        [rotamers_stats,msg]=get_rotamers(coor,site,calc_opt,rot_lib,Bfac);
                    end;
                end;

                add_msg_board(msg.text);
                potentials=rotamers_stats.all_potentials;
                if hMain.statistics && ~isempty(potentials) && sum(potentials.pop_rot)>=1e-6,
                    wr_rotamer_statistics(stat_file,resadr,rot_lib.calibration.T,rotamers_stats,threshold,calc_opt.T,calc_opt.forgive,rot_lib.label);
                    if hMain.rotamer_PDB,
                        mydir=pwd;
                        if ~isempty(PDB_path),
                            cd(PDB_path);
                        end;
                        wr_rotamer_pdb(resadr,rotamers_stats,T,rot_lib.label,rot_lib.MD2PDB,library)
                        cd(mydir);
                    end;
                end;
            end;
        else
            resadr=mk_address(cindices);
            add_msg_board(sprintf('Rotamers of %s previously computed in site scan.',resadr));
            potentials=my_rotamers.potentials;
        end;

        if isempty(potentials),
            add_msg_board(msg.text);
            add_msg_board(sprintf('Rotamer potentials are missing. Residue %s not labeled.',mk_address(cindices)));
            set(hMain.figure,'Pointer','arrow');
            drawnow;
        elseif sum(potentials.pop_rot)<1e-6,
            add_msg_board(msg.text);
            add_msg_board(sprintf('### Warning ### All populations are zero. Residue %s not labeled.',mk_address(cindices)));
            set(hMain.figure,'Pointer','arrow');
            drawnow;
        else
            respoi=respoi+1; % one more residue successfully labeled
            if ~computed,
                calc_positions(respoi).label=rot_lib.label;
                calc_positions(respoi).indices=cindices;
                calc_positions(respoi).NOpos=rotamers_stats.NOall;
                calc_positions(respoi).partition_function=rotamers_stats.all_potentials.partition_function;
                calc_positions(respoi).T=T;
                calc_positions(respoi).potentials=potentials;
                calc_positions(respoi).attach_frame = attach_frame;
                if isfield(rot_lib,'exclude'),
                    calc_positions(respoi).exclude = rot_lib.exclude;
                else
                    calc_positions(respoi).exclude = ':';
                end;
                [populations,rotamer_numbers]=sort(potentials.pop_rot,2,'descend');
                cumul=cumsum(populations);
                rotnum=length(cumul(cumul<1-threshold))+1;
                % for debugging
                %     figure(1); clf;
                %     plot(populations,'k.');
                %     hold on;
                %     plot(cumul,'r');
                %     plot([rotnum,rotnum],[0,1],'g:');
                pop=rotamers_stats.NOall(:,4);
                pop=pop/sum(pop);
                xmean=sum(rotamers_stats.NOall(:,1).*pop);
                ymean=sum(rotamers_stats.NOall(:,2).*pop);
                zmean=sum(rotamers_stats.NOall(:,3).*pop);
                dx=(rotamers_stats.NOall(:,1)-xmean);
                dy=(rotamers_stats.NOall(:,2)-ymean);
                dz=(rotamers_stats.NOall(:,3)-zmean);
                nNO=length(dx);
                if nNO > 1,
                    rmsd=sqrt(nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1));
                else
                    rmsd = 0;
                end;
                add_msg_board(sprintf('Standard deviation of label position: %4.1f Å',rmsd));
                add_msg_board(sprintf('%i rotamer(s) cover(s) %3.1f%% of all population',rotnum,100*(1-threshold)));
                add_msg_board(sprintf('%i out of %i objects completed',kres,mres));
                drawnow;
                rotamer_numbers=rotamer_numbers(1:rotnum);
                populations=populations(1:rotnum);
                populations=populations/sum(populations); % renormalization of populations

                for rotamer=1:length(rotamer_numbers),
                    ecoor=rot_lib.library(rotamer_numbers(rotamer)).ecoor;
                    [mm,nn]=size(ecoor);
                    coor=zeros(length(rot_lib.MD2PDB),3);
                    for kk=1:length(rot_lib.MD2PDB),
                        if rot_lib.MD2PDB(kk)>0 && rot_lib.MD2PDB(kk)<=mm,
                            coor(kk,:)=ecoor(rot_lib.MD2PDB(kk),2:4);
                        else
                            coor(kk,:)=[NaN,NaN,NaN];
                        end;
                    end;
                    calc_positions(respoi).rotamers(rotamer).coor=coor;
                    calc_positions(respoi).rotamers(rotamer).pop=populations(rotamer);
                end;
            else
                calc_positions(respoi).label=my_rotamers.label;
                calc_positions(respoi).indices=my_rotamers.indices;
                calc_positions(respoi).NOpos=my_rotamers.NOpos;
                calc_positions(respoi).partition_function=my_rotamers.partition_function;
                calc_positions(respoi).T=my_rotamers.T;
                calc_positions(respoi).potentials=my_rotamers.potentials;
                calc_positions(respoi).rotamers=my_rotamers.rotamers;
                calc_positions(respoi).attach_frame = attach_frame;
                if isfield(rot_lib,'exclude'),
                    calc_positions(respoi).exclude = rot_lib.exclude;
                else
                    calc_positions(respoi).exclude = ':';
                end;
            end;
        end;
        
        
        set(hMain.figure,'Pointer','arrow');

    end;
    if msg.error==5,
        break;
        add_msg_board('Rotamer computation cancelled.');
        calc_positions=[];
    end;
end;

drawnow;
if mres>5 && comp_status, status_figure(1); end;

% add the reference, if it does not yet exist
rotamer_ref=true;
id=tag2id('Polyhach:2010_rotamers',third_party.tags,[],'|');
if ~isempty(id),
    if isfield(model,'auto_references'),
        if ~isempty(find(id==model.auto_references, 1)),
            rotamer_ref=false;
        end;
    else
        model.auto_references=[];
    end;
    if rotamer_ref,
        if ~isfield(model,'references'),
            model.references(1)=third_party.references(id);
        elseif isempty(model.references)
            model=rmfield(model,'references');
            model.references(1)=third_party.references(id);
        else
            model.references(end+1)=third_party.references(id);
        end;
        model.auto_references(end+1)=id;
    end;
end;

function [ecoor,site,Bfac]=ecoor_structure(indices,waterflag,sidechain,frame,noclash)
% returns xyz coordinates for a structure and a labeling site description 
% for a residue within this structure
% if there are several models, the mean coordinates are returned
%
% indices    index vector of at least length 4, addresses the site
% waterflag  (optional), true: water atoms are included, false or missing:
%            water atoms are excluded, default: false
% sidechain  (optional), true: sidechains of native amino acids are 
%            included, false: only backbone atoms and cofactors, 
%            default: true; for nucleotides the sidechain (nucleobase) 
%            is always included
% frame      (optional) three atoms defining the attachment frame,  
%            MMM tag string, for instance for peptides
%            ':C:N:CA:'; % y' axis atom ypax, x axis atom xax, origin orig
%            defaults to the peptide case
% noclash    atoms of the attachment residue that do not need to be tested 
%            for clashes, defaults to *, which means not clash test
%
% ecoor      extended coordinate array, column 1: atomic numbers, columns
%            2-4: xyz coordinates of all atoms in the structure (possibly
%            excluding water)
% site       site description, struct with fields
%            site.res_atoms  vector of integer, indices into ecoor for all
%                            atoms of the labeling site
%            three indices for attachment frame atoms, for instance for
%            peptides:
%            site.N          integer, index into ecoor for backbone N of
%                            site
%            site.C          integer, index into ecoor for backbone
%                            carbonyl C of site
%            site.CA         integer, index into ecoor for C^alpha of
%                            site
%            otherwise, dynamic field names created from frame are used
% Bfac       B factors corresponding to extended coordinates
%
% if there are several locations for the site, only the first one is
% considered for the backbone coordinates

global model

if nargin<2,
    waterflag=false;
end;

if nargin<3 || isempty(sidechain),
    sidechain=true;
end;

if ~exist('frame','var') || isempty(frame),
    frame = ':C:N:CA:';
end;

if ~exist('noclash','var'),
    noclash = '*';
end;

ypax_name = id2tag(1,frame);
xax_name = id2tag(2,frame);
orig_name = id2tag(3,frame);

max_atoms=100000;
ecoor=zeros(max_atoms,4);
Bfac=zeros(1,max_atoms);
translation=zeros(1,max_atoms);
site.res_atoms = [];
site.ypax = [];
site.xax = [];
site.orig = [];

% get all atom indices of the residue into the MMM object structure and
% check whether they need to be excluded from clash tests, index_array has
% all atom numbers that must be excluded from clash tests
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
% if info.type ~= 1 && info.type ~= 2,
%     add_msg_board(sprintf('Residue %s is neither an amino acid nor a nucleotide. Not labeled.',mk_address(indices,1)));
%     ecoor=[];
%     return;
% end;
poi=0;
index_array=zeros(5*length(info.atom_numbers),1);
for k=1:length(info.atom_numbers),
    numbers=info.atom_numbers{k};
    atags=info.atom_tags;
    atag=id2tag(k,atags);
    bid=tag2id(atag,':N:CA:C:O:');
    if noclash == '*',
        is_clash = false;
    elseif isempty(noclash),
        is_clash = true;
    else
        id_noclash = tag2id(atag,noclash);
        is_clash = isempty(id_noclash);
    end;
%     if is_clash,
%         fprintf(1,'Atom %s is tested for clashes.\n',atag);
%     else
%         fprintf(1,'Atom %s is not tested for clashes.\n',atag);
%     end;
    if ~is_clash,
        if sidechain || ~isempty(bid) || info.type == 2,
            [mm,nn]=size(numbers);
            index_array(poi+1:poi+mm)=numbers(:,1);
            poi=poi+mm;
        end;
    end;
end;
index_array=index_array(1:poi);

% check if the attachment frame exists and store indices
id=tag2id(xax_name,info.atom_tags);
if isempty(id),
    add_msg_board(sprintf('Attachment frame atom %s of residue %s is missing. Not labeled.',xax_name,mk_address(indices,1)));
    ecoor=[];
    return;
end;
numbers=info.atom_numbers{id};
xax = numbers(1,1);

id=tag2id(ypax_name,info.atom_tags);
if isempty(id),
    add_msg_board(sprintf('Attachment frame atom %s of residue %s is missing. Not labeled.',ypax_name,mk_address(indices,1)));
    ecoor=[];
    return;
end;
numbers=info.atom_numbers{id};
ypax = numbers(1,1);

id=tag2id(orig_name,info.atom_tags);
if isempty(id),
    add_msg_board(sprintf('Attachment frame atom %s of residue %s is missing. Not labeled.',orig_name,mk_address(indices,1)));
    ecoor=[];
    return;
end;
numbers=info.atom_numbers{id};
orig = numbers(1,1);

poi=0;

nc=length(model.structures{indices(1)}); % number of chains
for k=1:nc,
    nmod=length(model.structures{indices(1)}(k).xyz);
    if indices(3)<=nmod,
        chindices=[indices(1),k,indices(3)];
    else
        chindices=[indices(1),k,1];
    end;
    info=model.structures{chindices(1)}(chindices(2)).residues{chindices(3)}.info;
    elements=model.structures{chindices(1)}(chindices(2)).isotopes(:,1);
    xyz_ch=model.structures{chindices(1)}(chindices(2)).xyz{chindices(3)};
    Bfac_ch=model.structures{chindices(1)}(chindices(2)).Bfactor{chindices(3)};
    for kk=1:length(info),
        if waterflag || ~strcmpi(info(kk).name,'HOH'),
            isaminoacid=info(kk).type==1;
            ishet=info(kk).hetflag;
            for kkk=1:length(info(kk).atom_numbers), % loop over all atoms in one rotamer
                anum=info(kk).atom_numbers{kkk};
                atags=info(kk).atom_tags;
                atag=id2tag(kkk,atags);
                bid=tag2id(atag,':N:CA:C:O:');
                elm = info(kk).elements(kkk);
                if sidechain || ~isempty(bid) || ~isaminoacid || ishet,
                    [ma,na]=size(anum);
                    for k4=1:ma, % loop over all rotamers
                        poi=poi+1;
                        % ecoor(poi,1)=elements(anum(k4,1));
                        ecoor(poi,1)= elm;
                        if ecoor(poi,1)==0,
                            disp('Aber hallo');
                        end;
                        ecoor(poi,2:4)=xyz_ch(anum(k4,1),:);
                        Bfac(poi)=Bfac_ch(anum(k4,1));
                        if k==indices(2),
                            translation(poi)=anum(k4,1);
                        end;
                    end;
                end;
            end;
        end;
    end;
end;
ecoor=ecoor(1:poi,:);
translation=translation(1:poi);
Bfac=Bfac(1:poi);

% translate atom numbers of labeling site
site.res_atoms=index_array;
[m,n]=size(index_array);
for k=1:m,
    poi=find(translation==index_array(k));
    if isempty(poi) || numel(poi)~=1,
        ecoor=[];
        site.res_atoms=[];
        Bfac=[];
        add_msg_board(sprintf('Structure corrupted. Atoms of residue %s missing or ambiguous. Residue is not labelled.',mk_address(indices,1)));
        return;
    end;
    site.res_atoms(k)=poi;
end;
poi=find(translation == xax);
if isempty(poi) || numel(poi)~=1,
    ecoor=[];
    site.res_atoms=[];
    Bfac=[];
    add_msg_board(sprintf('Structure corrupted. Attachment frame atom %s of residue %s missing or ambiguous. Residue is not labelled.',xax_name,mk_address(indices,1)));
    return;
end;
site.xax = poi;
poi=find(translation == ypax);
if isempty(poi) || numel(poi)~=1,
    ecoor=[];
    site.res_atoms=[];
    Bfac=[];
    add_msg_board(sprintf('Structure corrupted. Attachment frame atom %s of residue %s missing or ambiguous. Residue is not labelled.',ypax_name,mk_address(indices,1)));
    return;
end;
site.ypax = poi;
poi=find(translation == orig);
if isempty(poi) || numel(poi)~=1,
    ecoor=[];
    site.res_atoms=[];
    Bfac=[];
    add_msg_board(sprintf('Structure corrupted. Attachment frame atom %s of residue %s missing or ambiguous. Residue is not labelled.',orig_name,mk_address(indices,1)));
    return;
end;
site.orig = poi;
