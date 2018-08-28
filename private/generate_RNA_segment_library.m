function generate_RNA_segment_library(RRM,RNA,opt,sel)
% function generate_RNA_segment_library(RRM,RNA,opt,sel)
%
% uses a set of RNA segment PDB files and the structure of an RNA
% recognition motif (RRM) with a short RNA piece bound  
% creates a library of segments fitted to the RRM
% user is prompted for identifying the files
%
% the protein rigid-body model must already be loaded in MMM and must be
% the only selected object
%
% the library is for later use in RigiFlex and contains
%   a paradigm PDB file (for the first entry)
%   coordinate sets for the other entries
%   coordinates of the P atom of the first and last residue for each entry
%
% first step: restraint-based rigid-body docking of the RNA to the protein
% second step: AMBER optimization of the RNA by calling Tinker
%
% the set of RNA segment PDB files can be generated by the Rosetta server
% Rosie FarFar: http://rosie.graylab.jhu.edu/documentation/rna_denovo
%
% all PDB files of the RNA segments must be in the same directory and this
% directory must not contain any other file with extension .pdb
%
% RRM   indices of the RRM on chain level, will be treated as inactive in
%       Tinker minimization of the RNA segment
% RNA   indices of the RNA template describing the binding mode, only
%       backbone structure will be used for initial superposition
% opt   options, structure with fields
%       .align  [nt,2] matrix of nt aligned nucleotides between RNA
%               template and RNA segments
% sel   optional vector of model numbers that should be considered
%
% G. Jeschke, 18.8.2017

global model

sugar_mode = true;

clash_test = opt.clash_test;
replace_mode = opt.replace_mode;

forgive = 1; % scaling factor for van-der-Waals radii in clash tests
clash_threshold = 1.5*forgive; % a uniform van-der-Waals radius of 1.5 � is assumed for heavy atoms

if sugar_mode
    nnt = 10;
else
    nnt = 6;
end;

R = 8.314; % universal gas constant
T = 298; % temperature

tinker_options.active = true;
tinker_options.solvation = 'still';
tinker_options.tolerance = 0.1;
tinker_options.forcefield = 'amber99';
if isfield(opt,'initial') && opt.initial
    tinker_options.initial = true;
end;
% tinker_options.stepmax = 0.02;

[msg,RRM_coor] = get_chain(RRM,'xyz_heavy'); % heavy atom coordinates of the RRM

if msg.error
    add_msg_board(sprintf('Heavy-atom coordinates of the RRM could not be retrieved (%s). Aborting',msg.text));
    return
end

RNA_adr = mk_address(RNA);
[nt,~] = size(opt.align);
temp_coor = zeros(nnt*nt,3);
for k = 1:nt,
    if sugar_mode
        adr1=sprintf('%s%i.P,O5'',C5'',C4'',C3'',O3'',C3'',O3'',C1'',C2'',O2'',O4''',RNA_adr,opt.align(k,1));
    else
        adr1=sprintf('%s%i.P,O5'',C5'',C4'',C3'',O3''',RNA_adr,opt.align(k,1));
    end
    [indices,message]=resolve_address(adr1);
    [m,~]=size(indices);
    if sugar_mode && m ~= 10
        add_msg_board('ERROR: A nucleotide in the RNA template does not exist or misses backbone atoms.');
        add_msg_board(sprintf('ERROR in resolve_address was %s',message.text));
        return
    end;
    if ~sugar_mode && m ~= 6
        add_msg_board('ERROR: A nucleotide in the RNA template does not exist or misses backbone atoms.');
        add_msg_board(sprintf('ERROR in resolve_address was %s',message.text));
        return
    end;
    for ka = 1:m
        [msg,cc0]=get_atom(indices(ka,1:5),'coor');
        if msg.error
            add_msg_board('ERROR: Atom coordinates not found in template.');
            return
        end
        temp_coor(nnt*(k-1)+ka,:) = cc0;
    end
end;

mydir = pwd;
[fname,pname]=uigetfile('*.pdb','Select paradigm RNA segment PDB file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('RNA library generation cancelled by user');
    return
else
    [msg,snum]=add_pdb(fullfile(pname,fname));
    if msg.error
        add_msg_board(sprintf('ERROR in reading paradigm PDB file: %s. Aborting.',msg.text));
        return
    end
    RNA_stag = sprintf('[%i]',snum);
    RNA_snum = snum;
    cd(pname);
    if ~exist('sel','var') || isempty(sel)
        RNA_list = dir('*.pdb');
    else
        for ks = 1:length(sel)
            mtag = sprintf('%i',sel(ks));
            preamb = 'R_000000';
            RNA_list(ks).name = sprintf('%s%s.pdb',preamb(1:8-length(mtag)),mtag);
        end
    end
end
set(gcf,'Pointer','watch');
tic,
paradigm_coor = model.structures{RNA_snum}(1).xyz{1};
segment_adr = mk_address([RNA_snum 1]);
% RNA_list = RNA_list(1:4);
rmsd_vec = zeros(1,length(RNA_list));
rmsd_vec_new = zeros(1,length(RNA_list));
energies = zeros(1,length(RNA_list));
clashes = zeros(1,length(RNA_list));
for k = 1:length(RNA_list) 
    RNA_name = RNA_list(k).name;
    add_msg_board(sprintf('Processing segment %s.',RNA_name));
    structure = rd_pdb(RNA_name);
    RNA_coor = structure(1).xyz{1}; % the segment files should have a single chain and only a single model
    % test whether the current segment conformation has the same number of
    % atoms as the paradigm
    [mp,np] = size(paradigm_coor);
    [mc,nc] = size(RNA_coor);
    if mp ~= mc || np ~= nc
        add_msg_board(sprintf('ERROR: Segment %s does not have the same number of atoms as the paradigm. Aborting.',RNA_name));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    model.structures{RNA_snum}(1).xyz{1} = RNA_coor;
    
    segment_coor = zeros(6*nt,3);
    for knt = 1:nt,
        if sugar_mode
            adr2=sprintf('%s%i.P,O5'',C5'',C4'',C3'',O3'',C3'',O3'',C1'',C2'',O2'',O4''',segment_adr,opt.align(knt,2));
        else
            adr2=sprintf('%s%i.P,O5'',C5'',C4'',C3'',O3''',segment_adr,opt.align(knt,2));
        end
        [indices,message]=resolve_address(adr2);
        [m,~]=size(indices);
        if ~sugar_mode && m ~= 6
            add_msg_board('ERROR: A nucleotide in the segment does not exist or misses backbone atoms.');
            add_msg_board(sprintf('ERROR in resolve_address was %s',message.text));
            return
        end;
        if sugar_mode && m ~= 10
            add_msg_board('ERROR: A nucleotide in the segment does not exist or misses backbone atoms.');
            add_msg_board(sprintf('ERROR in resolve_address was %s',message.text));
            return
        end;
        for ka = 1:m
            [msg,cc0]=get_atom(indices(ka,1:5),'coor');
            if msg.error
                add_msg_board('ERROR: Atom coordinates not found in segment.');
                return
            end
            segment_coor(nnt*(knt-1)+ka,:) = cc0;
        end
    end
    [rmsd,~,transmat]=rmsd_superimpose(temp_coor,segment_coor);
    rmsd_vec(k) = rmsd;
    
    
    add_msg_board(sprintf('Rigid-body superposition of RNA segment onto binding motif template.'));
    RNA_coor_si  = affine_trafo_coor(RNA_coor,transmat);
    model.structures{RNA_snum}(1).xyz{1} = RNA_coor_si;
    
    mdist = min_dist(RRM_coor,RNA_coor_si);
    if mdist > clash_threshold || ~clash_test
    
        if replace_mode
            replace_residues([RNA_snum,1],RNA,fliplr(opt.align)); % fliplr necessary since replacements should be in second column
        end;
        add_msg_board(sprintf('Optimization of the RNA segment by Tinker (amber99).'));
        clear molecule
        clear inactive
        molecule{1} = [RNA_snum,1];
        molecule{2} = RRM;
        active{1} = [RNA_snum,1];
        % active{1} = RRM;

        [msg,snum,energy] = optimize_by_tinker('seg2temp',molecule,active,tinker_options);

        [msg,RNA_coor_opt] = get_chain([snum 2],'xyz_heavy'); % heavy atom coordinates of the RRM

        if msg.error
            add_msg_board(sprintf('Heavy-atom coordinates of the RNA segment could not be retrieved (%s). Aborting',msg.text));
            return
        end

        mdist = min_dist(RRM_coor,RNA_coor_opt);
        if mdist < clash_threshold
            clashes(k) = 1;
        else
            clashes(k) = 0;
        end

        energies(k) = energy;
        segment_adr_new = mk_address([snum 2]);
        segment_coor_new = zeros(nnt*nt,3);
        for knt = 1:nt
            if sugar_mode
                adr2=sprintf('%s%i.P,O5'',C5'',C4'',C3'',O3'',C3'',O3'',C1'',C2'',O2'',O4''',segment_adr_new,opt.align(knt,2));
            else
                adr2=sprintf('%s%i.P,O5'',C5'',C4'',C3'',O3''',segment_adr_new,opt.align(knt,2));
            end
            [indices,message]=resolve_address(adr2);
            [m,~]=size(indices);
            if m ~= nnt
                add_msg_board('ERROR: A nucleotide in the optimized segment does not exist or misses backbone atoms.');
                add_msg_board(sprintf('ERROR in resolve_address was %s',message.text));
                return
            end;
            for ka = 1:m
                [msg,cc0]=get_atom(indices(ka,1:5),'coor');
                if msg.error
                    add_msg_board('ERROR: Atom coordinates not found in optimized segment.');
                    return
                end
                segment_coor_new(nnt*(knt-1)+ka,:) = cc0;
            end
        end;
        rmsd = sqrt(sum(sum((temp_coor-segment_coor_new).^2))/nnt*nt);
        rmsd_vec_new(k) = rmsd;
    else
        clashes(k) = 1;
    end;
end
toc,
set(gcf,'Pointer','arrow');
cd(mydir);

clash_vec = find(clashes == 1);
non_clash_vec = find(clashes == 0);

figure(1); clf;
plot(rmsd_vec,'k.');
set(gca,'FontSize',16);
xlabel('Segment number');
ylabel('Backbone RMSD [�]');
hold on;
plot(non_clash_vec,rmsd_vec_new(non_clash_vec),'ro');
plot(clash_vec,rmsd_vec_new(clash_vec),'r*');

figure(2); clf;
plot(non_clash_vec,energies(non_clash_vec)/1000,'k.');
hold on
plot(clash_vec,energies(clash_vec)/1000,'r*');
set(gca,'FontSize',16);
xlabel('Segment number');
ylabel('Energy [kJ/mol]');

energies(clash_vec) = 1e18;
energies = energies - min(energies);
pops = exp(-energies/(R*T));
pops = pops / sum(pops);

figure(3); clf;
plot(non_clash_vec,pops(non_clash_vec),'k.');
hold on
plot(clash_vec,pops(clash_vec),'r.');
set(gca,'FontSize',16);
xlabel('Segment number');
ylabel('Population');

