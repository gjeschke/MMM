function success = make_loop_and_save(restraints,restrain,aux)

global model
global general
global Ramachandran
global residue_defs

parnum = 100; % number of trials performed in the parfor loop

disp_update = 200;

auxiliary = false;

secvec = zeros(1,length(restrain));
for k = 1:length(restrain)
    secvec(k) = restrain(k).secondary;
end

if ~isempty(restraints.Na_indices)
    Nanchor_res = restraints.Na_indices(4);
else
    Nanchor_res = [];
end

if ~isempty(restraints.Ca_indices)
    Canchor_res = restraints.Ca_indices(4);
else
    Canchor_res = [];
end

% Make auxiliary restraints by triangle bound smoothing
% find all beacon residues
beacons = zeros(1,1000);
bref = zeros(1000,2);
bpoi = 0;
for k = 1:length(restrain),
    for kr = 1: length(restrain(k).r_beacon),
        res = restrain(k).r_beacon(kr).resb;
        if isempty(find(beacons == res)),
            bpoi = bpoi+1;
            beacons(bpoi) = res;
            bref(bpoi,1) = k;
            bref(bpoi,2) = kr;
        end
    end
end
ni = length(restrain);
vert = ni + bpoi;
beacons = beacons(1:bpoi);
bref = bref(1:bpoi,:);
lower_bounds = zeros(vert);
upper_bounds = 1e6*ones(vert);
missing = ones(vert);
for k = 1:length(restrain)
    for kr = 1: length(restrain(k).r_beacon)
        if strcmp(restrain(k).r_beacon(kr).type,'Gaussian')
            b = find(beacons == restrain(k).r_beacon(kr).resb) + ni;
            rmean = restrain(k).r_beacon(kr).par1;
            sigr = restrain(k).r_beacon(kr).par2;
            lower_bounds(k,b) = rmean - 2*sigr;
            lower_bounds(b,k) = rmean - 2*sigr;
            upper_bounds(k,b) = rmean + 2*sigr;
            upper_bounds(b,k) = rmean + 2*sigr;
            missing(k,b) = 0;
            missing(b,k) = 0;
        end
    end
    for kr = 1: length(restrain(k).r_intern)
        if strcmp(restrain(k).r_intern(kr).type,'Gaussian')
            b = restrain(k).r_intern(kr).site;
            rmean = restrain(k).r_intern(kr).par1;
            sigr = restrain(k).r_intern(kr).par2;
            lower_bounds(k,b) = rmean - 2*sigr;
            lower_bounds(b,k) = rmean - 2*sigr;
            upper_bounds(k,b) = rmean + 2*sigr;
            upper_bounds(b,k) = rmean + 2*sigr;
            missing(k,b) = 0;
            missing(b,k) = 0;
        end
    end
end
[~,~,err] = triangle_bound_smoothing(lower_bounds,upper_bounds);
if err
    return;
end

if auxiliary
    for k1 = 1:vert-1,
        for k2 = k1+1:vert,
            setbound = false;
            if abs(lb(k1,k2) - lower_bounds(k1,k2)) > eps,
                if missing(k1,k2),
                    if lb(k1,k2) > 2,
                        setbound = true;
                    else
                        setbound = false;
                    end
                end
            end
            if abs(ub(k1,k2) - upper_bounds(k1,k2)) > eps,
                if missing(k1,k2),
                    if ub(k1,k2) < 150,
                        setbound = true;
                    end
                end
            end
            if setbound,
                if k2 > ni,
                    fprintf(1,'Setting auxiliary beacon restraint:');
                    kr = length(restrain(k1).r_beacon)+1;
                    if kr > 1,
                        label_type =  restrain(k).r_beacon(1).label_type;
                        label_T =  restrain(k).r_beacon(1).label_T;
                    else
                        label_type =  restrain(k).r_intern(1).label_type;
                        label_T =  restrain(k).r_intern(1).label_T;
                    end
                    b = k2-ni;
                    kb = bref(b,1);
                    krb = bref(b,2);
                    restrain(k1).r_beacon(kr).xyz = restrain(kb).r_beacon(krb).xyz;
                    restrain(k).r_beacon(kr).label_type = label_type;
                    restrain(k).r_beacon(kr).label_T = label_T;
                    restrain(k).r_beacon(kr).bindices = restrain(kb).r_beacon(krb).bindices;
                    restrain(k).r_beacon(kr).resb = restrain(kb).r_beacon(krb).resb;
                    restrain(k).r_beacon(kr).type = 'bounds';
                    restrain(k).r_beacon(kr).par1 = -lb(k1,k2);
                    restrain(k).r_beacon(kr).par2 = -ub(k1,k2);
                    resa = restraints.res1 + k1 - 1;
                    fprintf(1,'(%i,%i): [%4.2f,%4.2f] Å\n',resa,restrain(k).r_beacon(kr).resb,lb(k1,k2),ub(k1,k2));
                else
                    fprintf(1,'Setting auxiliary internal restraint:');
                    kr = length(restrain(k1).r_intern)+1;
                    if kr > 1,
                        label_type =  restrain(k).r_intern(1).label_type;
                        label_T =  restrain(k).r_intern(1).label_T;
                    else
                        label_type =  restrain(k).r_beacon(1).label_type;
                        label_T =  restrain(k).r_beacon(1).label_T;
                    end
                    restrain(k).r_intern(kr).site = k2;
                    restrain(k).r_intern(kr).label_type = label_type;
                    restrain(k).r_intern(kr).label_T = label_T;
                    restrain(k).r_intern(kr).resb = restraints.res1 + k2 - 1;
                    restrain(k).r_intern(kr).type = 'bounds';
                    restrain(k).r_intern(kr).par1 = -lb(k1,k2);
                    restrain(k).r_intern(kr).par2 = -ub(k1,k2);
                    resa = restraints.res1 + k1 - 1;
                    fprintf(1,'(%i,%i): [%4.2f,%4.2f] Å\n',resa,restrain(k).r_beacon(kr).resb,lb(k1,k2),ub(k1,k2));
                end
            end
        end
    end
end

min_approach = 1.2;

directory = pwd;

display_terminus = false;
attach_it = true;

rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs

pmodel = aux.p_model;
pthr = exp(-erfinv(pmodel)^2);

max_models = restraints.ensemble;

all_p_model = zeros(1,max_models);

snum_vec = zeros(1,max_models);

ntrials = 20000000; % number of Monte Carlo trials
max_seconds = 3600*restraints.max_time; % maximum runtime in seconds


load([general.Ramachandran 'Ramachandran_disordered']);

closed_loop = false;
free_standing = false;

sequence = restraints.sequence;
if isempty(restraints.anchorN) && ~isempty(restraints.anchorC) % reverse model
    reverse = true;
    template_indices = restraints.Ca_indices(1:3);
    sequence = [sequence restraints.Cseq(1)];
elseif ~isempty(restraints.anchorN) % loop anchored at N terminus
    reverse = false;
    template_indices = restraints.Na_indices(1:3);
    sequence = [restraints.Nseq(2) sequence];
    if ~isempty(restraints.anchorC) % closed loop
        sequence = [sequence restraints.Cseq];
        closed_loop = true;
        if attach_it
            if sum(abs(restraints.Na_indices(1:3)-restraints.Ca_indices(1:3))) > 0
                add_msg_board('Warning: N-terminal and C-terminal anchors are in different chains or chain models.');
                add_msg_board('Loop models cannot be attached.');
                attach_it = false;
            end
        end
    end
else % free-standing loop
    reverse = false;
    closed_loop = false;
    free_standing = true;
    template_indices = [];
    if attach_it
        add_msg_board('Warning: No anchor residues defined. Loop is not attached, but models are collected into a single structure.');
    end
%     attach_it = false;
end

if attach_it && ~free_standing % make template copies in current structure
    for modnum = 2:max_models
        copy_structure(template_indices(1),'+mod',[],modnum,template_indices(1));
    end
end


res1 = restraints.res1;
resend = restraints.rese;
terminus = resend;
N_terminus = res1;
if reverse
    terminus = res1;
end

if ~isempty(restraints.anchorN)
    display_start = res1-1;
else
    display_start = res1;
end

if ~isempty(restraints.anchorC)
    display_end = resend+1;
else
    display_end = resend;
end

if ~free_standing
    adr=mk_address(restraints.template);
    aux.fname_bas = strcat(adr,'_loop');
else
    aux.fname_bas = 'free_loop';
end

fname = sprintf('%s_%i_%i',aux.fname_bas,res1,resend);

min_prob = pthr^aux.n_restraints;

add_msg_board(sprintf('Cumulative probability threshold for %i restraints is %6.4f.',aux.n_restraints,min_prob));

if ~free_standing
    [msg,prot_coor] = get_structure(restraints.template,'xyz');
    if msg.error
        add_msg_board('ERROR: Coordinates of template structure could not be retrieved.');
    end
else
    prot_coor = [];
end

success = 0;
err_count=zeros(1,11);
Ram_fixed = 0;
Ram_fix_clash = 0;
resax = res1:resend;
res_stat = zeros(1,length(resax));
set(gcf,'Pointer','watch');
distributions = cell(aux.n_restraints+aux.b_restraints);
restraint_distr = cell(aux.n_restraints);
descriptors = cell(aux.n_restraints);
monitor_distr = cell(aux.n_monitor);
monitor_descr = cell(aux.n_monitor);
rax = get_distribution;
for k = 1:aux.n_restraints + aux.b_restraints
    distributions{k} = zeros(1,length(rax));
    restraint_distr{k} = zeros(1,length(rax));
end
for k = 1:aux.n_monitor
    monitor_distr{k} = zeros(1,length(rax));
end
for k = 1:aux.n_restraints + aux.n_monitor
    aux.exp_distr{k} = [];
    aux.exp_r{k} = [];
end

fullname = fullfile(aux.save_path,sprintf('%s_%s.log',aux.save_name,datestr(now,30)));
fid_report = fopen(fullname,'wt');

mname_stub = fullfile(aux.save_path,aux.save_name);
add_msg_board(sprintf('Saving models to %s_m#.pdb\n',mname_stub));

p_coor = cell(parnum);
p_restrain = cell(parnum);
p_errcode = zeros(1,parnum);
p_cumprob = zeros(1,parnum);
p_kres = zeros(1,parnum);
drawnow;
tic;
kMC = 1;
p_anchorC = restraints.anchorC;
p_anchorCn = restraints.anchorCn;
p_anchorN = restraints.anchorN;
p_anchorNp = restraints.anchorNp;
n_restraints = aux.n_restraints;

rescodes = zeros(1,length(sequence));
for k = 1:length(sequence)
    rescodes(k) = strfind(residue_defs.single_letter_code,sequence(k));
end

Rama_res.me = Ramachandran.me;
Rama_res.ephi = Ramachandran.ephi;
Rama_res.epsi = Ramachandran.epsi;
Rama_res.allowed_P = Ramachandran.allowed_P;
Rama_res.allowed_G = Ramachandran.allowed_G;
Rama_res.allowed_gen = Ramachandran.allowed_gen;

while kMC <= ntrials
    parfor kp = 1:parnum % parfor
        if reverse
            [coor,errcode,restrain1,cumprob,kres] = mk_loop_model_reverse(sequence, p_anchorC, p_anchorCn, prot_coor, restrain, Rama_res, rescodes, min_prob,n_restraints);
        else
            [coor,errcode,restrain1,cumprob,kres] = mk_loop_model(sequence, p_anchorN, p_anchorC, p_anchorNp, p_anchorCn, prot_coor, restrain, Rama_res, rescodes, min_prob, n_restraints);
            kres = kres-1;
        end
        p_coor{kp} = coor; 
        p_errcode(kp) = errcode;
        p_restrain{kp} = restrain1;
        p_cumprob(kp) = cumprob;
        p_kres(kp) = kres;
    end
    for kp = 1:parnum
        kMC = kMC + 1;
        coor = p_coor{kp};
        errcode = p_errcode(kp);
        restrain1 = p_restrain{kp};
        cumprob = p_cumprob(kp);
        kres = p_kres(kp);
%         for k = 1:length(restrain),
%             restrain(k).secondary = secvec(k);
%         end
        res_stat(kres) = res_stat(kres)+1;
        runtime = toc;
        if errcode == -1
            Ram_fixed = Ram_fixed + 1;
            errcode = 0;
        end
        if errcode == -4
            Ram_fixed = Ram_fixed + 1;
            Ram_fix_clash = Ram_fix_clash + 1;
            errcode = 4;
        end
        if aux.n_restraints == 0
            p_model = 1;
        else
            p_model = erf(sqrt(-log(cumprob)/aux.n_restraints));
        end
        err_count(errcode+1) = err_count(errcode+1) + 1;

        if ~errcode
            tpm = runtime/err_count(1);
            fprintf(1,'Time per model: %8.1f s\n',tpm);
            success = success + 1;
            if success == 1
                bb0 = coor;
            elseif  isempty(restraints.anchorN) && isempty(restraints.anchorC)
                [rms,coor] = rmsd_superimpose(bb0,coor);
                add_msg_board(sprintf('Model superimposes onto first model with rmsd of %4.1f Å',rms));
            end
            loopname = write_pdb_backbone(coor,restraints.sequence,fname,success,res1,directory);
            pmodel = make_SCWRL4_sidegroups(loopname,directory);
            [pclash,iclash] = check_decorated_loop(pmodel,prot_coor,res1,resend,min_approach,directory);

            if pclash
    %             fprintf(1,'Loop with sidegroups clashes with protein\n');
                err_count(10) = err_count(10) + 1;
                success = success - 1;
            elseif iclash
    %             fprintf(1,'Loop with sidegroups clashes with itself\n');
                err_count(11) = err_count(11) + 1;
                success = success - 1; 
            else
                all_p_model(success) = p_model;
                if attach_it
                    if free_standing && success == 1
                        [msg,snum]=add_pdb(pmodel);
                        template_indices = [snum 1 1];
                        model_indices = template_indices;
                        for modnum = 2:max_models
                            copy_structure(template_indices(1),'+mod',[],modnum,template_indices(1));
                        end
                    else
                        structure=rd_pdb(pmodel);
                        model_indices = template_indices;
                        model_indices(3) = success;
                        if free_standing
                            replace_model(model_indices,structure);
                        else
                            insert_residues(res1,resend,structure,model_indices,Nanchor_res,Canchor_res);
                        end
                        if isfield(model,'selected')
                            model = rmfield(model,'selected');
                        end
                        model.selected{1} = model_indices;
                        wr_pdb_selected(sprintf('%s_m%i',aux.save_name,success),aux.pdbid);
                    end
                end
            end
        end
        if mod(kMC,disp_update) == 0
            ftr = (1 - kMC/ntrials)*max_seconds;
            fti = max_seconds - runtime;
            fmo = runtime*(max_models-success)/success;
            time_left = min([fti fmo ftr]);
            hours = floor(time_left/3600);
            minutes = round((time_left-3600*floor(time_left/3600))/60);
            if minutes == 60
                hours = hours + 1;
                minutes = 0;
            end
            fprintf(1,'%i h %i min estimated run time to completion\n',hours,minutes); 
        end
        if success >= max_models
            break
        end
    end
    if success >= max_models
        break
    end
%     if mod(kMC,10000)
%         fprintf(1,'%6.3f ms per trial\n',1000*runtime/kMC);
%     end
    if runtime >= max_seconds
        add_msg_board('Warning: Ensemble generation stopped as maximum allotted time was exceeded.');
        break
    end
end

runtime = toc;
hours = floor(runtime/3600);
minutes = round((runtime-3600*floor(runtime/3600))/60); 

fprintf(fid_report,'\n--- Final statistics ---\n\n');
fprintf(fid_report,'Total runtime: %i h %i min.\n',hours,minutes);
fprintf(fid_report,'%5.2f%% of maximum number of MC trials spent.\n',100*kMC/ntrials);
fprintf(fid_report,'Ensemble with %i models generated.\n',success);
fprintf(fid_report,'%i backbone models were generated.\n',err_count(1));
fprintf(fid_report,'%5.2f%% of all trials had restraint violations.\n',100*err_count(6)/kMC);
fprintf(fid_report,'This is a success rate of %6.1f ppm.\n',1e6*(1-err_count(6)/kMC));
fprintf(fid_report,'%5.2f%% of all trials had internal loop backbone clashes.\n',100*(err_count(3)+err_count(8))/kMC);
fprintf(fid_report,'%5.2f%% of all trials hand backbone clashes with the protein.\n',100*(err_count(5)+err_count(7))/kMC);
fprintf(fid_report,'The sidegroup clash threshold was %4.1f Å.\n',min_approach);
fprintf(fid_report,'%5.2f%% of the backbone models lead to internal sidegroup clashes.\n',100*err_count(11)/err_count(1));
fprintf(fid_report,'%5.2f%% of the backbone models lead to sidegroup clashes with the protein.\n',100*err_count(10)/err_count(1));
if closed_loop
    fprintf(1,'%5.2f%% of all trials lead to loops that could not be closed.\n',100*err_count(2)/kMC);
    fprintf(1,'%5.2f%% of all trials lead to loops for which the link Ramachandran angles could not be fixed.\n',100*err_count(4)/kMC);
end

hours = floor(max_seconds/3600);
minutes = round((max_seconds-3600*floor(max_seconds/3600))/60);

fprintf(fid_report,'\nThe run had a time limit of: %i h %i min.\n',hours,minutes);
fprintf(fid_report,'The run had a limit of %i Monte Carlo trials.\n',ntrials);
fprintf(fid_report,'%i models had been requested with an ensemble probability of %4.2f.\n',max_models,aux.p_model);
fprintf(fid_report,'The number of Gaussian-type restraints was %i.\n',aux.n_restraints);

fclose(fid_report);




if kMC == ntrials
    add_msg_board('Warning: Ensemble generation stopped as maximum allotted number of Monte Carlo trials was exceeded.');
end
if attach_it && success > 0
    consolidate_chain(template_indices(1:2));
end

function loopname = write_pdb_backbone(coor,sequence,fname0,model,res1,directory)

my_dir = pwd;
cd(directory);

residues='ALAARGASNASPCYSGLNGLUGLYHISILELEULYSMETPHEPROSERTHRTRPTYRVAL';
oneletter='ARNDCQEGHILKMFPSTWYV';

fname = sprintf('%s_m%i',fname0,model);

loopname = [fname '.pdb'];
wfile=fopen(loopname,'w');
for k = 1:length(sequence)
    respoi=strfind(oneletter,sequence(k));
    residue=residues(1+3*(respoi-1):3+3*(respoi-1));
    N = coor(4*k-3,:);
    CA = coor(4*k-2,:);
    C = coor(4*k-1,:);
    O = coor(4*k,:);
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           N\n','ATOM  ',4*k-3,'  N   ',residue,k+res1-1,N(1),N(2),N(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM  ',4*k-2,'  CA  ',residue,k+res1-1,CA(1),CA(2),CA(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM  ',4*k-1,'  C   ',residue,k+res1-1,C(1),C(2),C(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           O\n','ATOM  ',4*k,'  O   ',residue,k+res1-1,O(1),O(2),O(3));
end
fclose(wfile);

cd(my_dir);

function [outname,status,result] = make_SCWRL4_sidegroups(inname,directory)
%
% Attaches or corrects sidegroups using SCWRL4
% the program SCWRL4 needs to be on the current Matlab path
%
% inname   name of the PDB file to which sidegroups should be attached
% outname  name of the output PDB file with SCWRL4 sidegroups
%

my_dir = pwd;
cd(directory);

poi = strfind(inname,'.pdb');
outname = [inname(1:poi-1) '_SCWRL4.pdb'];

s=which('scwrl4.exe');
cmd=[s ' -i ' inname ' -o ' outname];
[status,result]=dos(cmd);

cd(my_dir);

function [pclash,iclash,approach_prot,approach_loop] = check_decorated_loop(loopname,prot_coor,res1,resend,min_approach,directory)

my_dir = pwd;
cd(directory);

if ~exist('min_approach','var')
    min_approach = 1.2; 
end

approach_prot = -1;
approach_loop = -1;
pclash = 1;
iclash = 1;
loop_coor = zeros(5000,3);
l_res_assign = zeros(5000,3);
fid=fopen(loopname);
if fid==-1
    add_msg_board(sprintf('Warning: Loop structure PDB file %s missing. Rejected.',loopname));
    cd(my_dir);
    return;
end
poi = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if length(tline) >= 6
        record=tline(1:6);
        resnum = str2double(tline(23:26));
        if strcmpi(record,'ATOM  ') || strcmpi(record,'HETATM')
            if length(tline) < 78 || tline(78)~='H'
                if ~strcmpi(tline(14),'H') && resnum ~= res1 && resnum ~= resend
                    poi = poi +1;
                    l_res_assign(poi) = resnum;
                    valstr = [tline(31:38) ' ' tline(39:46) ' ' tline(47:54)];
                    loop_coor(poi,:) = str2num(valstr);
                end
            end
        end
    end
end
fclose(fid);

cd(my_dir);

loop_coor = loop_coor(1:poi,:);
l_res_assign = l_res_assign(1:poi);

pclash = 0;
iclash = 0;

[m1,~] = size(loop_coor); % get sizes of the coordinates arrays
[m2,~] = size(prot_coor);

if m2 > 0
    a2 = repmat(sum(loop_coor.^2,2),1,m2);
    b2 = repmat(sum(prot_coor.^2,2),1,m1).';
    pair_dist = sqrt(abs(a2 + b2 - 2*loop_coor*prot_coor.'));
    min_dist = min(min(pair_dist));

    approach_prot = min_dist;
    if min_dist < min_approach
    %    fprintf(2,'Minimum sidegroup distance to protein is %6.2f Å\n',min_dist);
       pclash = 1;
       cd(my_dir);
       return
    end
end

min_dist = 1e6;
% test for minimum distance within loop
we_clash = zeros(1,2);
for k1 = 1:poi-1
    for k2 = k1+1:poi
        if abs(l_res_assign(k1)-l_res_assign(k2))>1
            approach = norm(loop_coor(k1,:) - loop_coor(k2,:));
            if approach < min_dist
                min_dist = approach;
                we_clash = [k1,k2];
            end
        end
    end
end
approach_loop = min_dist;
if min_dist < min_approach
%    fprintf(2,'Minimum heavy-atom distance: %6.2f Å at (%i,%i)\n',min_dist,we_clash);
    iclash = 1;
end

function insert_residues(res1,rese,structure,template_indices,Nanchor_res,Canchor_res)

% Bfactor is set to twice the maximum Bfactor of existing residues

global model
global residue_defs


snum = template_indices(1);
cnum = template_indices(2);
mnum = template_indices(3);

if ~isempty(Nanchor_res),
    model.structures{snum}(cnum).residues{mnum}.info(Nanchor_res).terminal = 0;
end

Bmax=max(model.info{snum}.B_range);

newresidues=res1:rese;
[ma,na]=size(model.structures{snum}(cnum).xyz{mnum});
[ma0,maxconn]=size(model.structures{snum}(cnum).conn);
newiso=[model.structures{snum}(cnum).isotopes; zeros(20*length(newresidues),2,'single')];
newxyz=[model.structures{snum}(cnum).xyz{mnum}; zeros(20*length(newresidues),3)];
newBfac=[model.structures{snum}(cnum).Bfactor{mnum}, zeros(1,20*length(newresidues))];
newBtens=[model.structures{snum}(cnum).Btensor{mnum}; zeros(20*length(newresidues),6,'int32')];
sequence=model.structures{snum}(cnum).sequence;
restags=model.structures{snum}(cnum).residues{mnum}.residue_tags;

newrnum=zeros(1,length(newresidues));
nr=length(model.structures{snum}(cnum).residues{mnum}.info);
rpoi=0;
for k=1:length(newresidues),
    rnum0=newresidues(k);
    rnum=0;
    for kk=1:length(structure(1).residues{1}.info),
        if structure(1).residues{1}.info(kk).number == rnum0,
            rnum=kk;
            break;
        end
    end
    if rnum==0,
        add_msg_board(sprintf('Warning: Residue %i not found in modelled loop',rnum0));
        continue;
    end
    tag=structure(1).residues{1}.info(rnum).name;
    id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
    sequence(rnum0)=id;
    if mnum==1 && isfield(model.structures{snum}(cnum),'seqexist'),
        model.structures{snum}(cnum).seqexist(rnum0)=1;
    end
    nr=nr+1;
    rpoi=rpoi+1;
    newrnum(rpoi)=nr;
    % model.structures{snum}(cnum).residues{mnum}.info(nr)=structure(1).residues{1}.info(rnum);
    restag=sprintf('%i:',rnum0);
    restags=strcat(restags,restag);
    pointers=structure(1).residues{1}.info(rnum).atom_numbers;
    if isempty(Canchor_res) && k == length(newresidues),
        model.structures{snum}(cnum).residues{mnum}.info(nr).terminal = 1;
    end
    for anum=1:length(pointers), % loop over atoms
        pointer=structure(1).residues{1}.info(rnum).atom_numbers{anum};
        [loc,n]=size(pointer);
        for lnum=1:loc, % loop over locations
            poi=pointer(lnum,1); % actual coordinate set number
            ma=ma+1;
            pointer(lnum,1)=ma;
            newiso(ma,:)=structure(1).isotopes(poi,:);
            newxyz(ma,:)=structure(1).xyz{1}(poi,:);
            newBfac(ma)=2*Bmax;
            newBtens(ma,:)=structure(1).Btensor{1}(poi,:);
        end
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_numbers{anum}=pointer;
        model.structures{snum}(cnum).residues{mnum}.info(nr).name=tag;
        model.structures{snum}(cnum).residues{mnum}.info(nr).type=structure(1).residues{1}.info(rnum).type;
        model.structures{snum}(cnum).residues{mnum}.info(nr).secondary=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).hetflag=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).connected=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).number=structure(1).residues{1}.info(rnum).number;
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_tags=structure(1).residues{1}.info(rnum).atom_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).elements=structure(1).residues{1}.info(rnum).elements;
        model.structures{snum}(cnum).residues{mnum}.info(nr).location_tags=structure(1).residues{1}.info(rnum).location_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).insertion_code=structure(1).residues{1}.info(rnum).insertion_code;
    end
end
newiso=newiso(1:ma,:);
newxyz=newxyz(1:ma,:);
newBfac=newBfac(1:ma);
newBtens=newBtens(1:ma,:);

if mnum==1,
    model.structures{snum}(cnum).sequence=sequence;
    model.structures{snum}(cnum).isotopes=newiso;
    model.structures{snum}(cnum).conn=[model.structures{snum}(cnum).conn; zeros(ma-ma0,maxconn)];
end

model.structures{snum}(cnum).xyz{mnum}=newxyz;
model.structures{snum}(cnum).Bfactor{mnum}=newBfac;
model.structures{snum}(cnum).Btensor{mnum}=newBtens;
model.structures{snum}(cnum).residues{mnum}.residue_tags=restags;

% make internal bonds in new residues
for k=1:length(newrnum),
    if newrnum(k) > 0,
        model.structures{snum}(cnum)=mk_internal_bonds(model.structures{snum}(cnum),newrnum(k),residue_defs);
    else
        disp('Aber Hallo!');
    end
end

% sort residues by number
info=model.structures{snum}(cnum).residues{mnum}.info;
numbers=length(info);
for k=1:length(info),
    numbers(k)=info(k).number;
end
info0=info;
[sorted,oldnumbers]=sort(numbers);
tags=':';
for k=1:length(oldnumbers),
    info(k)=info0(oldnumbers(k));
    tag=id2tag(oldnumbers(k),model.structures{snum}(cnum).residues{mnum}.residue_tags);
    tags=[tags tag ':'];
end
model.structures{snum}(cnum).residues{mnum}.residue_tags=tags;
model.structures{snum}(cnum).residues{mnum}.info=info;

function replace_model(indices,structure)

global model

snum = indices(1);
cnum = indices(2);
modnum = indices(3);

model.structures{snum}(cnum).atoms{modnum} = structure(1).atoms{1};
model.structures{snum}(cnum).residues{modnum} = structure(1).residues{1};
model.structures{snum}(cnum).xyz{modnum} = structure(1).xyz{1};
model.structures{snum}(cnum).Bfactor{modnum} = structure(1).Bfactor{1};
model.structures{snum}(cnum).Btensor{modnum} = structure(1).Btensor{1};
        
