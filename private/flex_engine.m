function diagnostics = flex_engine(restraints,restrain,options,handles)
% Make a flexible section for a RigiFlex model
%
% function diagnostics = flex_engine(restraints,options,handles)
%
% restraints    restraint definition for only this flexible section
% restrain      processed restraints that are defined on a per-residue
%               basis
% options       running options
% handles       handles to the RigiFlex GUI, can be missing for modelling
%               on a remote server handles.text_time_left

global general
global Ramachandran
global residue_defs

maxatoms = 50000;
disp_update = 200; % controls, how often progress of computation is displayed
parnum = 100; % parallelization

min_approach = options.min_approach;

diagnostics.success = 0;
diagnostics.runtime = 0;
runtime = 0;

% if handles of the GUI figure are supplied, the engine will update the MMM
% GUI during the computation, else this is considered as a run on a remote
% server that does not have access to the GUI
if exist('handles','var')
    interactive = true;
else
    interactive = false;
end;

secvec = zeros(1,length(restrain));
for k = 1:length(restrain),
    secvec(k) = restrain(k).secondary;
end;

if ~isempty(restraints.Na_indices)
    Nanchor_res = restraints.Na_indices(4);
else
    Nanchor_res = [];
end;

if ~isempty(restraints.Ca_indices)
    Canchor_res = restraints.Ca_indices(4);
else
    Canchor_res = [];
end;

directory = general.tmp_files;

if options.deterministic
    rng(13);
else
    rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs
end;

pmodel = restraints.prob;
pthr = exp(-erfinv(pmodel)^2);

max_models = restraints.models;

all_p_model = zeros(1,max_models);

snum_vec = zeros(1,max_models);

ntrials = options.max_trials; % number of Monte Carlo trials
max_seconds = 3600*options.max_time; % maximum runtime in seconds

time_left = sprintf('%i h %i min estimated run time to completion',floor(max_seconds/3600),round((max_seconds-3600*floor(max_seconds/3600))/60));
set(handles.text_time_left,'String',time_left);
set(handles.text_time_left,'ForegroundColor',[0.75,0,0]);

load([general.Ramachandran 'Ramachandran_disordered']);

closed_loop = false;
free_standing = false;
set(handles.text_SANS_fail,'String','n.a.');
set(handles.text_SAXS_fail,'String','n.a.');

sequence = restraints.sequence;
if isempty(restraints.Na_indices) && ~isempty(restraints.Ca_indices) % reverse model
    reverse = true;
    template_indices = restraints.Ca_indices(1:3);
    sequence = [sequence restraints.Cseq(1)];
elseif ~isempty(restraints.Na_indices) % loop anchored at N terminus
    reverse = false;
    template_indices = restraints.Na_indices(1:3);
    sequence = [restraints.Nseq(2) sequence];
    if ~isempty(restraints.anchorC) % closed loop
        sequence = [sequence restraints.Cseq];
        closed_loop = true;
    end;
else % free-standing loop
    reverse = false;
    closed_loop = false;
    free_standing = true;
    template_indices = [];
end;

res1 = restraints.start;
resend = restraints.end;
terminus = resend;
N_terminus = res1;
if reverse
    terminus = res1;
end;

fname = sprintf('%s_r%i_f%i_%i_%i',options.fname_bas,options.rm,options.fd,res1,resend);

min_prob = pthr^options.n_restraints;

if interactive
    add_msg_board(sprintf('Cumulative probability threshold for %i restraints is %6.4f.',options.n_restraints,min_prob));
end;

if ~free_standing
    [msg,all_chain_coor] = get_object(options.template,'xyz');
    prot_coor = zeros(maxatoms,3);
    poi = 0;
    for kc = 1:length(all_chain_coor)
        xyz = all_chain_coor{kc};
        [m,~] = size(xyz);
        prot_coor(poi+1:poi+m,:) = xyz;
        poi = poi+m;
    end;
    prot_coor = prot_coor(1:poi,:);
    if msg.error
        if interactive
            add_msg_board('ERROR: Coordinates of template structure could not be retrieved.');
            add_msg_board(msg);
        end
        return
    end;
else
    prot_coor = [];
end;

success = 0;
err_count=zeros(1,11);
Ram_fixed = 0;
Ram_fix_clash = 0;
resax = res1:resend;
res_stat = zeros(1,length(resax));
distributions = cell(options.n_restraints);
restraint_distr = cell(options.n_restraints);
descriptors = cell(options.n_restraints);
monitor_distr = cell(options.n_monitor);
monitor_descr = cell(options.n_monitor);
rax = get_distribution;
for k = 1:options.n_restraints
    distributions{k} = zeros(1,length(rax));
    restraint_distr{k} = zeros(1,length(rax));
end;
for k = 1:options.n_monitor
    monitor_distr{k} = zeros(1,length(rax));
end;

% fullname = fullfile(handles.save_path,sprintf('%s_%s.log',handles.save_name,datestr(now,30)));
fullname = fullfile(handles.save_path,handles.report_name);
fid_report = fopen(fullname,'at');

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
n_restraints = options.n_restraints;

rescodes = zeros(1,length(sequence));
for k = 1:length(sequence)
    rescodes(k) = strfind(residue_defs.single_letter_code,sequence(k));
end;

Rama_res.me = Ramachandran.me;
Rama_res.ephi = Ramachandran.ephi;
Rama_res.epsi = Ramachandran.epsi;
Rama_res.allowed_P = Ramachandran.allowed_P;
Rama_res.allowed_G = Ramachandran.allowed_G;
Rama_res.allowed_gen = Ramachandran.allowed_gen;

while kMC <= ntrials
    parfor kp = 1:parnum % ### parfor
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
        res_stat(kres) = res_stat(kres)+1;
        runtime = toc;
        if errcode == -1
            Ram_fixed = Ram_fixed + 1;
            errcode = 0;
        end;
        if errcode == -4
            Ram_fixed = Ram_fixed + 1;
            Ram_fix_clash = Ram_fix_clash + 1;
            errcode = 4;
        end;
        if options.n_restraints == 0
            p_model = 1;
        else
            p_model = erf(sqrt(-log(cumprob)/options.n_restraints));
        end;
        err_count(errcode+1) = err_count(errcode+1) + 1;

        if ~errcode
            tpm = runtime/err_count(1);
            set(handles.text_time_per_model,'String',sprintf('%8.1f',tpm));
            success = success + 1;
            if success == 1
                bb0 = coor;
            elseif  isempty(restraints.anchorN) && isempty(restraints.anchorC)
                [rms,coor] = rmsd_superimpose(bb0,coor);
                if interactive
                    add_msg_board(sprintf('Model superimposes onto first model with rmsd of %4.1f Å',rms));
                end
            end;
            loopname = write_pdb_backbone(coor,restraints.sequence,fname,success,res1,directory);
            [pmodel,status,result] = make_SCWRL4_sidegroups(loopname,directory);
    
            [pclash,iclash] = check_decorated_loop(pmodel,prot_coor,res1,resend,min_approach,directory);

            if pclash
                err_count(10) = err_count(10) + 1;
                success = success - 1;
            elseif iclash
                err_count(11) = err_count(11) + 1;
                success = success - 1; 
            else
                all_p_model(success) = p_model;
                if success == 1
                    [~,snum] = add_pdb(pmodel);
                    template_indices = [snum 1 1];
                    model_indices = template_indices;
                    for modnum = 2:max_models
                        copy_structure(template_indices(1),'+mod',[],modnum,template_indices(1));
                    end;
                else
                    structure = rd_pdb(pmodel);
                    model_indices = template_indices;
                    model_indices(3) = success;
                    replace_model(model_indices,structure);
                end;
                [distributions,restraint_distr,descriptors,monitor_distr,monitor_descr] = mk_report_distributions(fid_report,model_indices,restrain1,options.monitor,p_model,distributions,restraint_distr,descriptors,monitor_distr,monitor_descr,res1);
            end
        end;
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
            end;
            time_left = sprintf('%i h %i min estimated run time to completion',hours,minutes); 
            success_distr = (1-success/kMC)*res_stat/sum(res_stat);
            if reverse
                success_distr = fliplr(success_distr);
            end;
            normalize = 1- cumsum(success_distr);
            normalize(normalize==0)=1;
            for ks = 2:length(success_distr)
                success_distr(ks) = success_distr(ks)/normalize(ks-1);
            end;
            if reverse
                success_distr = fliplr(success_distr);
            end;
            if interactive
                set(handles.text_time_left,'String',time_left);
                axes(handles.axes_multi_plot);
                cla;
                plot(resax,success_distr);
                ylabel('Rejection probability');
                title('Rejection distribution along the loop');
                set(handles.text_percent_time,'String',sprintf('%5.2f',100*runtime/max_seconds));
                set(handles.text_max_trials,'String',sprintf('%5.2f',100*kMC/ntrials));
                set(handles.text_success,'String',sprintf('%i',success));
                set(handles.text_dmg_fail,'String',sprintf('%i',err_count(1)));
                set(handles.text_auxiliary_fail,'String',sprintf('%5.2f',100*err_count(6)/kMC));
                set(handles.text_core_fail,'String',sprintf('%5.2f',100*(err_count(3)+err_count(8))/kMC));
                set(handles.text_linker_fail,'String',sprintf('%5.2f',100*(err_count(5)+err_count(7))/kMC));
                set(handles.text_clash_fail,'String',sprintf('%5.2f',100*err_count(11)/err_count(1)));
                set(handles.text_xlink_fail,'String',sprintf('%5.2f',100*err_count(10)/err_count(1)));
                if closed_loop,
                    set(handles.text_SANS_fail,'String',sprintf('%5.2f',100*err_count(2)/kMC));
                    set(handles.text_SAXS_fail,'String',sprintf('%5.2f',100*err_count(4)/kMC));
                end;
                handles.text_time_per_model.String = sprintf('%4.0f',runtime/success);
            end
            drawnow;
        end;
        if success >= max_models,
            break
        end;
    end;
    if success >= max_models,
        break
    end;
%     if mod(kMC,10000)
%         fprintf(1,'%6.3f ms per trial\n',1000*runtime/kMC);
%     end
    if runtime >= max_seconds,
        add_msg_board('Warning: Ensemble generation stopped as maximum allotted time was exceeded.');
        break
    end;
end;

if success > max_models
    success = max_models;
end

diagnostics.resax = resax;
success_distr = (1-success/kMC)*res_stat/sum(res_stat);
if reverse,
    success_distr = fliplr(success_distr);
end;
normalize = 1- cumsum(success_distr);
normalize(normalize==0)=1;
for ks = 2:length(success_distr),
    success_distr(ks) = success_distr(ks)/normalize(ks-1);
end;
if reverse,
    success_distr = fliplr(success_distr);
end;
diagnostics.rax = 10*rax;

diagnostics.success_distr = success_distr;

diagnostics.runtime = runtime;

if ~exist('model_indices','var')
    diagnostics.success = 0;
    diagnostics.snum = 0;
    diagnostics.time_per_model = inf;
    return
end
delete_empty_models(model_indices(1),success);
diagnostics.success = success;
diagnostics.time_per_model = runtime/success;

diagnostics.distributions = distributions;
diagnostics.restraint_distr = restraint_distr;
diagnostics.descriptors = descriptors;
diagnostics.monitor_distr = monitor_distr;
diagnostics.monitor_descr = monitor_descr;
diagnostics.snum_vec = snum_vec;
diagnostics.snum = snum;

if kMC == ntrials,
    fprintf(fid_report,'Warning: Ensemble generation stopped as maximum allotted number of Monte Carlo trials was exceeded.\n');
end;
hours = floor(runtime/3600);
minutes = floor((runtime-3600*hours)/60);
seconds = round(runtime - 3600*hours - 60*minutes);
if success > 0
    fprintf(fid_report,'Generation of %i models took %i h %i min %i s\n',success,hours,minutes,seconds);
    fprintf(fid_report,'%4.0f s/model\n',runtime/success);
    fprintf(fid_report,'Models were stored in structure %i and in PDB files starting with %s.\n',diagnostics.snum,fname);
else
    fprintf(fid_report,'No model found within %i h %i min %i s\n',hours,minutes,seconds);
end;
fprintf(fid_report,'%5.2f%% of allowed number of trials was used.\n',100*kMC/ntrials);
fprintf(fid_report,'%i backbone models were generated.\n',err_count(1));
fprintf(fid_report,'%5.2f%% of all trials failed by restraint violation.\n',100*err_count(6)/kMC);
fprintf(fid_report,'%5.2f%% of all trials failed by internal loop clashes.\n',100*(err_count(3)+err_count(8))/kMC);
fprintf(fid_report,'%5.2f%% of all trials failed by clashes of the loop with the rigid bodies.\n',100*(err_count(5)+err_count(7))/kMC);
fprintf(fid_report,'%5.2f%% of all backbones failed by sidechain clashes within the loop.\n',100*err_count(11)/err_count(1));
fprintf(fid_report,'%5.2f%% of all backbones failed by sidechain clashes with rigid bodies.\n',100*err_count(10)/err_count(1));
if closed_loop,
    fprintf(fid_report,'%5.2f%% of all trials did not reach terminal anchor residue.\n',100*err_count(2)/kMC);
    fprintf(fid_report,'%5.2f%% of all trials failed Ramachandran angles at terminal anchor.\n',100*err_count(4)/kMC);
end;

if options.n_restraints > 0
    fprintf(fid_report,'--- Distribution shifts and overlaps ---\n\n');
    fprintf(fid_report,'Shift [Å]\tOverlap\n');
    for kd = 1:options.n_restraints,
        [overlap,shift] = get_overlap(diagnostics.rax,diagnostics.distributions{kd},diagnostics.restraint_distr{kd});
        fprintf(fid_report,'%5.2f\t%5.3f\n',shift,overlap);
    end;
end
fclose(fid_report);

if interactive

    if kMC == ntrials,
        add_msg_board('Warning: Ensemble generation stopped as maximum allotted number of Monte Carlo trials was exceeded.');
    end;
    set(handles.text_percent_time,'String',sprintf('%5.2f',100*runtime/max_seconds));
    set(handles.text_max_trials,'String',sprintf('%5.2f',100*kMC/ntrials));
    set(handles.text_success,'String',sprintf('%i',success));
    set(handles.text_dmg_fail,'String',sprintf('%i',err_count(1)));
    set(handles.text_auxiliary_fail,'String',sprintf('%5.2f',100*err_count(6)/kMC));
    set(handles.text_core_fail,'String',sprintf('%5.2f',100*(err_count(3)+err_count(8))/kMC));
    set(handles.text_linker_fail,'String',sprintf('%5.2f',100*(err_count(5)+err_count(7))/kMC));
    set(handles.text_clash_fail,'String',sprintf('%5.2f',100*err_count(11)/err_count(1)));
    set(handles.text_xlink_fail,'String',sprintf('%5.2f',100*err_count(10)/err_count(1)));
    if closed_loop,
        set(handles.text_SANS_fail,'String',sprintf('%5.2f',100*err_count(2)/kMC));
        set(handles.text_SAXS_fail,'String',sprintf('%5.2f',100*err_count(4)/kMC));
    end;
    handles.text_time_per_model.String = sprintf('%4.0f',runtime/success);
end
drawnow;

function loopname = write_pdb_backbone(coor,sequence,fname0,model,res1,directory)

my_dir = pwd;
cd(directory);

residues='ALAARGASNASPCYSGLNGLUGLYHISILELEULYSMETPHEPROSERTHRTRPTYRVAL';
oneletter='ARNDCQEGHILKMFPSTWYV';

fname = sprintf('%s_m%i',fname0,model);

loopname = [fname '.pdb'];
wfile=fopen(loopname,'w');
for k = 1:length(sequence),
    respoi=strfind(oneletter,sequence(k));
    residue=residues(1+3*(respoi-1):3+3*(respoi-1));
    N = coor(4*k-3,:);
    CA = coor(4*k-2,:);
    C = coor(4*k-1,:);
    O = coor(4*k,:);
    fprintf(wfile,'%s%3i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           N\n','ATOM    ',4*k-3,'  N   ',residue,k+res1-1,N(1),N(2),N(3));
    fprintf(wfile,'%s%3i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM    ',4*k-2,'  CA  ',residue,k+res1-1,CA(1),CA(2),CA(3));
    fprintf(wfile,'%s%3i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM    ',4*k-1,'  C   ',residue,k+res1-1,C(1),C(2),C(3));
    fprintf(wfile,'%s%3i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           O\n','ATOM    ',4*k,'  O   ',residue,k+res1-1,O(1),O(2),O(3));
end;
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

if ~exist('min_approach','var'),
    min_approach = 2.00; % ### should be 2.0
end;

approach_prot = -1;
approach_loop = -1;
pclash = 1;
iclash = 1;
loop_coor = zeros(5000,3);
l_res_assign = zeros(5000,3);
fid=fopen(loopname);
if fid==-1,
    add_msg_board(sprintf('Warning: Loop structure PDB file %s missing. Rejected.',loopname));
    cd(my_dir);
    return;
end;
poi = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if length(tline) >= 6
        record=tline(1:6);
        resnum = str2double(tline(23:26));
        if strcmpi(record,'ATOM  ') || strcmpi(record,'HETATM'),
            if length(tline) < 78 || tline(78)~='H',
                if ~strcmpi(tline(14),'H') && resnum ~= res1 && resnum ~= resend,
                    poi = poi +1;
                    l_res_assign(poi) = resnum;
                    valstr = [tline(31:38) ' ' tline(39:46) ' ' tline(47:54)];
                    loop_coor(poi,:) = str2num(valstr);
                end;
            end;
        end;
    end;
end;
fclose(fid);

cd(my_dir);

loop_coor = loop_coor(1:poi,:);
l_res_assign = l_res_assign(1:poi);

pclash = 0;
iclash = 0;

[m1,~] = size(loop_coor); % get sizes of the coordinates arrays
[m2,~] = size(prot_coor);

if m2 > 0,
    a2 = repmat(sum(loop_coor.^2,2),1,m2);
    b2 = repmat(sum(prot_coor.^2,2),1,m1).';
    pair_dist = sqrt(abs(a2 + b2 - 2*loop_coor*prot_coor.'));
    min_dist = min(min(pair_dist));

    approach_prot = min_dist;
    if min_dist < min_approach,
    %    fprintf(2,'Minimum sidegroup distance to protein is %6.2f Å\n',min_dist);
       pclash = 1;
       cd(my_dir);
       return
    end;
end;

min_dist = 1e6;
% test for minimum distance within loop
for k1 = 1:poi-1,
    for k2 = k1+1:poi,
        if abs(l_res_assign(k1)-l_res_assign(k2))>1,
            approach = norm(loop_coor(k1,:) - loop_coor(k2,:));
            if approach < min_dist,
                min_dist = approach;
            end;
        end;
    end;
end;
approach_loop = min_dist;
if min_dist < min_approach,
%     fprintf(2,'Minimum distance of two heavy atoms in loop with sidegroups is: %6.2f Å\n',min_dist);
    iclash = 1;
end;

function delete_empty_models(snum,success)

global model

chains=length(model.structures{snum});
for cnum=1:chains,
    for k=1:success,
        atoms{k} = model.structures{snum}(cnum).atoms{k};
    end;
    model.structures{snum}(cnum).atoms = atoms;
    for k=1:success,
        residues{k} = model.structures{snum}(cnum).residues{k};
    end;
    model.structures{snum}(cnum).residues = residues;    
    for k=1:success,
        xyz{k} = model.structures{snum}(cnum).xyz{k};
    end;
    model.structures{snum}(cnum).xyz = xyz;
    for k=1:success,
        Bfactor{k} = model.structures{snum}(cnum).Bfactor{k};
    end;
    model.structures{snum}(cnum).Bfactor = Bfactor;
    for k=1:success,
        Btensor{k} = model.structures{snum}(cnum).Btensor{k};
    end;
    model.structures{snum}(cnum).Btensor = Btensor;
end;

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
        
function [distributions,restraint_distr,descriptors,monitor_distr,monitor_descr] = mk_report_distributions(fid_report,indices,restrain,monitor,p_model,distributions,restraint_distr,descriptors,monitor_distr,monitor_descr,res1)

adr = mk_address(indices);

rpr = zeros(1,length(restrain));
for k = 1:length(restrain)
    rpr(k) = rpr(k) + length(restrain(k).r_beacon);
    rpr(k) = rpr(k) + length(monitor(k).r_beacon);
    rpr(k) = rpr(k) + length(restrain(k).r_intern);
    rpr(k) = rpr(k) + length(monitor(k).r_intern);
    rpr(k) = rpr(k) + length(restrain(k).depth);
    rpr(k) = rpr(k) + length(monitor(k).depth);
    rpr(k) = rpr(k) + length(restrain(k).oligomer);
    rpr(k) = rpr(k) + length(monitor(k).oligomer);
end

if sum(rpr) == 0 % there are no restraints and no potential restraints to be monitored
    return
end

fprintf(fid_report,'\n--- Restraint fulfillment for model %i with weight %5.3f ---\n\n',indices(3),p_model);

poi = 0;
dispoi = 0;
mdispoi = 0;
for k = 1:length(restrain),
    resnum = res1+k-1;
    radr = sprintf('%s%i',adr,resnum);
    if rpr(k) > 0 % make section for a residue only if it has restraints or potential restraints to be monitored
        fprintf(fid_report,'\nResidue %i\n\n',resnum);
    end;
    rindices = resolve_address(radr);
    for kr = 1:length(restrain(k).r_beacon),
        NO_pos1 = get_NO_pos(rindices,restrain(k).r_beacon(kr).label1,298);
        % [rmsd2,xyz1]=NOpos_rmsd(NO_pos1);
        NO_pos2 = get_NO_pos(restrain(k).r_beacon(kr).bindices,restrain(k).r_beacon(kr).label2,298);
        % [rmsd2,xyz2]=NOpos_rmsd(NO_pos2);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).r_beacon(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
%                 figure(poi); clf;
                rarg = (rax -restrain(k).r_beacon(kr).par1)/restrain(k).r_beacon(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Beacon restraint %i-%s with <r_{exp}> = %5.2f',resnum,restrain(k).r_beacon(kr).resb,restrain(k).r_beacon(kr).par1);
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Beacon restraint [%5.2f, %5.2f] Å to residue %s fulfilled at <r> = %5.2f\n',...
                    restrain(k).r_beacon(kr).par1,restrain(k).r_beacon(kr).par2,restrain(k).r_beacon(kr).resb,rmean);
            case 'bounds'
                if rmean >= restrain(k).r_beacon(kr).par1 && rmean <= restrain(k).r_beacon(kr).par2,
                    fprintf(fid_report,'Beacon restraint [%5.2f, %5.2f] Å to residue %s fulfilled at <r> = %5.2f\n',...
                        restrain(k).r_beacon(kr).par1,restrain(k).r_beacon(kr).par2,restrain(k).r_beacon(kr).resb,rmean);
                else
                    fprintf(fid_report,'Beacon restraint [%5.2f, %5.2f] Å to residue %s violated at <r> = %5.2f\n',...
                        restrain(k).r_beacon(kr).par1,restrain(k).r_beacon(kr).par2,restrain(k).r_beacon(kr).resb,rmean);
                end;
        end;
    end;

    for kr = 1:length(monitor(k).r_beacon),
        NO_pos1 = get_NO_pos(rindices,monitor(k).r_beacon(kr).label1,298);
        NO_pos2 = get_NO_pos(monitor(k).r_beacon(kr).bindices,monitor(k).r_beacon(kr).label2,298);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Beacon monitor %i-%s with <r_{sim}> = %5.2f',resnum,monitor(k).r_beacon(kr).resb,rmean);
        fprintf(fid_report,'Monitored distance to residue %s is <r> = %5.2f\n',monitor(k).r_beacon(kr).resb,rmean);
    end;

    for kr = 1:length(restrain(k).r_intern),
        r2adr = sprintf('%s%i',adr,restrain(k).r_intern(kr).resb);
        r2indices = resolve_address(r2adr);
        NO_pos1 = get_NO_pos(rindices,restrain(k).r_intern(kr).label1,298);
        % [rmsd2,xyz1]=NOpos_rmsd(NO_pos1);
        NO_pos2 = get_NO_pos(r2indices,restrain(k).r_intern(kr).label2,298);
        % [rmsd2,xyz2]=NOpos_rmsd(NO_pos2);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).r_intern(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
%                 figure(poi); clf;
                rarg = (rax -restrain(k).r_intern(kr).par1)/restrain(k).r_intern(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Internal restraint %i-%s with <r_{exp}> = %5.2f',resnum,restrain(k).r_intern(kr).resb,restrain(k).r_intern(kr).par1);

%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] Å to residue %s fulfilled at <r> = %5.2f\n',...
                    restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
            case 'bounds'
                if rmean >= restrain(k).r_intern(kr).par1 && rmean <= restrain(k).r_intern(kr).par2,
                    fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] Å to residue %s fulfilled at <r> = %5.2f\n',...
                        restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
                else
                    fprintf(fid_report,'Internal restraint [%5.2f, %5.2f] Å to residue %s violated at <r> = %5.2f\n',...
                        restrain(k).r_intern(kr).par1,restrain(k).r_intern(kr).par2,restrain(k).r_intern(kr).resb,rmean);
                end;
        end;
    end;
    
    for kr = 1:length(monitor(k).r_intern),
        r2adr = sprintf('%s%i',adr,monitor(k).r_intern(kr).resb);
        r2indices = resolve_address(r2adr);
        NO_pos1 = get_NO_pos(rindices,monitor(k).r_intern(kr).label1,298);
        NO_pos2 = get_NO_pos(r2indices,monitor(k).r_intern(kr).label2,298);
        [rax,sim_distr]=get_distribution(NO_pos1,NO_pos2,0.05);
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Internal restraint %i-%i',resnum,restrain(k).r_intern(kr).resb);
        fprintf(fid_report,'Monitored distance to residue %i is <r> = %5.2f\n',monitor(k).r_intern(kr).resb,rmean);
    end;

    for kr = 1:length(restrain(k).depth),
        NO_pos = get_NO_pos(rindices,restrain(k).depth(kr).label,298);
        [rax,sim_distr]=get_distribution_z(NO_pos,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).depth(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                rarg = (rax -restrain(k).r_intern(kr).par1)/restrain(k).r_intern(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Depth restraint %i with <r_{exp}> = %5.2f',resnum,restrain(k).r_depth(kr).par1);

%                 figure(poi); clf;
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] Å fulfilled at <r> = %5.2f\n',...
                    restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
            case 'bounds'
                if rmean >= restrain(k).depth(kr).par1 && rmean <= restrain(k).depth(kr).par2,
                    fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] Å fulfilled at <r> = %5.2f\n',...
                        restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
                else
                    fprintf(fid_report,'Depth restraint [%5.2f, %5.2f] Å violated at <r> = %5.2f\n',...
                        restrain(k).depth(kr).par1,restrain(k).depth(kr).par2,rmean);
                end;
        end;
    end;

    for kr = 1:length(monitor(k).depth),
        NO_pos = get_NO_pos(rindices,monitor(k).depth(kr).label,298);
        [rax,sim_distr]=get_distribution_z(NO_pos,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Monitored depth %i is <r_{sim}> = %5.2f',resnum,rmean);
        fprintf(fid_report,'Monitored depth is <r> = %5.2f\n',rmean);
    end;

    for kr = 1:length(restrain(k).oligomer),
        NO_pos = get_NO_pos(rindices,restrain(k).oligomer(kr).label,298);
        [rax,sim_distr]=get_distribution_oligomer(NO_pos,restrain(k).oligomer(kr).n,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        switch restrain(k).oligomer(kr).type
            case 'Gaussian'
                poi = poi + 1;
                dispoi = dispoi + 1;
                distributions{dispoi} = distributions{dispoi} + p_model*sim_distr;
                rarg = (rax -restrain(k).oligomer(kr).par1)/restrain(k).oligomer(kr).par2;
                exp_distr = exp(-rarg.^2);
                exp_distr = exp_distr/sum(exp_distr);
                restraint_distr{dispoi} = restraint_distr{dispoi} + p_model*exp_distr;
                descriptors{dispoi} = sprintf('Oligomer restraint %i (n = %i) with <r_{exp}> = %5.2f and <r_{sim}> = %5.2f',resnum,restrain(k).oligomer(kr).n,restrain(k).oligomer(kr).par1,rmean);
%                 figure(poi); clf;
%                 plot(rax,exp_distr,'k');
%                 hold on;
%                 plot(rax,sim_distr,'r');
%                 ma = max([max(exp_distr),max(sim_distr)]);
%                 axis([0,80,-0.1*ma,1.1*ma]);
                fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] Å fulfilled at <r> = %5.2f\n',...
                    restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
            case 'bounds'
                if rmean >= restrain(k).oligomer(kr).par1 && rmean <= restrain(k).oligomer(kr).par2,
                    fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] Å fulfilled at <r> = %5.2f\n',...
                        restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
                else
                    fprintf(fid_report,'Oligomer restraint [%5.2f, %5.2f] Å violated at <r> = %5.2f\n',...
                        restrain(k).oligomer(kr).par1,restrain(k).oligomer(kr).par2,rmean);
                end;
        end;
    end;

    for kr = 1:length(monitor(k).oligomer),
        NO_pos = get_NO_pos(rindices,monitor(k).oligomer(kr).label,298);
        [rax,sim_distr]=get_distribution_oligomer(NO_pos,monitor(k).oligomer(kr).n,0.05);        
        sim_distr = sim_distr/sum(sim_distr);
        rax = 10*rax;
        rmean = sum(rax.*sim_distr);
        mdispoi = mdispoi + 1;
        monitor_distr{mdispoi} = monitor_distr{mdispoi} + p_model*sim_distr;
        monitor_descr{mdispoi} = sprintf('Monitored oligomer distance %i (n = %i) is <r_{sim}> = %5.2f',resnum,monitor(k).oligomer(kr).n,rmean);
        fprintf(fid_report,'Monitored oligomer distance is <r> = %5.2f\n',rmean);
    end;

end;

function NO_pos = get_NO_pos(indices,label,T)

global model
global label_defs
global hMain

if strcmpi(label,'CA'),
    adr = sprintf('%s.CA',mk_address(indices));
    [~,xyz] = get_object(adr,'coor');
    NO_pos = [xyz 1];
    return
end;

NO_pos = [];
if isfield(model,'sites'),
    for k0=1:length(model.sites),
        for k1=1:length(model.sites{k0}),
            for k=1:length(model.sites{k0}(k1).residue),
                if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0,
                    id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                    if strcmpi(label,label_defs.residues(id).short_name),
                        if T == model.sites{k0}(k1).residue(k).T,
                            NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

if isempty(NO_pos),
    adr = mk_address(indices);
    command=sprintf('rotamers %s %s %i',adr,label,T);
    hMain.store_undo=false;
    hMain.dynamic_rotamers=false;
    cmd(hMain,command);
end;

for k0=1:length(model.sites),
    for k1=1:length(model.sites{k0}),
        for k=1:length(model.sites{k0}(k1).residue),
            if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0,
                id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                if strcmpi(label,label_defs.residues(id).short_name),
                    if T == model.sites{k0}(k1).residue(k).T,
                        NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                    end;
                end;
            end;
        end;
    end;
end;


function [overlap,shift] = get_overlap(rax,distr1,distr2)

distr1 = distr1/sum(distr1);
distr2 = distr2/sum(distr2);
shift = sum(rax.*distr1) - sum(rax.*distr2); 
% fprintf(1,'Sim. mean distance: %5.2f\n',sum(rax.*distr1));
% fprintf(1,'Restraint mean distance: %5.2f\n',sum(rax.*distr2));
diff = distr1-distr2;
overlap = 1 - sqrt(sum(diff.^2));
    
