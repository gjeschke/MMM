function [secdefs,RNA,ensemble,msg] = process_rna_domain_restraints(restraints,rba)
%
% restraints    restraint structure as read in by rd_restraints_rigiflex,
%               only the field RNA is used, at least one binding motif must
%               be defined
% rba           structure number and model number of the rigid-body
%               arrangement
%
% secdefs       section definitions, array of structures with fields
%               .nts        first and last nucleotide of the segment
%               .anchori    [6,3] array of initial anchor coordinates
%               .anchore    [3,3] array of target anchor coordinates
%               .ntoffset   offset for nucleotide numbers
% RNA           definition of the RNA, structure with fields
%               .sequence   sequence string
%               .nta        number of first nucleotide
% ensemble      number of models in the ensemble, optional, defaults to 1
%
% G. Jeschke, 02.05.2018


msg.error = 0;
msg.text = 'OK';

secdefs = [];
RNA = [];
ensemble = 1;

if ~isfield(restraints,'RNA')
    msg.error = 1;
    msg.text = 'No RNA restraints defined. Aborting.';
    return
end

if ~isfield(restraints.RNA,'bind') || isempty(restraints.RNA.bind)
    msg.error = 2;
    msg.text = 'No RNA binding motif defined. Aborting.';
    return
end

RNA = restraints.RNA;

if isfield(RNA,'models')
    ensemble = RNA.models;
end

if isfield(restraints,'stemlibs')
    for kb = 1:length(RNA.bind)
        achain1 = sep_anchor_address(RNA.bind(kb).anchora);
        achain2 = sep_anchor_address(RNA.bind(kb).anchore);
        if ~strcmp(achain1,achain2)
            msg.error = 3;
            msg.text = sprintf('RNA binding motif %i has inconsistent chain information. Aborting.',kb);
            return
        end
        for kl = 1:length(restraints.stemlibs)
            if strcmpi(restraints.stemlibs{kl}.chaintag,achain1)
                RNA.bind(kb).nta = restraints.stemlibs{kl}.nta;
                RNA.bind(kb).nte = restraints.stemlibs{kl}.nte;
                RNA.bind(kb).anchora = sprintf('%s%i',restraints.stemlibs{kl}.chaintag,restraints.stemlibs{kl}.nta);
                RNA.bind(kb).anchore = sprintf('%s%i',restraints.stemlibs{kl}.chaintag,restraints.stemlibs{kl}.nte);
            end
        end
    end
end

[restrain,RNA,~,~,~,~,msg] = process_RNA_restraints(RNA);
if msg.error
    return
end

seqlen = length(RNA.sequence);

% analyze topology (only forward segments)
% the first forward segment starts at the end of the first binding motif


knt = 1;
while isempty(restrain(knt).motif)
    knt = knt + 1;
end

while ~isempty(restrain(knt+1).motif)
    knt = knt + 1;
end
fwdpoi = knt;

type = -1;
if restrain(fwdpoi + 1).stem
    stem_start = fwdpoi+1;
else
    loop_start = fwdpoi + 1;
end

% find all sections (between anchor points in rigid bodies) and subsections
% (stem or loop sections) forward from the first binding motif
section = 1;
sections_fwd{1} = zeros(100,3);
segment = 0;
while fwdpoi < seqlen
    fwdpoi = fwdpoi + 1;
    newtype = restrain(fwdpoi).stem;
    if type == 1 && (newtype == 0 || restrain(fwdpoi).known)
        segment = segment + 1;
        sections_fwd{section}(segment,:) = [type,restrain(stem_start).nt,restrain(fwdpoi-1).nt];
        % fprintf(1,'Inserting stem from nucleotide %i to %i\n',restrain(stem_start).nt,restrain(fwdpoi-1).nt);
        restrain = consolidate_stem(restrain,RNA,restrain(stem_start).stem_id);
        loop_start = fwdpoi;
    end
    if type == 0 && (newtype == 1 || restrain(fwdpoi).known)
        segment = segment + 1;
        sections_fwd{section}(segment,:) = [type,restrain(loop_start).nt,restrain(fwdpoi-1).nt];
        % fprintf(1,'Inserting loop from nucleotide %i to %i\n',restrain(loop_start).nt,restrain(fwdpoi-1).nt);
        stem_start = fwdpoi;
    end
    if type ~= 2 && restrain(fwdpoi).known
        sections_fwd{section} = sections_fwd{section}(1:segment,:);
        section = section + 1;
        segment = 0;
        sections_fwd{section} = zeros(100,3);
    end
    if type == 2 && ~restrain(fwdpoi).known
        if newtype == 0
            loop_start = fwdpoi;
        else
            stem_start = fwdpoi;
        end
    end
    type = newtype;
    if restrain(fwdpoi).known
        type = 2;
    end
end
if type == 0
    segment = segment + 1;
    % fprintf(1,'Inserting loop from nucleotide %i to %i\n',restrain(loop_start).nt,restrain(fwdpoi).nt);
    sections_fwd{section}(segment,:) = [type,restrain(loop_start).nt,restrain(fwdpoi).nt];
elseif type == 1
    segment = segment + 1;
    % fprintf(1,'Inserting stem from nucleotide %i to %i\n',restrain(stem_start).nt,restrain(fwdpoi).nt);
    restrain = consolidate_stem(restrain,RNA,restrain(stem_start).stem_id);
    sections_fwd{section}(segment,:) = [type,restrain(stem_start).nt,restrain(fwdpoi).nt];
end

if segment > 0
    sections_fwd{section} = sections_fwd{section}(1:segment,:);
else
    sections_fwd = sections_fwd(1:section-1);
end

% assemble the information for modelling of all forward subsections
for ks = 1:length(sections_fwd)
    csection = sections_fwd{ks};
    [mss,~] = size(csection); % number of subsections
    for kss = 1:mss
        secdefs(ks,kss).nts = csection(kss,2:3);
        if kss == 1 % there is an initial anchor in a binding motif
            secdefs(ks,kss).anchori = get_initial_anchor(csection(kss,2),RNA,rba);
        else
            secdefs(ks,kss).anchori = []; % anchor must be resolved during modelling
        end
        if kss == mss % there is (possibly) a target anchor in a binding motif
            secdefs(ks,kss).anchore = get_target_anchor(csection(kss,3),RNA,rba);
        else
            secdefs(ks,kss).anchore = []; % anchor can possibly be resolved during modelling
        end
        secdefs(ks,kss).ntoffset = csection(kss,2)-1;
        if csection(kss,1) == 1 % only for stems
            [lib,libtype,defines] = get_stem_lib(csection(kss,2:3),RNA);
            secdefs(ks,kss).lib = lib;
            secdefs(ks,kss).libtype = libtype;
            secdefs(ks,kss).defines = defines;
            defined_by = stem_defined(secdefs,ks,kss,sections_fwd);
            if ~isempty(defined_by)
                secdefs(ks,kss).defined_by = [0 defined_by];
            else
                secdefs(ks,kss).defined_by = [];
            end
        end
    end
end

poi = 0;
m = length(secdefs);
for k = 1:m
    if secdefs(k).nts(end) < RNA.nte
        poi = poi + 1;
    end
end
secdefs = secdefs(1:poi);

function [restrain,restraints,monitor,cancelled,number,number_monitor,msg] = process_RNA_restraints(restraints,stemlibs)
% Processes restraints for a single flexible domain

msg.error = 0;
msg.text = 'Restraints processed';

cancelled=false;

restrain = [];
monitor = [];
number = 0;
number_monitor = 0;

if ~isfield(restraints,'sequence')
    msg.text = 'RNA sequence specification is missing.'; 
    restrain=[];
    msg.error = 11;
    return;
end;

for k = 1:length(restraints.sequence)
    restrain(k).stem = 0;
    restrain(k).stem_id = [];
    restrain(k).label = [];
    restrain(k).r_beacon = [];
    restrain(k).r_intern = [];
    restrain(k).oligomer = [];
    restrain(k).panchor = [];
    restrain(k).nanchor = [];
    restrain(k).nt = restraints.nta + k -1;
    restrain(k).motif = [];
    restrain(k).known = 0;
end

monitor = restrain;


% binding motif information, there is at least one binding motif

for k = 1:length(restraints.bind)
    ind1 = resolve_address(restraints.bind(k).anchora);
    ind2 = resolve_address(restraints.bind(k).anchore);
    if length(ind1)~= 4 || length(ind2)~=4 || sum(abs(ind1(1:3)-ind2(1:3)))
        msg.text = sprintf('Binding motif specification %s-%s cannot be resolved.',restraints.bind(k).anchora,restraints.bind(k).anchore);
        msg.error = 17;
        restrain=[];
        return;
    end
    m = ind2(4) - ind1(4);
    if m ~= restraints.bind(k).nte - restraints.bind(k).nta
        msg.text = sprintf('Binding motif %s does not match RNA range %i to %i.',adr,restraints.bind(k).nta,restraints.bind(k).nte);
        msg.error = 18;
        restrain=[];
        return;
    end
    indc = ind1;
    for knt = 0:m
        indc(4) = ind1(4) + knt;
        restrain(restraints.bind(k).nta+knt-restraints.nta+1).motif = indc;
        restrain(restraints.bind(k).nta+knt-restraints.nta+1).known = 1;
    end
end

function NO = get_relative_label(libname)

load(libname);
midNO = rot_lib.usefull_atoms.midNO;
pops = rot_lib.calibration.pop;
NO = zeros(1,3);
for k = 1:length(rot_lib.library)
    coor = rot_lib.library(k).ecoor;
    NO = NO + pops(k)*(coor(midNO(1),2:4) + coor(midNO(2),2:4))/2;
end
NO = NO/(sum(pops));

function res = correct_section_address(adr,nta)

if adr(1) == 'R'
    resstr = adr(2:end);
else
    resstr = adr;
end;
res = str2double(resstr) + nta - 1;
if isnan(res),
    res = [];
end;
if res - floor(res) > eps,
    res = [];
end;

function [restrain,number] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,rmean,sigr,res1,number,label1,label2,bindices,resb)

scale_units = 1; % rd_restraints_rigiflex already converts to Å

grace = 0.5; % 5 Å uncertainty of label position

k = res_loop - res1 + 1;
kr = length(restrain(k).r_beacon)+1;
restrain(k).label = NO;
restrain(k).r_beacon(kr).xyz = xyz_beacon;
restrain(k).r_beacon(kr).label1 = label1;
restrain(k).r_beacon(kr).label2 = label2;
restrain(k).r_beacon(kr).bindices = bindices;
restrain(k).r_beacon(kr).resb = resb;
if rmean > 0 && sigr > 0
    restrain(k).r_beacon(kr).type = 'Gaussian';
    restrain(k).r_beacon(kr).par1 = rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = sqrt(sigr^2 + grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_beacon(kr).type = 'bounds';
    restrain(k).r_beacon(kr).par1 = -rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = -sigr*scale_units;
end;

function [restrain,number] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,rmean,sigr,res1,number,label1,label2)

scale_units = 1; % rd_restraints_rigiflex already converts to Å

grace = 0.5; % 5 Å uncertainty of label position

if resa < resb, % restraint must be stored at the later site
    exch = resa;
    resa = resb;
    resb = exch;
    exch = NO1;
    NO1 = NO2;
    NO2 = exch;
end;
k = resa - res1 + 1;
k2 = resb - res1 + 1;
kr = length(restrain(k).r_intern)+1;
restrain(k).label = NO1;
restrain(k2).label = NO2;
restrain(k).r_intern(kr).site = k2;
restrain(k).r_intern(kr).label1 = label1;
restrain(k).r_intern(kr).label2 = label2;
restrain(k).r_intern(kr).resb = resb;
if rmean > 0 && sigr > 0
    restrain(k).r_intern(kr).type = 'Gaussian';
    restrain(k).r_intern(kr).par1 = rmean*scale_units;
    restrain(k).r_intern(kr).par2 = sqrt(sigr^2 + 2*grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_intern(kr).type = 'bounds';
    restrain(k).r_intern(kr).par1 = -rmean*scale_units;
    restrain(k).r_intern(kr).par2 = -sigr*scale_units;
end

function [restrain,number] = mk_oligomer_restraint(restrain,NO,res,n,rmean,sigr,res1,number,label)

scale_units= 1; % rd_restraints_rigiflex already converts to Å

k = res - res1 + 1;
kr = length(restrain(k).oligomer)+1;
restrain(k).oligomer(kr).label_type = label;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).oligomer(kr).site = 'label';
else
    restrain(k).oligomer(kr).site = 'CA';
end;
if rmean > 0 && sigr > 0,
    restrain(k).oligomer(kr).type = 'Gaussian';
    restrain(k).oligomer(kr).par1 = rmean*scale_units;
    restrain(k).oligomer(kr).par2 = sigr*scale_units;
    number = number + 1;
else
    restrain(k).oligomer(kr).type = 'bounds';
    restrain(k).oligomer(kr).par1 = -rmean*scale_units;
    restrain(k).oligomer(kr).par2 = -sigr*scale_units;
end;
restrain(k).oligomer(kr).n = n;

function restrain = consolidate_stem(restrain,restraints,stem_id)

for kn = restraints.stems(stem_id).C5pf:restraints.stems(stem_id).C3pf
    restrain(kn - restraints.nta + 1).known = 1;
end
for kn = restraints.stems(stem_id).C5pb:restraints.stems(stem_id).C3pb
    restrain(kn - restraints.nta + 1).known = 1;
end

function [lib,libtype,defines] = get_stem_lib(nts,RNA)

lib = '?';
libtype = [0,0];

for k = 1:length(RNA.stems)
    if RNA.stems(k).C5pf == nts(1)
        libtype(1) = 1;
        lib = RNA.stems(k).lib;
        defines = [RNA.stems(k).C5pb,RNA.stems(k).C3pb];
    end
    if RNA.stems(k).C3pf == nts(2)
        libtype(2) = 1;
    end
    if RNA.stems(k).C5pb == nts(1)
        libtype(1) = 2;
        lib = RNA.stems(k).lib;
        defines = [RNA.stems(k).C5pf,RNA.stems(k).C3pf];
    end
    if RNA.stems(k).C3pb == nts(2)
        libtype(2) = 2;
    end
end

function anchor = get_target_anchor(nt,RNA,rba)

anchor = [];

for k = 1:length(RNA.bind)
    stag = mk_address_parts(rba(1));
    if nt+1 == RNA.bind(k).nta % target anchor
        ntadr = RNA.bind(k).anchora;
        poi = strfind(ntadr,')');
        chain = ntadr(1:poi);
        ntnum = ntadr(poi+1:end);        
        indices = resolve_address(sprintf('[%s]%s{%i}%s',stag,chain,rba(2),ntnum));
        anchore = get_anchor_pseudo_torsion(indices);
        anchor = anchore(2:4,:);
    end
end

function anchor = get_initial_anchor(nt,RNA,rba)

anchor = [];

for k = 1:length(RNA.bind)
    stag = mk_address_parts(rba(1));
    if nt-1 == RNA.bind(k).nte % initial anchor
        ntadr = RNA.bind(k).anchore;
        poi = strfind(ntadr,')');
        chain = ntadr(1:poi);
        ntnum = ntadr(poi+1:end);
        indices = resolve_address(sprintf('[%s]%s{%i}%s',stag,chain,rba(2),ntnum));
        anchor = get_anchor_pseudo_torsion(indices,true);
    end
end

    
function defined_by = stem_defined(sec_def,ks,kss,sections)

defined_by = [];

nts = sec_def(ks,kss).nts;
for ls = 1:ks-1
    csection = sections{ls};
    [ms,~] = size(csection);
    for lss = 1:ms
        defines = sec_def(ls,lss).defines;
        if length(defines) == 2
            match = sum(nts-defines);
            if match ~= 0
                match = sum(fliplr(nts)-defines);
            end
            if match == 0
                defined_by = [ks,kss];
            end
        end
    end
end
for lss = 1:kss-1
    defines = sec_def(ks,lss).defines;
    if length(defines) == 2
        match = sum(nts-defines);
        if match ~= 0
            match = sum(fliplr(nts)-defines);
        end
        if match == 0
            defined_by = [ks,kss];
        end
    end
end

function terminal_anchor = get_terminal_anchor(cecoor,catomtags,lastnt)

terminal_anchor = 1e10*ones(4,3);
lastnti = find(cecoor(:,1) == lastnt);
lastntatoms = catomtags(lastnti);
for k = 1:length(lastnti)
    if strcmp(lastntatoms{k},'P')
        terminal_anchor(2,:) = cecoor(lastnti(k),2:4);
    end
    if strcmp(lastntatoms{k},'C4''')
        terminal_anchor(3,:) = cecoor(lastnti(k),2:4);
    end
    if strcmp(lastntatoms{k},'O3''')
        terminal_anchor(4,:) = cecoor(lastnti(k),2:4);
    end
end
prevnti = find(cecoor(:,1) == lastnt-1);
prevntatoms = catomtags(prevnti);
for k = 1:length(prevnti)
    if strcmp(prevntatoms{k},'C4''')
        terminal_anchor(1,:) = cecoor(prevnti(k),2:4);
    end
end
if max(max(terminal_anchor))>1e9
    terminal_anchor = [];
end

function ecoor = correct_section(ecoor,correction,environ,fixed,options)

ecoor = anchor_correction(ecoor,correction);

% try to repair clashes, if requested
if ~isempty(environ)
    coor = ecoor(:,2:4);
    % options.max_shift = 2;
    [coor,min_dist,cycles,max_shift] = clash_repair(coor,environ,fixed,options);
    ecoor(:,2:4) = coor;
    fprintf(1,'Clash repair required %i cycles and lead to a maximum shift of %6.2f Å\n',cycles,max_shift);
    if min_dist < 0.995*options.clash_thr || max_shift > 1.005*options.max_shift
        pair_dist = get_all_pair_dist(coor,environ);
        [mdv,mdp] = min(pair_dist);
        [min_dist,mdp2] = min(mdv);
        mdp1 = mdp(mdp2);
        fprintf(1,'RNA loop atom %i of nt %i clashed with environment atom %i at distance %6.3f Å.\n',mdp1,ecoor(mdp1,1),mdp2,min_dist);
        ecoor = [];
    end
end

function [chaintag,resnum]  = sep_anchor_address(address)

poi = strfind(address,')');
if ~isempty(poi)
    chaintag = address(1:poi);
else
    chaintag = '';
    poi = 0;
end
resnum = str2double(address(poi+1:end));
