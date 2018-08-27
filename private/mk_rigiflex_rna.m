function msg = mk_rigiflex_rna(restraints,rba,molecule,environ,ensemble,maxtime)
%
% restraints    restraint structure as read in by rd_restraints_rigiflex,
%               only the field RNA is used, at least one binding motif must
%               be defined
% rba           structure number and model number of the rigid-body
%               arrangement
% molecule      number of the RNA molecule to be generated, optional,
%               defaults to 1
% environ       coordinate array of environment atoms
% ensemble      number of models (per segment) in the ensemble, optional,
%               defaults to 1
% maxtime       maximum time per segment per model, optional, defaults to 
%               1 h
%
% G. Jeschke, 31.12.2017

clash_thr = 1.0;
max_nt_shift = 2;
max_shift_loop = 2.5;

msg.error = 0;
msg.text = 'RNA generated';

if ~exist('molecule','var')
    molecule = 1;
end

if ~exist('ensemble','var')
    ensemble = 1;
end

if ~exist('maxtime','var')
    maxtime = 1;
end

if ~isfield(restraints,'RNA')
    msg.error = 1;
    msg.text = 'No RNA restraints defined. Aborting.';
    return
end

if ~isfield(restraints.RNA(molecule),'bind') || isempty(restraints.RNA(molecule).bind)
    msg.error = 2;
    msg.text = 'No RNA binding motif defined. Aborting.';
    return
end

RNA = restraints.RNA(molecule);

[restrain,RNA,monitor,cancelled,number,number_monitor,msg] = process_RNA_restraints(RNA);
if msg.error
    return
end

seqlen = length(RNA.sequence);
offset = RNA.nta - 1;

% analyze topology
% the first forward segment starts at the end of the first binding motif
% the first backward segment starts at the beginning of the first binding
% motif or at the beginning of a stem section generated in the forward pass

knt = 1;
while isempty(restrain(knt).motif)
    knt = knt + 1;
end
backpoi = knt;

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

if backpoi > 1
% find all sections (between anchor points in rigid bodies) and subsections
% (stem or loop sections) backward from the first binding motif

    type = -1;
    if restrain(backpoi-1).known
        type = 2;
    elseif restrain(backpoi - 1).stem
        stem_start = backpoi-1;
    else
        loop_start = backpoi-1;
    end

end

section = 1;
sections_bck{1} = zeros(100,3);
segment = 0;
while backpoi > 1
    backpoi = backpoi - 1;
    newtype = restrain(backpoi).stem;
    if type == 1 && (newtype == 0 || restrain(backpoi).known)
        segment = segment + 1;
        sections_bck{section}(segment,:) = [type,restrain(stem_start).nt,restrain(backpoi+1).nt];
        % fprintf(1,'Inserting stem from nucleotide %i to %i\n',restrain(stem_start).nt,restrain(backpoi+1).nt);
        restrain = consolidate_stem(restrain,RNA,restrain(stem_start).stem_id);
        loop_start = backpoi;
    end
    if type == 0 && (newtype == 1 || restrain(backpoi).known)
        segment = segment + 1;
        sections_bck{section}(segment,:) = [type,restrain(loop_start).nt,restrain(backpoi+1).nt];
        % fprintf(1,'Inserting loop from nucleotide %i to %i\n',restrain(loop_start).nt,restrain(backpoi+1).nt);
        stem_start = backpoi;
    end
    if type ~= 2 && restrain(backpoi).known
        sections_bck{section} = sections_bck{section}(1:segment,:);
        section = section + 1;
        segment = 0;
        sections_bck{section} = zeros(100,3);
    end
    if type == 2 && ~restrain(backpoi).known
        if newtype == 0
            loop_start = backpoi;
        else
            stem_start = backpoi;
        end
    end
    type = newtype;
    if restrain(backpoi).known
        type = 2;
    end
end
if type == 0
    segment = segment + 1;
    % fprintf(1,'Inserting loop from nucleotide %i to %i\n',restrain(loop_start).nt,restrain(backpoi).nt);
    sections_bck{section}(segment,:) = [type,restrain(loop_start).nt,restrain(backpoi).nt];
elseif type == 1
    segment = segment + 1;
    % fprintf(1,'Inserting stem from nucleotide %i to %i\n',restrain(stem_start).nt,restrain(backpoi).nt);
    restrain = consolidate_stem(restrain,RNA,restrain(stem_start).stem_id);
    sections_bck{section}(segment,:) = [type,restrain(stem_start).nt,restrain(backpoi).nt];
end

if segment > 0
    sections_bck{section} = sections_bck{section}(1:segment,:);
else
    sections_bck = sections_bck(1:section-1);
end

% assemble the information for modelling of all forward subsections
for ks = 1:length(sections_fwd)
    csection = sections_fwd{ks};
    [mss,~] = size(csection); % number of subsections
    for kss = 1:mss
        secdef_fwd(ks,kss).nts = csection(kss,2:3);
        if kss == 1 % there is an initial anchor in a binding motif
            secdef_fwd(ks,kss).anchori = get_initial_anchor(csection(kss,2),restraints.RNA,rba);
        else
            secdef_fwd(ks,kss).anchori = []; % anchor must be resolved during modelling
        end
        if kss == mss % there is (possibly) a target anchor in a binding motif
            secdef_fwd(ks,kss).anchore = get_target_anchor(csection(kss,3),restraints.RNA,rba);
        else
            secdef_fwd(ks,kss).anchore = []; % anchor can possibly be resolved during modelling
        end
        secdef_fwd(ks,kss).ntoffset = csection(kss,2)-1;
        if csection(kss,1) == 1 % only for stems
            [lib,libtype,defines] = get_stem_lib(csection(kss,2:3),restraints.RNA);
            secdef_fwd(ks,kss).lib = lib;
            secdef_fwd(ks,kss).libtype = libtype;
            secdef_fwd(ks,kss).defines = defines;
            defined_by = stem_defined(secdef_fwd,ks,kss,sections_fwd);
            if ~isempty(defined_by)
                secdef_fwd(ks,kss).defined_by = [0 defined_by];
            else
                secdef_fwd(ks,kss).defined_by = [];
            end
        end
    end
end

% assemble the information for modelling of all backward subsections
for ks = 1:length(sections_bck)
    csection = sections_bck{ks};
    [mss,~] = size(csection); % number of subsections
    for kss = 1:mss
        secdef_bck(ks,kss).nts = csection(kss,2:3);
        if kss == 1 % there is a target anchor in a binding motif
            secdef_bck(ks,kss).anchore = get_target_anchor(csection(kss,2),restraints.RNA,rba);
        else
            secdef_bck(ks,kss).anchore = []; % anchor must be resolved during modelling
        end
        if kss == mss % there is (possibly) an initial anchor in a binding motif
            secdef_bck(ks,kss).anchori = get_initial_anchor(csection(kss,3),restraints.RNA,rba);
        else
            secdef_bck(ks,kss).anchori = []; % anchor can possibly be resolved during modelling
        end
        secdef_bck(ks,kss).ntoffset = csection(kss,2)-1;
        if csection(kss,1) == 1 % only for stems
            [lib,libtype,defines] = get_stem_lib(fliplr(csection(kss,2:3)),restraints.RNA);
            secdef_bck(ks,kss).lib = lib;
            secdef_bck(ks,kss).libtype = libtype;
            secdef_bck(ks,kss).defines = defines;
            ms = length(sections.fwd) + 1;
            defined_by = stem_defined(secdef_fwd,ms,0,sections_fwd);
            if ~isempty(defined_by)
                secdef_bck(ks,kss).defined_by = [0 defined_by];
            else
                defined_by = stem_defined(secdef_bck,ks,kss,sections_bck);
                if ~isempty(defined_by)
                    secdef_bck(ks,kss).defined_by = [1 defined_by];
                else
                    secdef_bck(ks,kss).defined_by = [];
                end
            end
        end
    end
end

% for the moment, only the first forward segment is generated

load nuclib_5

for ks = 1:length(sections_fwd)


    csection = sections_fwd{ks};
    [subsegments,~] = size(csection);
    
    % somewhat dirty workaround that catches a bug for trailing RNA loops
    if csection(end,end) - RNA.nta + 2 > length(RNA.sequence)
        break
    end

    modnum = 0;
    seg_trial = 0;
    tstart_seq = tic;
    runtime = toc(tstart_seq);
    
    while modnum < ensemble && runtime < 3600*maxtime 
        seg_trial = seg_trial + 1;
        maxsubtrials = 200*ones(1,subsegments);
        subtrials = zeros(1,subsegments);

        counters = ones(subsegments,2);
        scramble = cell(1,subsegments);
        scrambles = struct('previous',scramble,'next',scramble);
        timing = zeros(1,10000);
        levels = zeros(1,10000);
        timepoi = 0;
        atomtags = cell(1,10000);
        ecoor = zeros(10000,4);
        coorpoi = 0;
        dcoorpoi = zeros(1,subsegments);
        dseqpoi = zeros(1,subsegments);
        seqpoi = 0;
        sequence = char(1,1000);
        seqnum = zeros(1,1000);
        stem_def_kss = cell(1,subsegments); % to avoid overhead on repeated loading
        stems_kss = cell(1,subsegments);
        ecodes_kss = cell(1,subsegments); % to avoid overhead on recomputation
        termini = zeros(1,subsegments);
        for kss = 1:subsegments
            termini(kss) = secdef_fwd(ks,kss).nts(2);
            switch csection(kss,1)
                case 0
                    if secdef_fwd(ks,kss).nts(2) == secdef_fwd(ks,kss).nts(1)
                        maxsubtrials(kss) = 10; % single nt
                    else
                        maxsubtrials(kss) = 100; % loop
                    end
                case 1
                    maxsubtrials(kss) = 100; % stem
            end
        end
    %     profile off
    %     profile on
        options.save = false;
        kss = 0;
        while kss < subsegments && runtime < 3600*maxtime
            kss = kss + 1;
    %         if kss <= 2
    %             fprintf(1,'\nAt subsegment %i/%i\n',kss,subsegments);
    %         elseif kss == 3
    %             fprintf(1,'.');
    %         end
            success = false;
            secdef = secdef_fwd(ks,kss);
            if kss == 1 && ~isempty(secdef.anchori) % store initial anchor
                anchori = secdef.anchori;
                if kss == 1
                    initial_nt = secdef.nts(1);
                    initial_nt_O3p = anchori(6,:);
                end
            end
            % fprintf('Loop index %i\n',kss);
            switch csection(kss,1)
                case 0 % loop
                    if secdef.nts(2) == secdef.nts(1) % single nt
                        if subtrials(kss) > maxsubtrials(kss)
                            counters(kss,:) = [1,1];
                            subtrials(kss) = 0;
                            kss = kss - 2;
                            if kss < 0
                                ecoor = []; % invalidate result
                                break
                            else
                                coorpoi = coorpoi - dcoorpoi(kss+1);
                                seqpoi = seqpoi - dseqpoi(kss+1);
                                continue
                            end
                        end
                        %                 fprintf(1,'Insert single nt %i\n',secdef.nts(1));
                        if ~isempty(secdef.anchori) % initialized single nt
                            anchori = secdef.anchori;
                            if kss == 1
                                fixed = anchori(3,:);
                            end
                        else
                            anchori = get_terminal_anchor(ecoor(1:coorpoi,:),atomtags(1:coorpoi),termini(kss-1));
                        end
                        if ~isempty(secdef.anchore) % targeted single nt
                            anchore = secdef.anchore;
                        else
                            anchore = [];
                        end
                        base = RNA.sequence(secdef.nts(1)-RNA.nta+1);
                        seq3 = RNA.sequence(secdef.nts(1)-RNA.nta:secdef.nts(1)-RNA.nta+2);                    
                        counter = counters(kss,:);
                        err = 0;
                        while ~success && ~err
                            % for a targeted loop, correction is applied to the
                            % whole segment, this eases reaching the target anchor
                            options.clash_thr = clash_thr;
                            if ~isempty(anchore)
                                options.anchor0 = [initial_nt - secdef.nts(1) initial_nt_O3p];
                            elseif isfield(options,'anchor0')
                                options = rmfield(options,'anchor0');
                            end
                            [cecoor,catomtags,counter,~,err] = get_link_nt(fragments,non_clash_table,...
                                anchori,anchore,seq3,secdef.nts(1)-1,[],counter,options);
                            if ~isempty(cecoor) && isfield(options,'anchor0')
                                ecoor0 = [ecoor(1:coorpoi,:); cecoor];
                                fixed = [options.anchor0(2:4); anchore(1,:)];
                                ecoor0 = correct_section(ecoor0,correction,environ,fixed,options);
                                if ~isempty(ecoor0)
                                    ecoor(1:coorpoi,:) = ecoor0(1:coorpoi,:);
                                    cecoor = ecoor0(coorpoi+1:end,:);
                                else
                                    cecoor = [];
                                    err = -5;
                                end
                            end
                            if err == 0
                                cecoor = cecoor(2:end-2,:); % remove C4'(i-1) and P(i+1), C4'(i+1)
                                [m,~] = size(cecoor);
                                ecoor(coorpoi+1:coorpoi+m,:) = cecoor;
                                atomtags(coorpoi+1:coorpoi+m) = catomtags(2:end-2);
                                options.clash_thr = clash_thr;
                                options.max_shift = max_nt_shift;
                                [coor,min_dist,~,max_shift] = clash_repair(ecoor(1:coorpoi+m,2:4),environ,fixed,options);
                                subtrials(kss) = subtrials(kss) + 1;
                                if min_dist >= clash_thr && max_shift < max_nt_shift
                                    ecoor(1:coorpoi+m,2:4) = coor;
                                    dcoorpoi(kss) = m;
                                    sequence(seqpoi+1) = base;
                                    seqnum(seqpoi+1) = secdef.nts(1);
                                    dseqpoi(kss) = 1;
                                    seqpoi = seqpoi+1;
                                    coorpoi = coorpoi + m;
                                    counters(kss,:) = counter;
                                    %                             fname = 'test_RNA_2';
                                    %                             fname = wr_pdb_RNA(fname,ecoor(1:coorpoi,:),atomtags(1:coorpoi),sequence(1:seqpoi),seqnum(1:seqpoi));
                                    %                             fprintf(1,'Test RNA written to %s\n',fname);
                                    success = true;
                                end
                            end
                        end
                        if ~success
                            counters(kss,:) = [1,1];
                            subtrials(kss) = 0;
                            kss = kss - 2;
                            if kss < 0
                                ecoor = []; % invalidate result
                                break
                            else
                                coorpoi = coorpoi - dcoorpoi(kss+1);
                                seqpoi = seqpoi - dseqpoi(kss+1);
                                continue
                            end
                        end
                    else % loop with more than one nt
                        if subtrials(kss) > maxsubtrials(kss)
                            counters(kss,:) = [1,1];
                            subtrials(kss) = 0;
                            kss = kss - 2;
                            if kss < 0
                                ecoor = []; % invalidate result
                                break
                            else
                                coorpoi = coorpoi - dcoorpoi(kss+1);
                                seqpoi = seqpoi - dseqpoi(kss+1);
                                continue
                            end
                        end
                        if ~isempty(secdef.anchori) % initialized loop
                            anchori = secdef.anchori;
                            if kss == 1
                                fixed = anchori(3,:);
                            end
                        else
                            anchori = get_terminal_anchor(ecoor(1:coorpoi,:),atomtags(1:coorpoi),termini(kss-1));
                        end
                        if ~isempty(secdef.anchore) % targeted loop
                            anchore = secdef.anchore;
                            if isempty(ecodes_kss{kss})
                                ecodes = get_anchor_fragments(anchore,'back',shortfrag);
                                ecodes_kss{kss} = ecodes;
                            else
                                ecodes = ecodes_kss{kss};
                            end
                        else
                            anchore = [];
                            ecodes = [];
                        end
                        % sequence includes initial anchor and next nt
                        options.maxtime = [];
                        options.clash_thr = clash_thr;
                        options.max_shift = max_shift_loop;
                        seq = RNA.sequence(secdef.nts(1)-RNA.nta:secdef.nts(2)-RNA.nta+2);
                        [acodes,transmat] = get_RNA_initial_anchor(fragments,anchori(1:3,:),upper(seq(1)));
                        % for a targeted loop, correction is applied to the
                        % whole segment, this eases reaching the target anchor
                        if ~isempty(anchore)
                            options.anchor0 = [initial_nt - secdef.nts(1) initial_nt_O3p];
                        elseif isfield(options,'anchor0')
                            options = rmfield(options,'anchor0');
                        end
                        subtrials(kss) = subtrials(kss) + 1;
                        fprintf(1,'Insert loop nt %i-%i (Trial %i, Subseqment %i, Subtrial : %i)\n',secdef.nts,seg_trial,kss,subtrials(kss));
                        [cecoor,catomtags,seq,~,correction,err] = mk_RNA_loop(seq,fragments,...
                            shortfrag,non_clash_table,anchore,acodes,transmat,ecodes,options,secdef.ntoffset,[environ; ecoor(1:coorpoi,2:4)]);
                        [m,~] = size(cecoor);
                        if ~isempty(cecoor) && isfield(options,'anchor0')
                            ecoor0 = [ecoor(1:coorpoi,:); cecoor];
                            fixed = [options.anchor0(2:4); anchore(1,:)];
                            ecoor0 = correct_section(ecoor0,correction,environ,fixed,options);
                            if ~isempty(ecoor0)
                                ecoor(1:coorpoi,:) = ecoor0(1:coorpoi,:);
                                cecoor = ecoor0(coorpoi+1:end,:);
                            else
                                cecoor = [];
                            end
                        end
                        if ~isempty(cecoor)
                            success = true;
                            ecoor(coorpoi+1:coorpoi+m,:) = cecoor;
                            atomtags(coorpoi+1:coorpoi+m) = catomtags;
                            dcoorpoi(kss) = m;
                            n = length(seq);
                            sequence(seqpoi+1:seqpoi+n) = seq;
                            seqnum(seqpoi+1:seqpoi+n) = secdef.nts(1):secdef.nts(2);
                            dseqpoi(kss) = n;
                            seqpoi = seqpoi+n;
                            coorpoi = coorpoi + m;
                            %                     fname = 'test_RNA_1';
                            %                     fname = wr_pdb_RNA(fname,ecoor(1:coorpoi,:),atomtags(1:coorpoi),sequence(1:seqpoi),seqnum(1:seqpoi));
                            %                     fprintf(1,'Test RNA written to %s\n',fname);
%                         else % reinitialize subsegment computation and go back one subsegment
%                             %                     if err > -3
%                             %                         if secdef.nts(1) == 42
%                             %                             fname = 'test_RNA_2';
%                             %                             fname = wr_pdb_RNA(fname,ecoor(1:coorpoi,:),atomtags(1:coorpoi),sequence(1:seqpoi),seqnum(1:seqpoi));
%                             %                             fprintf(1,'Test RNA written to %s\n',fname);
%                             %                         end
%                             %                     end
%                             counters(kss,:) = [1,1];
%                             kss = kss - 2;
%                             if kss < 0
%                                 ecoor = []; % invalidate result
%                                 break
%                             else
%                                 coorpoi = coorpoi - dcoorpoi(kss+1);
%                                 seqpoi = seqpoi - dseqpoi(kss+1);
%                                 continue
%                             end
                        end
                    end
                case 1 % stem
                    if isempty(secdef.defined_by)
                        if subtrials(kss) > maxsubtrials(kss)
                            subtrials(kss) = 0;
                            counters(kss,:) = [1,1];
                            scramble = scramble_RNA_stems(stems_kss{kss});
                            scrambles(kss).previous = scramble.previous;
                            scrambles(kss).next = scramble.next;
                            kss = kss - 2;
                            if kss < 0
                                ecoor = []; % invalidate result
                                break
                            else
                                coorpoi = coorpoi - dcoorpoi(kss+1);
                                seqpoi = seqpoi - dseqpoi(kss+1);
                                continue
                            end
                        end
                        %                 fprintf(1,'Insert stem nt %i-%i\n',secdef.nts);
                        if ~isempty(secdef.anchori) % initialized stem
                            anchori = secdef.anchori;
                            if kss == 1
                                fixed = anchori(3,:);
                            end
                        else
                            anchori = get_terminal_anchor(ecoor(1:coorpoi,:),atomtags(1:coorpoi),termini(kss-1));
                        end
                        if ~isempty(secdef.anchore) % targeted stem
                            anchore = secdef.anchore;
                            ecodes = get_anchor_fragments(anchore,'back',shortfrag);
                        else
                            anchore = [];
                            ecodes = [];
                        end
                        % stem libraries should be loaded only on first encounter
                        if isempty(stem_def_kss{kss})
                            load(secdef.lib);
                            stems_kss{kss} = stems;
                            stem_def_kss{kss} = stem_def;
                        else
                            stem_def = stem_def_kss{kss};
                            stems = stems_kss{kss};
                        end
                        counter = counters(kss,:);
                        % initialize random selection of previous and next
                        % fragments, if necessary, otherwise retrieve random
                        % fragment order
                        if isempty(scrambles(kss).previous)
                            options.scramble = scramble_RNA_stems(stems_kss{kss});
                            scrambles(kss).previous = options.scramble.previous;
                            scrambles(kss).next = options.scramble.next;
                        else
                            options.scramble.previous = scrambles(kss).previous;
                            options.scramble.next = scrambles(kss).next;
                        end
                        % attach RNA stem runs only a single time in a given
                        % attempt, as it either returns a valid solution or
                        % exhausts the solution space
                        if secdef.defines(2) < secdef.nts(1)
                            ntnums = [secdef.defines(1):secdef.defines(2) secdef.nts(1):secdef.nts(2)];
                        else
                            ntnums = [secdef.nts(1):secdef.nts(2) secdef.defines(1):secdef.defines(2)];
                        end
                        [cecoor,catomtags,~,counter,err] = attach_RNA_stem(fragments,...
                            stem_def,stems,counter,anchori(1:3,:),anchore,ecodes,[environ; ecoor(1:coorpoi,2:4)],options,ntnums);
                        counters(kss,:) = counter;
                        subtrials(kss) = subtrials(kss) + 1;
                        if ~isempty(cecoor) % store subsegment information and retrieve current initial anchor
                            success = true;
                            [m,~] = size(cecoor);
                            ecoor(coorpoi+1:coorpoi+m,:) = cecoor;
                            atomtags(coorpoi+1:coorpoi+m) = catomtags;
                            dcoorpoi(kss) = m;
                            coorpoi = coorpoi + m;
                            n = length(stem_def.seq);
                            sequence(seqpoi+1:seqpoi+n) = stem_def.seq;
                            seqnum(seqpoi+1:seqpoi+n) = ntnums;
                            dseqpoi(kss) = n;
                            seqpoi = seqpoi+n;
                        else % reinitialize subsegment computation and go back one subsegment
                            counters(kss,:) = [1,1];
                            scramble = scramble_RNA_stems(stems_kss{kss});
                            scrambles(kss).previous = scramble.previous;
                            scrambles(kss).next = scramble.next;
                            kss = kss - 2;
                            if kss < 0
                                ecoor = []; % invalidate result
                                break
                            else
                                coorpoi = coorpoi - dcoorpoi(kss+1);
                                seqpoi = seqpoi - dseqpoi(kss+1);
                                continue
                            end
                        end
                    else
                        if secdef.defined_by(1) == 0
                            deftype = 'fwd';
                        else
                            deftype = 'bck';
                        end
                        fprintf(1,'Stem nts %i-%i were defined by %s section %i, subsection %i\n',secdef.nts,deftype,secdef.defined_by(2:3));
                        error('Predefined stems not yet implemented');
                        success = true;
                    end
                otherwise
                    error('Unknown case (neither stem nor loop)');
            end
            if kss < 4
                timepoi = timepoi+1;
                timing(timepoi) = toc(tstart_seq);
                levels(timepoi) = kss;
            end
            runtime = toc(tstart_seq);
        end

        if success && kss == subsegments % ### one should also test for completion of whole RNA
            modnum = modnum + 1;
            ecoor = ecoor(1:coorpoi,:);
            atomtags = atomtags(1:coorpoi);
            sequence = sequence(1:seqpoi);
            seqnum = seqnum(1:seqpoi);
%             timing = timing(1:timepoi);
%             levels = levels(1:timepoi);
            fname = sprintf('RNA_sec%i_mod%i',ks,modnum);
            fname1 = wr_pdb_RNA(fname,ecoor,atomtags,sequence,seqnum);
            fprintf(1,'RNA written to %s\n',fname1);
            save(fname,'ecoor','atomtags','sequence','seqnum','timing','levels');
            if modnum == 1,
                [~,snum] = add_pdb(fname);
                template_indices = [snum 1 1];
            else
                copy_structure(template_indices(1),'+mod',[],modnum,template_indices(1));
                structure = rd_pdb(pmodel);
                model_indices = template_indices;
                model_indices(3) = modnum;
                replace_model(model_indices,structure);
            end;
        end

        % profile viewer
    end
end

function [restrain,restraints,monitor,cancelled,number,number_monitor,msg] = process_RNA_restraints(restraints)
% Processes restraints for a single flexible domain

global rotamer_libraries

msg.error = 0;
msg.text = 'Restraints processed';

cancelled=false;

restrain = [];
monitor = [];
number = 0;
number_monitor = 0;

if ~isfield(restraints,'sequence'),
    msg.text = 'RNA sequence specification is missing.'; 
    restrain=[];
    msg.error = 11;
    return;
end;

for k = 1:length(restraints.sequence),
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
end;

monitor = restrain;

if isfield(restraints,'DEER'),
    poi = 0;
    llist = cell(0);
    label = cell(0);
    for k = 1:length(restraints.DEER),
        [indices,message] = resolve_address(restraints.DEER(k).adr1);
        if message.error ~= 2 && message.error ~= 13, % this site exists and is thus a beacon
            if message.error,
                address = mk_address(indices,1);
                msg.text = sprintf('ERROR: Residue address %s is invalid for %s.',restraints.DEER(k).adr1,address); 
                msg.error = 12;
                restrain=[];
                return;
            end;
            poi = poi + 1;
            restraints.DEER(k).type1 = 1;
            indices(3) = modnum;
            restraints.DEER(k).indices1 = indices;
            llist{poi} = mk_address(indices);
            label{poi} = restraints.DEER(k).label1;
        else % this site is in the RNA to be modelled
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.DEER(k).label1) || strcmpi(rotamer_libraries(kr).tc,restraints.DEER(k).label1)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end;
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel1 = NO;
                    restraints.DEER(k).type1 = 0;
                end;
            end;
        end;
        [indices,message]=resolve_address(restraints.DEER(k).adr2);
        if message.error ~= 2 && message.error ~= 13, % this site exists and is thus a beacon
            if message.error,
                address = mk_address(indices,1);
                msg.text = sprintf('ERROR: Residue address %s is invalid for %s.',restraints.DEER(k).adr2,address); 
                msg.error = 13;
                restrain=[];
                return;
            end;
            poi = poi + 1;
            restraints.DEER(k).type2 = 1;
            indices(3) = modnum;
            restraints.DEER(k).indices2 = indices;
            llist{poi} = mk_address(indices);
            label{poi} = restraints.DEER(k).label2;
        else % this site is in the RNA to be modelled
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.DEER(k).label2) || strcmpi(rotamer_libraries(kr).tc,restraints.DEER(k).label2)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end;
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel2 = NO;
                    restraints.DEER(k).type2 = 0;
                end;
            end;
        end;
    end;
    labels = get_labels(llist,label);
    for k = 1:length(restraints.DEER),
        clabel1 = restraints.DEER(k).label1;
        clabel2 = restraints.DEER(k).label2;
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 0, % internal restraint
            resa = correct_section_address(restraints.DEER(k).adr1,restraints.nta);
            resb = correct_section_address(restraints.DEER(k).adr2,restraints.nta);
            if isempty(resa) || isempty(resb),
                msg.text = sprintf('ERROR: Residue address %s or %s is invalid.',restraints.DEER(k).adr1,restraints.DEER(k).adr2); 
                msg.error = 14;
                restrain=[];
                return;
            end;
            NO1 = restraints.DEER(k).NO_rel1;
            NO2 = restraints.DEER(k).NO_rel2;
            if resa < resb,
                exch = resa; resa = resb; resb = exch;
                exch = NO1; NO1 = NO2; NO2 = exch;
            end;
            if restraints.DEER(k).r ~= 0,
                [restrain,number] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,clabel1,clabel2);
            else
                [monitor,number_monitor] = mk_internal_restraint(monitor,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,clabel1,clabel2);
                number_monitor = number_monitor + 1;
            end;
        end;
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 ==0, % beacon restraint, first residue beacon
            indices = restraints.DEER(k).indices1;
            adr1 = mk_address(indices);
            for kl = 1:length(labels),
                if sum(abs(indices-labels(kl).indices)) == 0,
                    xyz_beacon = labels(kl).xyz;
                end;
            end;
            res_loop = correct_section_address(restraints.DEER(k).adr2,restraints.nta);
            if isempty(res_loop),
                msg.text = sprintf('ERROR: Residue address %s inside domain is invalid.',restraints.DEER(k).adr2); 
                msg.error = 15;
                restrain=[];
                return;
            end;
            NO = restraints.DEER(k).NO_rel2;
            if restraints.DEER(k).r ~= 0,
                [restrain,number] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,clabel1,clabel2,indices,adr1);
            else
                [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,clabel1,clabel2,indices,adr1);
                number_monitor = number_monitor + 1;
            end;
        end;
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 1, % beacon restraint, second residue beacon
            indices = restraints.DEER(k).indices2;
            adr2 = mk_address(indices);
            for kl = 1:length(labels),
                if sum(abs(indices-labels(kl).indices)) == 0,
                    xyz_beacon = labels(kl).xyz;
                end;
            end;
            res_loop = correct_section_address(restraints.DEER(k).adr1,restraints.nta);
            if isempty(res_loop),
                msg.text = sprintf('ERROR: Residue address %s inside domain is invalid.',restraints.DEER(k).adr1); 
                msg.error = 16;
                restrain=[];
                return;
            end;
            NO = restraints.DEER(k).NO_rel1;
            if restraints.DEER(k).r ~= 0,
               [restrain,number] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,clabel1,clabel2,indices,adr2);
            else
               [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,clabel1,clabel2,indices,adr2);
                number_monitor = number_monitor + 1;
            end;
        end;
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 == 1, % nonsensical restraint inside defined structure
            add_msg_board(sprintf('Warning: Restraint between residues %s and %s inside rigid bodies will be ignored.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
        end;
    end;
end;

if isfield(restraints,'oligomer'),
    for k = 1:length(restraints.oligomer),
        clabel = restraints.oligomer(k).label;
        if strcmpi(restraints.oligomer(k).label,'C5'),
            NO = [];
        else
            for kr = 1:length(rotamer_libraries),
                if strcmpi(rotamer_libraries(kr).label,restraints.oligomer(k).label),
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps,
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end;
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end;
            end;
        end;
        res = restraints.oligomer(k).num;
        if isempty(res) || isnan(res),
            msg.text = sprintf('ERROR: Residue number %i is invalid in oligomer restraint.',res); 
            msg.error = 17;
            restrain=[];
            return;
        end;
        if restraints.oligomer(k).r ~= 0,
            [restrain,number] = mk_oligomer_restraint(restrain,NO,res,restraints.oligomer(k).mult,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number,clabel);
        else
            [monitor,number_monitor] = mk_oligomer_restraint(monitor,NO,res,restraints.oligomer(k).mult,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number_monitor,clabel);
            number_monitor = number_monitor + 1;
        end;
    end;
end;

if isfield(restraints,'stems'),
    for k = 1:length(restraints.stems),
        for kn = restraints.stems(k).C5pf:restraints.stems(k).C3pf
            restrain(kn - restraints.nta + 1).stem = 1;
            restrain(kn - restraints.nta + 1).stem_id = k;
        end
        for kn = restraints.stems(k).C5pb:restraints.stems(k).C3pb
            restrain(kn - restraints.nta + 1).stem = 1;
            restrain(kn - restraints.nta + 1).stem_id = k;
        end
    end;
end;

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
    end;
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
for k = 1:length(rot_lib.library),
    coor = rot_lib.library(k).ecoor;
    NO = NO + pops(k)*(coor(midNO(1),2:4) + coor(midNO(2),2:4))/2;
end;
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
if rmean > 0 && sigr > 0,
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
