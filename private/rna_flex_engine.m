function diagnostics = rna_flex_engine(handles,secdef,RNA,environ,ensemble,maxtime,id)
% diagnostics = rna_flex_engine(handles,secdef,RNA,environ,ensemble,maxtime)
%
% handles       handles to GUI figure
% secdef        section definition, structure with fields
%               .nts        first and last nucleotide of the segment
%               .anchori    [6,3] array of initial anchor coordinates
%               .anchore    [3,3] array of target anchor coordinates
%               .ntoffset   offset for nucleotide numbers
% RNA           definition of the RNA, structure with fields
%               .sequence   sequence string
%               .nta        number of first nucleotide
% environ       coordinate array of environment atoms
% ensemble      number of models in the ensemble, optional, defaults to 1
% maxtime       maximum runtime per model, optional, defaults to 12 h
% id            identifier vector [1,2] for rigid-body arrangement and
%               section, defaults to empty vector, if empty, a generic
%               file name is used for the saved PDB file
%
% G. Jeschke, 02.05.2018

diagnostics.snum = [];
diagnostics.success = 0;
diagnostics.time_per_model = 0;

snum = [];

clash_thr = 1.0;
max_nt_shift = 2;
max_shift_loop = 2.5;

if ~exist('ensemble','var') || isempty(ensemble)
    ensemble = 1;
end

if ~exist('maxtime','var') || isempty(maxtime)
    maxtime = 12;
end

if ~exist('id','var')
    id = [];
end


load nuclib_5

modnum = 0;
tstart_seq = tic;
runtime = toc(tstart_seq);
failure = false;

while modnum < ensemble && runtime < 3600*maxtime && ~failure
    seg_trial = 0;
    
    counters = ones(1,2);
    atomtags = cell(1,10000);
    ecoor = zeros(10000,4);
    coorpoi = 0;
    seqpoi = 0;
    sequence = char(1,1000);
    seqnum = zeros(1,1000);
    termini = 0;
    options.save = false;
    kss = 0;
    success = false;
    while runtime < 3600*maxtime && ~success && ~failure
        seg_trial = seg_trial + 1;
        if ~isempty(secdef.anchori) % store initial anchor
            anchori = secdef.anchori;
            initial_nt = secdef.nts(1);
            initial_nt_O3p = anchori(6,:);
        end
        if secdef.nts(2) == secdef.nts(1) % single nt
            if ~isempty(secdef.anchori) % initialized single nt
                fixed = anchori(3,:);
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
                    if min_dist >= clash_thr && max_shift < max_nt_shift
                        ecoor(1:coorpoi+m,2:4) = coor;
                        sequence(seqpoi+1) = base;
                        seqnum(seqpoi+1) = secdef.nts(1);
                        seqpoi = seqpoi+1;
                        coorpoi = coorpoi + m;
                        counters(kss,:) = counter;
                        success = true;
                    end
                end
            end
        else % loop with more than one nt
            if ~isempty(secdef.anchori) % initialized loop
                anchori = secdef.anchori;
                fixed = anchori(3,:);
            else
                anchori = get_terminal_anchor(ecoor(1:coorpoi,:),atomtags(1:coorpoi),termini(kss-1));
            end
            if ~isempty(secdef.anchore) % targeted loop
                anchore = secdef.anchore;
                ecodes = get_anchor_fragments(anchore,'back',shortfrag);
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
            fprintf(1,'Insert loop nt %i-%i (Trial %i)\n',secdef.nts,seg_trial);
            [cecoor,catomtags,seq,~,correction,err] = mk_RNA_loop(seq,fragments,...
                shortfrag,non_clash_table,anchore,acodes,transmat,ecodes,options,secdef.ntoffset,[environ; ecoor(1:coorpoi,2:4)]);
            if err == -3
                failure = true; % too large distance per nucleotide
            end
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
                n = length(seq);
                sequence(seqpoi+1:seqpoi+n) = seq;
                seqnum(seqpoi+1:seqpoi+n) = secdef.nts(1):secdef.nts(2);
                seqpoi = seqpoi+n;
                coorpoi = coorpoi + m;
            end
        end
        runtime = toc(tstart_seq);
        handles.text_auxiliary_fail.String = sprintf('%5.2f',100*modnum/ensemble);
        handles.text_core_fail.String = sprintf('%5.2f',100*runtime/(3600*maxtime));
        drawnow
    end
    
    if success
        modnum = modnum + 1;
        ecoor = ecoor(1:coorpoi,:);
        atomtags = atomtags(1:coorpoi);
        sequence = sequence(1:seqpoi);
        seqnum = seqnum(1:seqpoi);
        fname = sprintf('RNA_mod%i',modnum);
        if ~isempty(id)
            fname = sprintf('RNA_rba%i_sec%i_mod%i',id(1),id(2),modnum);
        end
        fname1 = wr_pdb_RNA(fname,ecoor,atomtags,sequence,seqnum);
        fprintf(1,'RNA written to %s\n',fname1);
%         save(fname,'ecoor','atomtags','sequence','seqnum');
        if modnum == 1
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
end

diagnostics.snum = snum;
diagnostics.success = modnum;
if modnum > 0
    diagnostics.time_per_model = runtime/modnum;
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
