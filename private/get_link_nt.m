function [ecoor,atomtags,counter,code,err] = get_link_nt(fragments,nct,anchori,anchore,seq,ntoffset,environ,counter,options)
% function [ecoor,atomtags,counter,err] = get_link_nt(fragments,nct,anchori,anchore,seq,ntoffset,environ,counter,options)
% 
% tries to find nucleotide (RNA fragment) that links an initial anchor
% anchori and a target anchor anchore, to be used instead of mk_RNA_loop if
% a link is only one nucleotide long, at least one anchor must be present
%
% fragments     fragment library, as generated by recompile_library.m
% nct           talbe of non-clashing fragment pairs
% anchori       coordinates of the C4'(i-2), P(i-1), C4'(i-1), and C3'(i-1) 
%               atoms of the initial anchor nucleotide, defaults to empty 
%               array, which is interpreted as absence of an initial anchor
% anchore       coordinates of the P(i+1),C4'(i+1),P(i+2) atoms of the
%               final anchor nucleotide, defaults to empty array, which is
%               interpreted as absence of a final anchor
% seq           single letter codes from previous up to next nt, 
%               allowed nucleotides A,C,G,U
% ntoffset      nucleotide index in RNA of the preceding nucleotide,
%               defaults to 0
% environ       (ne,3) array of Cartesian coordinates of environment atoms,
%               if present and not empty, the model is tested for clashes 
%               with the environment
% counter       for an nt anchored only at one side, counter is an index
%               into the fragment library, the first fitting fragment with
%               this or a larger index is returned, useful for testing all
%               allowed fragments in an external loop, defaults to [1,1]
% options       computation options, struct with fields
%               .max_rot    maximum rotation of the loop end for matching 
%                           final anchor fragment, defaults to 20� per 4 nt
%               .anchor_acc RMSD accuracy of anchor coordinate fit,
%                           defaults to 0.5 �
%               .fit_tol    contribution of a nucleotide to the convergence
%                           radius for a final anchor P atom, defaults to
%                           0.2 �
%               .anchor0    vector (1,4) defining an initial anchor of a 
%                           longer RNA segment by its nucleotide number and
%                           P atom coordinates; if present and if the loop
%                           is an insert, deformation of the whole segment
%                           is allowed for reaching the target anchor, the
%                           nucleotide number is relative to the first nt
%                           of the loop
%
% ecoor         extended coordinates, first column: nucleotide index,
%               columns 2-4: coordinates, empty if function fails
% atomtags      cell array of atom tags corresponding to ecoor
% counter       update counter
% code          fragment codes
% err           error code, 0 if successful, else
%               -1  neither initial nor target anchor provided
%               -2  no valid single-letter code provided
%               -3  distance between initial and target anchor too long
%               -4  no fragment in the library fits acceptably well
%               -5  a fitting fragment was found, but it clashes with the
%                   environment
%               -6  no pair of inital and target fragments leads to
%                   sufficient C4' match of the nucleotide
%               -7  run out of initial anchor fragments to test
%               -8  run out of target anchor fragments to test
%
% G. Jeschke, 31.12.2017

ecoor = [];
atomtags = cell(0,1);
code = [];
err = 0;

anchor_acc = 0.5;
clash_thr = 1.5;
rPO = 1.61; % typical P-O bond length
maxdrPO = 0.2; % maximum allowed deviation of PO bond length

maxlen_per_nt = 7;
% fprintf(1,'Counter at (%i,%i)\n',counter);

if ~exist('counter','var') || isempty(counter) 
    counter = [1,1];
end

if ~exist('environ','var')
    environ = [];
end

if ~exist('ntoffset','var') || isempty(ntoffset) 
    ntoffset = 0;
end

if ~exist('seq','var') || isempty(seq)
    err = -2;
    return
end

if ~exist('anchore','var') 
    anchore = [];
end

if ~exist('anchori','var') 
    anchori = [];
end

if ~exist('options','var') || ~isfield(options,'max_rot')
    options.max_rot = 20; % maximum rotation of the loop end for matching final anchor fragment
end

if ~exist('options','var') || ~isfield(options,'anchor_acc')
    options.anchor_acc = 0.5; % RMSD accuracy of anchor coordinate fit
end

if ~exist('options','var') || ~isfield(options,'fit_tol')
    options.fit_tol = 0.2; % tolerance per nucleotide for convergence radius of target P atom
end

if isfield(options,'anchor0')
    len = length(seq) - options.anchor0(1) - 2;
%     disp('Aber hallo');
else
    len = length(seq);
end
fit_thresh = options.fit_tol*len;
max_rot = len*options.max_rot/4;


if isempty(anchore)
    targeted = false;
else
    targeted = true;
    Ptarget = anchore(1,:);
end

if isempty(anchori)
    initialized = false;
else
    initialized = true;
    O3pi = anchori(4,:);
    anchori = anchori(1:3,:);
    Pinitial = O3pi;
end

if ~targeted && ~initialized
    err = -1;
    return
end

if targeted && initialized
    P2dist = norm(anchori(2,:)-anchore(1,:));
    if P2dist > 2*maxlen_per_nt + fit_thresh
        err = -3;
        if P2dist < 5 
            fprintf(2,'nt not quite in reach of target anchor (%5.2f �)\n',P2dist);
        end
        return
    end
    fprintf(1,'Single nt within reach of target anchor (P2dist = %5.2f �)\n',P2dist);
end


% determine all allowed initial fragments
if initialized
    iallowed = zeros(1,length(fragments));
    iC4p = zeros(length(fragments),3);
    iP = zeros(length(fragments),3);
    ipoi = 0;
    for kf = 1:length(fragments)
        sfrag = fragments(kf).A.coor(fragments(kf).A.assign.previous(1:3),:);
        [rmsd, ~, transmat] = superimpose_3points(anchori,sfrag);
        if rmsd <= anchor_acc
            ipoi = ipoi + 1;
            iallowed(ipoi) = kf;
            sfrag = fragments(kf).A.coor(fragments(kf).A.assign.next(2:3),:);
            sfrag = [sfrag ones(2,1)]*transmat';
            iP(ipoi,:) = sfrag(1,1:3);
            iC4p(ipoi,:) = sfrag(2,1:3);
        end
    end
    iallowed = iallowed(1:ipoi);
    iP = iP(1:ipoi,:);
    iC4p = iC4p(1:ipoi,:);
    if counter(1) > ipoi
        err = -7;
        return
    end
end

% determine all allowed target fragments
if targeted
    eallowed = zeros(1,length(fragments));
    eC4p = zeros(length(fragments),3);
    epoi = 0;
    for kf = 1:length(fragments)
        sfrag = fragments(kf).A.coor(fragments(kf).A.assign.next,:);
        [rmsd, ~, transmat] = superimpose_3points(anchore,sfrag);
        if rmsd <= anchor_acc
            epoi = epoi + 1;
            eallowed(epoi) = kf;
            sfrag = fragments(kf).A.coor(fragments(kf).A.assign.previous(1),:);
            sfrag = [sfrag 1]*transmat';
            eC4p(epoi,:) = sfrag(1:3);
        end
    end
    eallowed = eallowed(1:epoi);
    eC4p = eC4p(1:epoi,:);
    if counter(1) > epoi
        err = -8;
        return
    end
end


coor = [];
if initialized
    if ~targeted
        success = false;
        for kf = counter(1):ipoi
            to_anchor = [anchori(3,:); iP(kf,:); iC4p(kf,:)];
            ifcode = iallowed(kf);
            for ks = counter(2):length(fragments)
                sfrag = fragments(ks).(base).coor(fragments(ks).(base).assign.previous(1:3),:);
                [rmsd, ~, transmat] = superimpose_3points(to_anchor,sfrag);
                if rmsd <= anchor_acc
                    coor0 = fragments(ifcode).(base).coor;
                    [m,~] = size(coor0);
                    coor0 = [coor0 ones(m,1)];
                    coor0 = coor0*transmat';
                    coor0 = coor0(:,1:3);
                    atomtags0 = fragments(ifcode).(base).atomtags;
                    rPOi = norm(coor0(fragments(ks).(base).assign.previous(2),:)-O3pi);
                    drPOi = abs(rPOi-rPO);
                    if drPOi < maxdrPO
                        coor = coor0;
                        atomtags = atomtags0;
                        code = [ifcode,ks];
                        success = true;
                        if ks < length(fragments)
                            counter = [kf ks+1];
                        else
                            counter = [kf+1 1];
                        end
                        break
                    end
                end
            end
            if success
                break
            end
        end
    end
else % only targeted
    success = false;
    for kf = counter(1):epoi
        to_anchor = [eC4p(kf,:); anchore(1:2,:)];
        ifcode = eallowed(kf);
        for ks = counter(2):length(fragments)
            sfrag = fragments(ks).(base).coor(fragments(ks).(base).assign.next,:);
            [rmsd, sfrag, transmat] = superimpose_3points(to_anchor,sfrag);
            if rmsd <= anchor_acc
                coor0 = fragments(ks).(base).coor;
                [m,~] = size(coor0);
                coor0 = [coor0 ones(m,1)];
                coor0 = coor0*transmat';
                coor0 = coor0(:,1:3);
                atomtags0 = fragments(ks).(base).atomtags;
                % find O3' atom of current fragment
                poiO3p = 0;
                for k = 1:length(atomtags0)
                    if strcmp(atomtags0{k},'O3''')
                        poiO3p = k;
                    end
                end
                rPOe = norm(sfrag(2,:)-coor0(poiO3p,:));
                drPOe = abs(rPOe-rPO);
                if drPOe <= maxdrPO % P-O bond with target nt has an allowed length
                    coor = coor0;
                    atomtags = atomtags0;
                    code = [ks ifcode];
                    success = true;
                    if ks < length(fragments)
                        counter = [kf ks+1];
                    else
                        counter = [kf+1 1];
                    end
                    break
                end
            end
        end
        if success
            break
        end
    end
end
    
if initialized && targeted
    min_rot = 1e6;
    min_approach = 1e6;
    % test all possible fragments
    pcode = strfind(slc,seq(1));
    base = seq(2);
    mycode = strfind(slc,base);
    ncode = strfind(slc,seq(3));
    for ki = 1:length(iallowed)
        pair = 4*(pcode-1) + mycode;
        non_clash_table = nct{pair};
        available = non_clash_table{iallowed(ki)};
        to_anchor = [anchori(3,:); iP(ki,:); iC4p(ki,:)];
        % try all available fragments
        for ka = 1:length(available)
            ks = available(ka);
            sfrag = fragments(ks).(base).coor(fragments(ks).(base).assign.previous(1:3),:);
            [rmsd, sfrag, transmat] = superimpose_3points(to_anchor,sfrag);
            if rmsd <= anchor_acc
                Pcurr = sfrag(2,:);
                if norm(Pcurr - Ptarget) <= maxlen_per_nt + fit_thresh
                    coor0 = fragments(ks).(base).coor;
                    [m,~] = size(coor0);
                    coor0 = [coor0 ones(m,1)];
                    coor0 = coor0*transmat';
                    coor0 = coor0(:,1:3);
                    coort = coor0(fragments(ks).(base).assign.next,:);
                    atomtags0 = fragments(ifcode).(base).atomtags;
                    pair = 4*(mycode-1) + ncode;
                    non_clash_table = nct{pair};
                    ecodes = non_clash_table{available(ka)}; % non-clashing target fragments
                    epoi = 0;
                    % check which non-clashing fragments are allowed as target
                    % fragments
                    for ke = 1:length(ecodes)
                        if min(abs(ecodes-eallowed)) < eps
                            epoi = epoi + 1;
                            ecodes(epoi) = ecodes(ke);
                        end
                    end
                    ecodes = ecodes(1:epoi);
                    for ke = 1:length(ecodes)
                        kf = ecodes(ke);
                        fanchor = fragments(kf).A.coor([fragments(kf).A.assign.previous(2:3) fragments(kf).A.assign.next(2)],:);
                        [~, ~, transmat0] = superimpose_3points(anchore,fanchor);
                        linkage = fragments(kf).A.coor(fragments(kf).A.assign.previous,:);
                        linkage = [linkage ones(3,1)]*transmat0';
                        [rmsd, ~, transmat_link] = superimpose_3points(linkage(:,1:3),coort);
                        EV = affine2EV(transmat_link);
                        approach = norm(transmat_link(1:3,4));
                        if approach < min_approach
                            min_approach = approach;
                        end
                        if abs(EV(4)) < min_rot && rmsd < 2*options.anchor_acc && approach < fit_thresh
                            err = 0;
                            coor = coor0;
                            atomtags = atomtags0;
                            correction.trans = transmat_link(1:3,4)';
                            correction.rot = EV;
                            min_rot = abs(EV(4));
                        end
                    end
                end
            end
        end
    end
end
    
if isempty(coor)
    err = -4;
    return
end

ecoor = [(ntoffset+1)*ones(m,1) coor];

if ~isempty(environ) && ~isfield(options,'anchor0')
    coor = ecoor(:,2:4);
    pair_dist = get_all_pair_dist(coor,environ);
    [mdv,mdp] = min(pair_dist);
    [min_dist,mdp2] = min(mdv);
    mdp1 = mdp(mdp2);
    if min_dist < clash_thr
        fprintf(1,'RNA nt atom %s of nt %i clashed with environment atom %i at distance %4.2f �.\n',...
            atomtags{mdp1},ecoor(mdp1,1),mdp2,min_dist);
        err = -5;
        ecoor = [];
    end
end

