function [ecoor,atomtags,code,counter,err] = attach_RNA_stem(fragments,...
    stem_def,stems,counter,anchori,anchore,ecodes,environ,options,ntnums)
% [ecoor,atomtags,code,counter,err] = 
% attach_RNA_stem(stem_def,stems,counter,anchori,anchore,ecodes,ntnums)
%
% attaches an RNA stem to an (optional) initial anchor and targets it to an
% (optional) final anchor
%
% fragments fragment library, as generated by recompile_library.m
% stem_def  stem definitions for stem library
% stems     stem coordinates and additional information for each stem
%           library member
% counter   (1,2) array of counters, first element is index into stems,
%           second element index into preceding fragment vector, defaults
%           to [1,1], used for exhaustive sampling of the stem library
%           under external control
% anchori   coordinates of the C4'(i-1), P(i), and C4'(i) atoms of the
%           initial anchor nucleotide, defaults to empty array, which is
%           interpreted as absence of an initial anchor
% anchore   coordinates of the P(j),C4'(j),P(j+1) atoms of the
%           final anchor nucleotide, defaults to empty array, which is
%           interpreted as absence of a final anchor
% ecodes    allowed fragment codes at final anchor, optional, defaults to
%           1:10000 (any fragment in any reasonable fragment library), it
%           is recommended to provide ecodes if anchore is provided
% environ   (ne,3) array of Cartesian coordinates of environment atoms, if
%           present and not empty, the model is tested for clashes with the
%           environment
% options   computation options
%           .maxtime    maximum run time, defaults to 1 h, should always be
%                       faster than that
%           .scramble   scrambled index vectors for allowed previous and
%                       next fragments to avoid statistical bias, optional
%           .clash_thr  clash threshold, defaults to 1.4 �
%           .anchor_acc threshold for acceptable anchor fit, defaults to
%                       0.5 � 
%           .max_rot    maximum rotation of the loop end for matching 
%                       final anchor fragment, defaults to 20�
%           .fit_tol    contribution of a base pair to the convergence
%                       radius for a final anchor, defaults to
%                       0.2 �
%           .max_shift  maximum shift allowed for an atom in clash repair,
%                       defaults to 2 �
%           .save       solution is saved as a PDB file, optional flag,
%                       defaults to false
% ntnums    numbers of the nucleotides, defaults to original numbers
%
% ecoor     extended coordinates (n,4) first column: nucleotide index,
%           columns 2-4: atom coordinates
% atomtags  cell string of PDB atom tags (length poi+2)
% code      assigned fragment index for preceding and following fragment
%           NaN if no assignment is made
% err       error code
%           -9  undefined error
%           -8  all stem/anchor fragment combinations have been exhausted
%           -7  too large atom shift in clash repair
%           -6  no matching code with final anchor
%           -5  allowed final anchor fragment does not superimpose well 
%               with anchor  
%           -4  allowed initial anchor fragment does not superimpose well 
%               with anchor  
%           -3  target P atom can never be reached
%           -2  target anchor was not reached
%           -1  model clashed with environment
%
% G. Jeschke, 26.12.2017

% initialize output
ecoor = [];
code = nan(1,2);
err = 0;

nat = length(stem_def.atomtags);
atomtags = stem_def.atomtags;

if ~exist('options','var') || ~isfield(options,'maxtime')
    options.maxtime = 1;
end

if isfield(options,'scramble')
    scramble = options.scramble;
else
    scramble = scramble_RNA_stems(stems);
end

if ~isfield(options,'clash_thr') || isempty(options.clash_thr)
    options.clash_thr = 1.4;
end

if ~isfield(options,'anchor_acc')
    options.anchor_acc = 0.5; % RMSD accuracy of anchor coordinate fit
end

if ~isfield(options,'max_rot')
    options.max_rot = 20; % maximum rotation of the stem end for matching final anchor fragment
end

if ~isfield(options,'fit_tol')
    options.fit_tol = 0.2; % tolerance per nucleotide for convergence radius of target P atom
end

if ~isfield(options,'max_shift')
    options.max_shift = 2.0; % maximum shift allowed for an atom in clash repair
end

if ~isfield(options,'save')
    options.save = false; % by default, the solution is not automatically saved to a PDB file
end

if ~exist('environ','var') || isempty(environ)
    clash_test = false;
else
    clash_test = true;
end

if ~exist('ecodes','var') || isempty(ecodes)
    ecodes = 1:10000;
end

if ~exist('anchore','var')
    anchore = [];
end

fixed = [];

if ~isempty(anchore)
    Ptarget = anchore(1,:); % target coordinate of P atom
    targeted = true;
else
    targeted = false;
end

if ~exist('anchori','var')
    anchori = [];
end

if ~isempty(anchori)
    initialized = true;
else
    initialized = false;
end

ntnum = stems(1).ecoor(:,1);
if exist('ntnums','var') && ~isempty(ntnums)
    for k = 1:length(ntnum)
        ntnum(k) = ntnums(ntnum(k));
    end
end

pindices = scramble.previous{counter(1)};
nindices = scramble.next{counter(1)};

corr.trans = [0,0,0]; % by default, no correction translation
corr.rot = eye(3); % by default, no correction rotation
corr.length = 0;

fit_thresh = options.fit_tol*stem_def.length; 

success = false;
failed = false;
correction = false;

max_min_dist = 0;
while ~success && ~failed
    success = true;
    err = 0;
    transmat = eye(4);
    ptal = [stems(counter(1)).ptal ones(3,1)];
    cstem = [stems(counter(1)).ecoor(:,2:4) ones(nat,1)];
    if initialized % the stem is attached to the initial anchor
        % select a fragment and extract all pseudo-torsion defining
        % atom coordinates
        ifcode = pindices(counter(2));
        code(1) = ifcode;
        sfrag = [fragments(ifcode).A.coor(fragments(ifcode).A.assign.previous(1:3),:);...
            fragments(ifcode).A.coor(fragments(ifcode).A.assign.next(2:3),:);];
            % amend to coordinate array for affine transformation and transform to
        % extract the coordinates that must fit the initial anchor
        to_anchor = sfrag(1:3,:);
        % find the affine transformation that superimposes the fragment to
        % the anchor
        [rmsd, ~, transmat] = superimpose_3points(anchori,to_anchor);
        % try the next allowed anchor fragment or next stem, if superposition is poor 
        if rmsd > options.anchor_acc
            err = -4;
            [counter,pindices,nindices,failed] = count_up(counter,stems,scramble);
            success = false;
            continue
        else
            success = true;
        end
        % transform anchor fragment to anchor coordinates
        sfrag = [sfrag ones(5,1)];
        sfrag = transmat*sfrag';
        sfrag = sfrag';
        % extract the coordinates that must fit the stem
        steminit = sfrag(3:5,1:3);
        % extract initial stem coordinates
        steminit0 = [stems(counter(1)).prev_coor(ifcode,:); stems(counter(1)).ptaf(1:2,:)];
        [rmsd, ~, transmat] = superimpose_3points(steminit,steminit0);
        if rmsd > options.anchor_acc
            err = -4;
            [counter,pindices,nindices,failed] = count_up(counter,stems,scramble);
            success = false;
            continue
        else
            success = true;
        end
        % transform stem so that it ataches to anchor
        cstem = cstem*transmat';
        % store coordinate of initial P atom
        Pinitial = sfrag(4,1:3);
        fixed = Pinitial;
        corr.Pinitial = Pinitial;
        % transform pseudo-atom coordinates of last nt of both strands
        ptal = ptal*transmat';
    end    
    minrmsd = 1e6;
    minrmsd2 = 1e6;
    if success && targeted % the stem must attach to a final anchor
        if initialized  % a frame transformation was already performed,
                        % we can only match by a small translation/rotation
                        % this is untested
                        error('Stem insertion with initial and final anchor is untested.\nPlease contact gjeschke@ethz.ch');
            correction = true;
            ontarget = false;
            efn = 0; % counter for possible final anchor fragments
            min_rot = options.max_rot+eps;
            while efn < length(nindices)
                efn = efn + 1;
                efcode = nindices(efn); % fragment index of anchor fragment
                [match,mpoi] = min(abs(ecodes-efcode));
                while match ~= 0 && efn < length(nindices)
                    efn = efn +1;
                    efcode = nindices(efn);
                    [match,mpoi] = min(abs(ecodes-efcode));
                end
                if match ~= 0
                    err = -6;
                    break
                end
                code(2) = efcode;
                % determine code of anchor fragment
                anchorcode = ecodes(mpoi);
                % pseudo-torsion defining atom coordinates of final
                % anchor fragment C4'(j),P(j+1),C4'(j+1),P(j+2)
                efrag = [fragments(anchorcode).A.coor(fragments(anchorcode).A.assign.previous(1:3),:);...
                    fragments(anchorcode).A.coor(fragments(anchorcode).A.assign.next(2),:);];
                % extract the ones that must superimpose with anchor
                % fragment
                fanchor = efrag(2:4,:);
                % find the affine transformation that superimposes the fragment to
                % the anchor
                [rmsd, ~, ftransmat] = superimpose_3points(anchore,fanchor);
                % try the next allowed anchor fragment or next stem, if superposition is poor
                if rmsd > options.anchor_acc
                    continue
                end
                % transform anchor fragment to anchor coordinates
                efrag = [efrag ones(4,1)];
                efrag = ftransmat*efrag';
                efrag = efrag';
                % extract the coordinates that must fit the stem fragment,
                % these are C4'(j),P(j+1),C4'(j+1)
                stemfinal = efrag(1:3,1:3);
                % extract corresponding stem fragment coordinates
                % C4'(j),P(j+1),C4'(j+1)
                stemfinal0 = ptal(:,1:3);
                [rmsd, ~, transmat_link] = superimpose_3points(stemfinal,stemfinal0);
                EV = affine2EV(transmat_link);
                % check if this correction is allowed
                if abs(EV(4)) < min_rot && rmsd < 2*options.anchor_acc && norm(transmat_link(1:3,4)) < fit_thresh
                    ontarget = true;
                    corr.trans = transmat_link(1:3,4)';
                    corr.rot = EV;
                    min_rot = abs(EV(4));
                    Pfinal = stemfinal(2,:);
                    corr.length = norm(Pfinal-corr.Pinitial);
                    fixed = [corr.Pinitial; Pfinal];
                end
            end
            if ~ontarget
                success = false;
            else
                success = true;
            end
        else % we can superimpose by simple frame transformation
            ontarget = false;
            efn = 0; % counter for possible final anchor fragments
            while ~ontarget && efn < length(nindices)
                efn = efn + 1;
                efcode = nindices(efn); % fragment index of anchor fragment
                [match,mpoi] = min(abs(ecodes-efcode));
                while match ~= 0 && efn < length(nindices)
                    efn = efn +1;
                    efcode = nindices(efn);
                    [match,mpoi] = min(abs(ecodes-efcode));
                end
                if match ~= 0
                    err = -6;
                    break
                end
                % determine code of anchor fragment
                anchorcode = ecodes(mpoi);
                % pseudo-torsion defining atom coordinates of final
                % anchor fragment C4'(j),P(j+1),C4'(j+1),P(j+2)
                efrag = [fragments(anchorcode).A.coor(fragments(anchorcode).A.assign.previous(1:3),:);...
                    fragments(anchorcode).A.coor(fragments(anchorcode).A.assign.next(2),:);];
                % extract the ones that must superimpose with anchor
                % fragment
                fanchor = efrag(2:4,:);
                % find the affine transformation that superimposes the fragment to
                % the anchor
                [rmsd, ~, ftransmat] = superimpose_3points(anchore,fanchor);
                if rmsd < minrmsd
                    minrmsd = rmsd;
                end
                % try the next allowed anchor fragment or next stem, if superposition is poor
                if rmsd > options.anchor_acc
                    continue
                end
                % transform anchor fragment to anchor coordinates
                efrag = [efrag ones(4,1)];
                efrag = ftransmat*efrag';
                efrag = efrag';
                % extract the coordinates that must fit the stem fragment,
                % these are C4'(j),P(j+1),C4'(j+1)
                stemfinal = [efrag(1,1:3); anchore(1:2,:)];
                escode = stems(counter(1)).efragment;
                % extract stem fragment coordinates
                % C4'(j-1),P(j),C4'(j),P(j+1),C4'(j+1)
                sfrag = [fragments(escode).A.coor(fragments(escode).A.assign.previous(1:3),:);...
                    fragments(escode).A.coor(fragments(escode).A.assign.next(2:3),:);];
                % extract the ones that must superimpose with the anchor
                % nucleotide, these are C4'(j),P(j+1),C4'(j+1)
                [rmsd, ~, stransmat] = superimpose_3points(stemfinal,sfrag(3:5,:));
                % try the next allowed anchor fragment or next stem, if superposition is poor
                if rmsd < minrmsd2
                    minrmsd2 = rmsd;
                end
                if rmsd > 2*options.anchor_acc
                    continue
                end
                % transform final stem fragment to stem frame
                sfrag = [sfrag ones(5,1)];
                sfrag = sfrag*stransmat';
                % extract the atoms that must superimpose with the final
                % stem nucleotide C4'(j-1),P(j),C4'(j)
                stemfinal = sfrag(1:3,1:3);
                fixed = sfrag(2,:);
                % compute the transformation matrix for the stem
                [~, ~, stransmat] = superimpose_3points(stemfinal,ptal(:,1:3));
                cstem = cstem*stransmat';
                ontarget = true;
            end
            if ~ontarget
                success = false;
            else
                success = true;
            end
        end
    end
    if success
        ecoor = [ntnum cstem(:,1:3)];
    end
    if correction
        ecoor = anchor_correction(ecoor,corr);
    end
    % try to repair clashes, if requested
    if success && clash_test
        coor = ecoor(:,2:4);
        options.max_shift = 2;
        options.clash_thr = 1.4;
        [coor,min_dist,cycles,max_shift] = clash_repair(coor,environ,fixed,options);
        ecoor(:,2:4) = coor;
%        fprintf(1,'Clash repair required %i cycles and lead to a maximum shift of %6.2f �\n',cycles,max_shift);
        if min_dist > max_min_dist
            max_min_dist = min_dist;
        end
        if min_dist < options.clash_thr
%             pair_dist = get_all_pair_dist(coor,environ);
%             [mdv,mdp] = min(pair_dist);
%             [min_dist,mdp2] = min(mdv);
%             mdp1 = mdp(mdp2);
%             fprintf(1,'RNA loop atom %s of nt %i clashed with environment atom %i at distance %6.3f �.\n',atomtags{mdp1},ecoor(mdp1,1),mdp2,min_dist);
            err = -1;
            success = false;
        end
        if max_shift > options.max_shift
%             fprintf(1,'Clash repair required a shift of %5.2f � exceeding the limit of %5.2f �\n',max_shift,options.max_shift);
            err = -7;
            success = false;
        end
    end
    % count up to the next trial
    [counter,pindices,nindices,failed] = count_up(counter,stems,scramble);
end

% fprintf(1,'Maximum closest approach was at %5.2f �\n',max_min_dist);

if failed
    err = -8;
    ecoor = [];
    code = nan(1,2);
end

if ~success
    ecoor = [];
    code = nan(1,2);
    if err == 0
        err = -9; % catch cases, where no error code was assigned
    end
    return
end

% fprintf(1,'Fixed atom has coordinates (%5.3f, %5.3f, %5.3f) �\n',fixed);
if options.save
    fname = 'test_RNA';
    fname = wr_pdb_RNA(fname,ecoor,atomtags,stem_def.seq);
    fprintf(1,'Test RNA written to %s\n',fname);
end


function [counter,pindices,nindices,failed] = count_up(counter,stems,scramble)

pindices = [];
nindices = [];
failed = false;

counter(2) = counter(2) + 1;
if counter(2) > length(scramble.previous{counter(1)})
    counter(2) = 1;
    counter(1) = counter(1) + 1;
end
if counter(1) > length(stems)
    failed = true;
else
    pindices = scramble.previous{counter(1)};
    nindices = scramble.next{counter(1)};
end