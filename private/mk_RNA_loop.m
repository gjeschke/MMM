function [ecoor,atomtags,seq,options,correction,err] = mk_RNA_loop(seq,fragments,...
    shortfrag,non_clash_table,anchor,acodes,transmat,ecodes,options,ntoffset,environ)
% function [ecoor,atomtags,seq,options] = mk_RNA_loop(seq,fragments,...
%    shortfrag,non_clash_table,anchor,acodes,transmat,ecodes,options)
%
% err   error code
%       -1  clashes could not be corrected
%       -2  run out of time
%       -3  distance between anchors too large for number of nucleotides

atomtags = cell(0,1);
ecoor = [];
max_rpnt = 7;
err = 0;
correction = [];

% Check whether the loop has two anchors, and if so, if the anchors are
% sufficiently close
if exist('transmat','var') && ~isempty(transmat) 
    if exist('anchor','var') && ~isempty(anchor)
        rpnt = norm(transmat(1:3,4)'-anchor(1,:))/(length(seq)-1); % distance per nt
        if rpnt > max_rpnt
            seq = seq(2:end-1);
%             fprintf(1,'Distance per nt is %4.1f Å\n',rpnt);
            err = -3;
            return
        end
    end
end

% if ntoffset > 40
%     profile viewer
%     disp('Aber hallo!');
% end

if ~exist('environ','var')
    environ = [];
end

if ~exist('ntoffset','var') || isempty(ntoffset)
    ntoffset = 0;
end

if ~exist('options','var') || ~isfield(options,'clash_thr')
    options.clash_thr = 1.4;
end

if ~isfield(options,'maxtime') || isempty(options.maxtime)
    options.maxtime = 0.1;
end

tstart = tic;
runtime = toc(tstart);

while isempty(ecoor) && runtime < 3600*options.maxtime
    
    [~,code,~,correction,lerr,statistics] = mk_RNA_loop_backbone(seq,...
        fragments,shortfrag,non_clash_table,anchor,acodes,transmat,ecodes,options);
    if lerr == 0
        % fprintf(1,'Successful RNA loop modelling after %5.1f s and %i trials\n',3600*statistics.runtime,statistics.trials);
    else
        fprintf(1,'RNA modelling failed after  %5.1f s and %i trials\n',3600*statistics.runtime,statistics.trials);
        for ke = 1:length(statistics.errors)
            if statistics.errors(ke) > 0.5/statistics.trials                
                fprintf(1,'error %i accounts for %6.2f%%\n',statistics.errnum(ke),100*statistics.errors(ke));
            end
        end
        if statistics.min_convg < 1e5
            fprintf(1,'Best convergence was %4.2f Å\n',statistics.min_convg);
        end
        if statistics.amin_rmsd < 1e5
            fprintf(1,'Best anchor rmsd was %4.2f Å\n',statistics.amin_rmsd);
        end
        if statistics.amin_rot < 1e5
            fprintf(1,'Minimum rotation was %4.2f degree\n',statistics.amin_rot);
        end
        if statistics.amin_shift < 1e5
            fprintf(1,'Minimum translation was %4.2f Å\n',statistics.amin_shift);
        end
        runtime = toc(tstart);
        continue
    end

    [ecoor,atomtags] = mk_RNA(fragments,seq(1:end-1),code(1:end-1),non_clash_table,ntoffset,transmat);
       
    if isempty(ecoor)
        % fprintf(1,'mk_RNA failed\n');
        runtime = toc(tstart);
        continue
    else
        valid = find(ecoor(:,1)>ecoor(1,1)); % atoms not in the initial anchor nt
        atomtags = atomtags(valid);
        ecoor = ecoor(valid,:);
    end

    if correction.corrected && ~isfield(options,'anchor0')
        ecoor = anchor_correction(ecoor,correction);
    end

        % try to repair clashes, if requested
    if ~isempty(environ) && ~isfield(options,'anchor0')
        coor = ecoor(:,2:4);
        if ~isempty(anchor)
            fixed = anchor(1,:);
        else
            fixed = zeros(0,3);
        end
        if ~isempty(transmat)
            f2 = transmat(1:3,4);
            fixed = [f2'; fixed];
        end
        options.max_shift = 2;
        [coor,min_dist,~,max_shift] = clash_repair(coor,environ,fixed,options);
        ecoor(:,2:4) = coor;
        % fprintf(1,'Clash repair required %i cycles and lead to a maximum shift of %6.2f Å\n',cycles,max_shift);
        if min_dist < options.clash_thr || max_shift > options.max_shift
%             pair_dist = get_all_pair_dist(coor,environ);
%             [mdv,mdp] = min(pair_dist);
%             [min_dist,mdp2] = min(mdv);
%             mdp1 = mdp(mdp2);
%             fprintf(1,'RNA loop atom %s of nt %i clashed with environment atom %i at distance %6.3f Å.\n',atomtags{mdp1},ecoor(mdp1,1),mdp2,min_dist);
            ecoor = [];
            err = -1;
        else
            err = 0;
        end
    end
    runtime = toc(tstart);
    
end

if isempty(ecoor) && runtime >= 3600*options.maxtime
    err = -2;
end

seq = seq(2:end-1);