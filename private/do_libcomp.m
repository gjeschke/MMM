function [err,msg] = do_libcomp(adr,lib1,lib2,infile,outfile)
% [err,msg] = do_libcomp(adr,lib1,lib2,infile,outfile)
%
% Compares in silico spin labeling with two or more different rotamer libraries
% libraries can be for the same label type or for different types
% either site address adr, lib1, and lib2 or infile must specify the task
% for two libraries, a difference analysis is performed
% for more than two libraries, specified in infile, mean value and standard
% deviation of the predictions are computed
%
% adr       MMM address of the site(s) to be compared, can also be a whole
%           structure, chain, or chain model, if empty, information
%           is read from infile
% lib1      first library file name (extension .mat can be present or
%           missing), if empty, information is read from infile
% lib2      second library file name (extension .mat can be present or
%           missing), if empty, information is read from infile
% infile    name of input specification file, can be missing or empty if
%           adr, lib1, and lib2 are specified, is ignored in that case
%           extension .dat is appended, if none exists
%           format: comment lines starting with % are allowed
%           # LIBCOMP lib1 lib2  % this key line must be present
%           adr1 % at least one address line must be present
%           adr2
%           ...
%           adrn
%           # END
% outfile   name of output file, extension .dat is appended, if none
%           exists, if missing or empty, output goes to message board
%
% err       error code, 0 no error
% msg       plain text error message
%
% G. Jeschke, 2014

global general

err = 0;
msg = 'No error.';

% if any of the first three arguments is missing, input specification comes
% from infile, otherwise an infile is written to the tmp subdirectory
if ~isempty(adr) && ~isempty(lib1) && ~isempty(lib2),
    infile = strcat(general.tmp_files,'libcomp_input.dat');
    fid = fopen(infile,'wt');
    if fid == -1,
        err = 1;
        msg = 'Input specification file could not be written to tmp directory';
        return
    end;
    fprintf(fid,'%% Temporary input file for libcomp\n\n');
    fprintf(fid,'# LIBCOMP %s %s\n',lib1,lib2);
    fprintf(fid,'%s\n',adr);
    fprintf(fid,'# END\n');
    fclose(fid);
else
    if ~exist('infile','var') || isempty(infile),
        err = 2;
        msg = 'Input specification for LIBCOMP was incomplete.';
        return
    end;
end;

% Automatic file extension for infile
if isempty(strfind(infile,'.')),
    infile = strcat(infile,'.dat');
end;

if exist('outfile','var') && ~isempty(outfile),
    % automatic extension for output file
    if isempty(strfind(outfile,'.')),
        outfile = strcat(outfile,'.dat');
    end;
    ofid = fopen(outfile,'wt');
    if ofid == -1,
        err = 3;
        msg = 'Output file could not be written';
        return
    end;
else
    ofid = 1;
end;

ifid=fopen(infile);
if ifid==-1,
    err = 4;
    msg = 'Input file could not be read';
    return;
end;

modus = 0;
while modus >= 0
    tline = fgetl(ifid);
    if ~ischar(tline) || modus<0, break, end
    if ~isempty(tline),
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k),
            if k(1)>1,
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end;
        end;
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#'),
            switch upper(char(args(2)))
                case 'LIBCOMP'
                    modus=1;
                    if length(args) < 4,
                        fclose(ifid);
                        if ofid ~=1,
                            fclose(ofid);
                        end;
                        err = 5;
                        msg = 'LIBCOMP keyword in input file must have at least two library name arguments';
                        return
                    end;
                    mprintf(ofid,'--- MMM rotamer library comparison ---\n\n');
                    for kla = 3:length(args),
                        libs{kla-2} = char(args(kla));
                        if ~strfind(libs{kla-2},'.'),
                            libs{kla-2} = strcat(libs{kla-2},'.mat');
                        end;
                        mprintf(ofid,sprintf('Lib%i: %s\n',kla-2,libs{kla-2}));
                    end;
                    softlinefeed(ofid);
                case 'END'
                    modus=-1;
                otherwise
                    modus=0;
                    mprintf(ofid,'Unknown keyword: %s in LIBCOMP input.\n',char(args(2)));
            end;
        elseif modus == 1,
            type = 4;
            adr = char(args(1));
            % extend structure addresses to all residues of this structure
            if adr(end) == ']',
                adr = strcat(adr,'(:){:}:');
                type = 1;
            end;
            % extend chain addresses to all residues of this chain
            if adr(end) == ')',
                adr = strcat(adr,'{:}:');
                type = 2;
            end;
            % extend chain model addresses to all residues of this chain
            % model
            if adr(end) == '}',
                adr = strcat(adr,':');
                type = 3;
            end;
            [rindices,message] = resolve_address(adr); 
            if message.error,
                mprintf(ofid,sprintf('Address %s could not be resolved.',adr));
                mprintf(ofid,sprintf('Error text: %s',message.text));
            else
                [m,n] = size(rindices);
                if n == 4,
                    mprintf(ofid,sprintf('Object %s has %i residues',adr,m));
                    softlinefeed(ofid);
                    for kla = 1:length(libs),
                        load(libs{kla});
                        labels{kla} = rotamer_populations(rindices,rot_lib);
                    end;
                    analyze_rotamers(rindices,labels,ofid,type);
                else
                    mprintf(ofid,sprintf('(Extended) address %s is not on residue level.',adr));
                    softlinefeed(ofid);
                end;
            end;
        end;
    end;
end;

fclose(ifid);

if ofid ~=1,
    fclose(ofid);
end;

function softlinefeed(ofid)

if ofid ~=1,
    fprintf(ofid,'\n');
end;

function analyze_rotamers(rindices,labels,ofid,type)

global model

[m,~] = size(rindices);
% check for all residues if labeled and make by-residue comparison
mprintf(ofid,'Residue-by-residue comparison');
softlinefeed(ofid);
pl = ones(1,length(labels));
comparison = zeros(m,4+4*length(labels));
poi = 0;
found = zeros(1,length(labels));
rmsd = zeros(1,length(labels));
for k = 1:m,
    adr = mk_address(rindices(k,:),1);
    break_it = false;
    for kla = 1:length(labels),
        if ~sum(abs(labels{kla}(pl(kla)).indices - rindices(k,:))), % labeling with lib1 was successful
            [lp,crmsd] = get_label_info(labels{kla},pl(kla));
            label_pos{kla} = lp;
            rmsd(kla) = crmsd;
            pl(kla) = pl(kla) + 1;
            if pl(kla) > length(labels{kla}),
                break_it = true;
            end;
            mprintf(ofid,sprintf('Lib%i: Residue %25s at [%5.2f,%5.2f,%5.2f] Å with r.m.s.d. %5.2f Å',...
                kla,adr,lp(1),lp(2),lp(3),crmsd));
            found(kla) = true;
        else
            mprintf(ofid,sprintf('Lib%i: Residue %25s could not be labeled.',kla,adr));
            found(kla) = false;
        end;
    end;
    if sum(found) == length(labels), % residue could be labeled with all libraries
        poi = poi + 1;
        comparison(poi,1:4) = rindices(k,:);
        for kla = 1:length(labels),
            comparison(poi,1+4*kla:3+4*kla) = label_pos{kla};
            comparison(poi,4+4*kla) = rmsd(kla);
        end;
    end;
    if break_it,
        break;
    end;
end;
comparison = comparison(1:poi,:);
softlinefeed(ofid);
if length(labels) == 2, % for two libraries, make a difference analysis
    mprintf(ofid,'Differences');
    diffs = zeros(1,poi);
    rmsd_diffs = zeros(1,poi);
    for k = 1:poi,
        adr = mk_address(comparison(k,1:4),1);
        diff = norm(comparison(k,5:7)-comparison(k,9:11));
        diffs(k) = diff;
        rmsd_diff = comparison(k,8)-comparison(k,12);
        rmsd_diffs(k) = rmsd_diff;
        mprintf(ofid,sprintf('%25s: %5.2f Å, r.m.s.d. diff: %6.3f',adr,diff,rmsd_diff));
    end;
    if type < 4, % analysis by chain model
        softlinefeed(ofid);
        mprintf(ofid,'Analysis by chain models');
        softlinefeed(ofid);
        % make list of unique chain models and their assignments
        cid = 1000000*comparison(:,1) + 1000*comparison(:,2) + comparison(:,3);
        ids = zeros(1,m);
        assignment = zeros(m,1);
        cmpoi = 0;
        for k = 1:length(cid),
            if ~min(abs(ids - cid(k))), % new chain model identifier
                cmpoi = cmpoi + 1;
                ids(cmpoi) = cid(k);
            end;
            [~,cpoi] = min(abs(ids - cid(k)));
            assignment(k) = cpoi;
        end;
        for k = 1:cpoi,
            this_indices = comparison((assignment == k),1:4);
            sid = this_indices(1,1);
            cid = this_indices(1,2);
            mid = this_indices(1,3);
            this_coor = model.structures{sid}(cid).xyz{mid};
            center = mean(this_coor);
            this_model_1 = comparison((assignment == k),5:7);
            this_model_2 = comparison((assignment == k),9:11);
            [mr,~] = size(this_model_1);
            this_model_1 = this_model_1 - reprowvector(center,mr);
            rad_gyr1 = sqrt(sum(sum(this_model_1.^2))/mr);
            this_model_2 = this_model_2 - reprowvector(center,mr);
            rad_gyr2 = sqrt(sum(sum(this_model_2.^2))/mr);
            this_diff = sqrt(sum(diffs(assignment == k).^2)/mr);
            this_rmsd_diff = mean(rmsd_diffs(assignment == k));
            adr = mk_address(this_indices(1,1:3));
            mprintf(ofid,sprintf('%s: Mean position difference %5.2f Å, mean r.m.s.d. difference: %5.2f Å',adr,this_diff,this_rmsd_diff));
            mprintf(ofid,sprintf('Label radii of gyration are Lib1: %6.2f Å, Lib2: %6.2f Å',rad_gyr1,rad_gyr2));
            softlinefeed(ofid);
        end;
    end;
else % analysis for more than two libraries
    mprintf(ofid,'Mean coordinates and prediction uncertainties');
    uncertain = zeros(poi,5);
    for k = 1:poi,
        uncertain(k,1) = comparison(k,1);
        adr = mk_address(comparison(k,1:4),1);
        mean_coor = zeros(1,3);
        for kla = 1:length(labels),
            mean_coor = mean_coor + comparison(k,1+4*kla:3+4*kla);
        end;
        mean_coor = mean_coor/length(labels);
        diff = 0;
        for kla = 1:length(labels),
            diff = diff + norm(comparison(k,1+4*kla:3+4*kla) - mean_coor)^2;
        end;
        rmsd_diff = sqrt(diff/length(labels));
        mprintf(ofid,sprintf('%25s: [%5.2f,%5.2f,%5.2f] Å with r.m.s.d. %5.2f Å',adr,mean_coor,rmsd_diff));
        uncertain(k,2:4) = mean_coor;
        uncertain(k,5) = rmsd_diff;
    end;
end;
save uncertainties_CASD uncertain

function [label_pos,rmsd] = get_label_info(calc_positions,poi)

pop = calc_positions(poi).NOpos(:,4);
coor = calc_positions(poi).NOpos(:,1:3);
label_pos = pop'*coor;
[mc,~] = size(coor);
diff = coor - reprowvector(label_pos,mc);
varvec = sum(diff.^2,2);
rmsd = sqrt(pop'*varvec);
