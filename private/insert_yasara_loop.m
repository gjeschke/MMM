function insert_yasara_loop(fname,Nanchors,Canchors,seqs,options)

if ~exist('options','var') || isempty(options)
    options.console = true;
    options.optimize = false;
end

if ~isfield(options,'console')
    options.console = true;
end

if ~isfield(options,'optimize')
    options.optimize = false;
end

yasara = which('yasara.exe');
if isempty(yasara)
    add_msg_board('ERROR: Yasara is not on Matlab path');
    return
end

yasara_path = fileparts(yasara); 

[fpath,bname,ext] = fileparts(fname);
if isempty(ext)
    ext = '.pdb';
end
if isempty(fpath)
    fpath = pwd;
end

ffname = fullfile(fpath,strcat(bname,ext));

if isfield(options,'fname') && ~isempty(options.fname)
    yname = options.fname;
elseif options.optimize
    yname = strcat(bname,'_yasara.pdb');
else
    yname = strcat(bname,'_loops.pdb');
end

mypath = pwd;
cd(yasara_path);

for kseg = 1:length(seqs) % loop over all segments
    
    seq = seqs{kseg};
    Nanchor = Nanchors{kseg};
    Canchor = Canchors{kseg};
    
    [atnums,slc,~,chains] = get_atnums(ffname,Nanchor,Canchor);
    if min(atnums) < 1
        add_msg_board('ERROR insert_yasara_loop: Not all anchor backbone atoms found.');
        return
    end
    
    fseq = [slc(1) seq slc(2)];
    
    fid = fopen('MMM.mcr','wt');
    
    fprintf(fid,'LoadPDB %s\n',ffname);
    fprintf(fid,'BuildLoop %i %i %i,Sequence=%s,%i %i %i,Structures=1,Mutate=All,Bumpsum=1.0\n'...
        ,atnums(1),atnums(2),atnums(3),fseq,atnums(4),atnums(5),atnums(6));
    % fprintf(fid,'SelectRes %i-%i\n',res(1),res(2));
    fprintf(fid,'SavePDB 1,MMM.pdb,Format=PDB,Transform=Yes\n'); % Selected
    fprintf(fid,'Exit\n');
    fclose(fid);
    
    cmd = sprintf('yasara MMM.mcr');
    if options.console
        cmd = strcat(cmd,' -con');
    end
    [s, w] = dos(cmd);
    if s~=0
        cd(mypath);
        add_msg_board(w);
        return
    end
    delete('MMM.mcr');
    ffname = 'MMM_corr.pdb';
    if ~options.optimize && kseg == length(seqs)
        correct_yasara('MMM',chains(2),chains(1),yname);
    else
        correct_yasara('MMM',chains(2),chains(1));
    end
    delete('MMM.pdb');
end


if options.optimize
    iname = which('minimization_server.mcr');
    ifid = fopen(iname,'r');
    ofid = fopen('MMM.mcr','wt');
    fprintf(ofid,'LoadPDB MMM_corr.pdb\n');
    while 1
        tline = fgetl(ifid);
        if ~ischar(tline), break, end
        fprintf(ofid,'%s\n',tline);
    end
    % fprintf(ofid,'SelectRes %i-%i\n',res(1),res(2));
    fprintf(ofid,'SavePDB 2,MMM.pdb,Format=PDB,Transform=Yes\n');
    fprintf(ofid,'Exit\n');
    fclose(ifid);
    fclose(ofid);
    cmd = sprintf('yasara MMM.mcr');
    if options.console
        cmd = strcat(cmd,' -con');
    end
    [s, w] = dos(cmd);
    if s~=0
        cd(mypath);
        add_msg_board(w);
        return
    end
    correct_yasara('MMM',chains(2),chains(1),yname);
    delete('MMM.pdb');
    delete('MMM.mcr');
end


cd(mypath);

function [atnums,slc,res,chains] = get_atnums(fname,Nanchor,Canchor)

global residue_defs

slc = '  ';
atnums = zeros(1,6);
ifid = fopen(fname);

p1 = strfind(Nanchor,'(');
p2 = strfind(Nanchor,')');
Nchain = Nanchor(p1+1:p2-1);
Nres = str2double(Nanchor(p2+1:end));

p1 = strfind(Canchor,'(');
p2 = strfind(Canchor,')');
Cchain = Canchor(p1+1:p2-1);
Cres = str2double(Canchor(p2+1:end));
res = [Nres Cres];
chains = [Nchain Cchain];

if ifid==-1
    add_msg_board('ERROR: PDB file does not exist');
    return;
end

offset = 0;

while 1
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if strcmpi(tline(1:3),'TER')
        offset = offset -1;
    end
    if length(tline) >= 4 && strcmp(tline(1:4),'ATOM')
        resnum = str2double(tline(23:26));
        if resnum == Nres
            if strcmpi(tline(22),Nchain)
                if strcmpi(tline(14:15),'N ')
                    atnums(1) = str2double(tline(7:11)) + offset;
                    tlc = tline(18:20);
                    rid = tag2id(tlc,upper(residue_defs.restags));
                    if ~isempty(rid)
                        slc(1) = residue_defs.single_letter_code(rid);
                    end
                end
                if strcmpi(tline(14:15),'CA')
                    atnums(2) = str2double(tline(7:11)) + offset;
                end
                if strcmpi(tline(14:15),'C ')
                    atnums(3) = str2double(tline(7:11)) + offset;
                end
            end
        end
        if resnum == Cres
            if strcmpi(tline(22),Cchain)
                if strcmpi(tline(14:15),'N ')
                    atnums(4) = str2double(tline(7:11)) + offset;
                    tlc = tline(18:20);
                    rid = tag2id(tlc,upper(residue_defs.restags));
                    if ~isempty(rid)
                        slc(2) = residue_defs.single_letter_code(rid);
                    end
                end
                if strcmpi(tline(14:15),'CA')
                    atnums(5) = str2double(tline(7:11)) + offset;
                end
                if strcmpi(tline(14:15),'C ')
                    atnums(6) = str2double(tline(7:11)) + offset;
                end
            end
        end
    end
end
fclose(ifid);