function optimize_by_yasara(fname,options)

if ~exist('options','var') || isempty(options)
    options.console = false;
end

if ~isfield(options,'console')
    options.console = false;
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
else
    yname = strcat(bname,'_yasara.pdb');
end

mypath = pwd;
cd(yasara_path);

iname = which('minimization_server.mcr');
ifid = fopen(iname,'r');
ofid = fopen('MMM.mcr','wt');
fprintf(ofid,'LoadPDB %s\n',ffname);
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
correct_yasara('MMM','','',yname);
delete('MMM.pdb');
delete('MMM.mcr');


cd(mypath);

