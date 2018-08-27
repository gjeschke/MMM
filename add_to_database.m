function add_to_database(hours)
% Adds files to the coarse-grained PDB data base for a given run time
% internet access to the PDB server is required
%
% hours     scheduled runtime in hours

fname='PDB_directory_20121130.txt';

% myname=which('add_to_database.m');
% mypath=fileparts(myname); 
dbpath='E:\coarse_PDB\';
load([dbpath 'dir_pointer']); % must define variables poi, total, empty, and pdbfiles

maxtime=hours*3600; % maximum time in seconds

fid=fopen(fname);
if fid==-1,
    fprintf(2,'ERROR: PDB directory could not be opened.\n');
    return;
end;

tline='';

nl=0;
while nl<poi,
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
end;

if ~ischar(tline) || length(tline)<7,
    fprintf(1,'Warning: Directory exhausted. Database appears to be complete.\n');
    fclose(fid);
    return
end;

tic;
runtime=0;
poi0=poi;
empty0=empty;
while runtime<=maxtime,
    tline = fgetl(fid);
    if ~ischar(tline) || length(tline)<7, break, end    
    poi=poi+1;
    pdbid=tline(4:7);
    pdbfiles(poi,:)=upper(pdbid);
    fname=[dbpath upper(pdbid) '_coarse'];
    fprintf(1,'Reading PDB file %s...\n',upper(pdbid));
    [protein,chain_tags,chain_ids]=rd_pdb_coarse(pdbid);
    save(fname,'protein','chain_tags','chain_ids');
    if isempty(protein),
        fprintf(2,'Structure %s returned empty records.\n',upper(pdbid));
        empty=empty+1;
    end;
    runtime=toc;
    h=floor(runtime/3600);
    min=floor((runtime-3600*h)/60);
    sec=runtime-3600*h-60*min;
    fprintf(1,'Runtime %i h %i min %4.1f. %5.2f%% completed.\n',h,min,sec,100*runtime/maxtime);
end;
fclose(fid);

save([dbpath 'dir_pointer'],'poi','total','empty','pdbfiles');

fprintf(1,'%i PDB structures were added to the data base, %i of them were empty.\n',poi-poi0,empty-empty0); 
fprintf(1,'%4.1f%% of the data base are now completed.\n',100*poi/total); 

if ~ischar(tline),
    fprintf(1,'Warning: Directory exhausted. Database appears to be complete.\n');
    return
end;

if runtime>maxtime,
    fprintf(1,'Program stopped normally after exceeding alloted runtime.\n');
end;
