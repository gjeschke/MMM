function scan_pdb_for_MTSL(hours)
% Adds files to the coarse-grained PDB data base for a given run time
% internet access to the PDB server is required
%
% hours     scheduled runtime in hours

fname='PDB_directory_20121203.txt';

old_basis='PDB_directory_20121130.txt';

if nargin<1,
    hours=1;
end;

fid=fopen(old_basis);
if fid==-1,
    fprintf(2,'ERROR: Old basis file could not be opened.\n');
    return;
end;

tags=':';
fprintf(1,'Scanning existing database for PDB identifiers...\n');
while 1,
    tline = fgetl(fid);
    if ~ischar(tline) || length(tline)<7, break, end    
    pdbid=tline(4:7);
    pdbid=upper(pdbid);
    tags=[tags pdbid ':'];
end;
fprintf(1,'Existing database scanned. Now adding new data.\n');

fclose(fid);

% myname=which('add_to_database.m');
% mypath=fileparts(myname); 
dbpath='E:\coarse_PDB\';

oname=[dbpath 'MTSL_records.txt'];
ofid=fopen(oname,'a');

maxtime=hours*3600; % maximum time in seconds

fid=fopen(fname);
if fid==-1,
    fprintf(2,'ERROR: PDB directory could not be opened.\n');
    return;
end;

tic;
runtime=0;
success=0;
poi=0;
while runtime<=maxtime,
    tline = fgetl(fid);
    if ~ischar(tline) || length(tline)<7, break, end    
    poi=poi+1;
    pdbid=tline(4:7);
    pdbid=upper(pdbid);
    k=strfind(tags,pdbid);
    if ~isempty(k),
        continue;
    end;
    fprintf(1,'%i> Checking PDB file %s...\n',success,pdbid);
    labeled = is_pdb_MTSL(pdbid);
    if isempty(labeled),
        fprintf(2,'Structure %s could not be read.\n',pdbid);
    elseif labeled,
        success=success+1;
        fprintf(1,'(Structure %s contains MTSL.\n',pdbid);
        fprintf(ofid,'%s\n',pdbid);
    end;
    if mod(poi,20)==0,
        runtime=toc;
        h=floor(runtime/3600);
        min=floor((runtime-3600*h)/60);
        sec=runtime-3600*h-60*min;
        fprintf(1,'Runtime %i h %i min %4.1f. %5.2f%% completed.\n',h,min,sec,100*runtime/maxtime);
    end;
end;
fclose(fid);
fclose(ofid);

runtime=toc;
h=floor(runtime/3600);
min=floor((runtime-3600*h)/60);
sec=runtime-3600*h-60*min;
fprintf(1,'\nTotal runtime %i h %i min %4.1f.\n',h,min,sec);

fprintf(1,'%i PDB structures were scanned, %i of them contained MTSL.\n',poi,success); 
fprintf(1,'Last PDB ID was %s.\n',pdbid); 

if ~ischar(tline),
    fprintf(1,'Warning: Directory exhausted. Database appears to be complete.\n');
    return
end;

if runtime>maxtime,
    fprintf(1,'Program stopped after exceeding alloted runtime.\n');
end;
