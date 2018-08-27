function fname = get_pdb_file(PDBID,silent)
% function fname = get_pdb_file
%
% retrieve a PDB file with given PDB identifier PDBID directly from the PDB
% ftp server 
% output is the full local file name (in the temporary directory)
% if download is not successful, the filename fname is empty
%
% if flag silent exists and is true, error output is suppressed
%
% G. Jeschke, 2009

global web_adr
global general
global queries

if ~exist('silent','var') || isempty(silent)
    silent = false;
end

PDB_ftp=web_adr.PDB;
basdir=queries.PDB_structures;

curr_dir=pwd;

cd(general.tmp_files);

fname='';

% add_msg_board(sprintf('Downloading PDB file %s from FTP server...',PDBID));
fname0=['pdb' lower(PDBID) '.ent.gz'];
try
    ftp_obj=ftp(PDB_ftp,'anonymous','anonymous');
catch errid
    if ~silent
        add_msg_board('ERROR: ftp connection to PDB could not be opened.');
        add_msg_board(errid.message);
    end
    cd(curr_dir);
    return;
end;
try
    cd(ftp_obj,basdir);
catch errid
    cd(curr_dir);
    if ~silent
        add_msg_board('ERROR: Change to PDB structure directory not successful.');
        add_msg_board(errid.message);
    end
    return;
end;
try
    binary(ftp_obj);
catch errid
    cd(curr_dir);
    if ~silent
        add_msg_board('ERROR: Change to ftp binary mode not successful.');
        add_msg_board(errid.message);
    end
    return;
end;
try    
    mget(ftp_obj,fname0);
catch errid
    cd(curr_dir);
    if ~silent
        add_msg_board(sprintf('ERROR: Download of PDB file with ID "%s" not successful.',PDBID));
        add_msg_board(errid.message);
    end
    return;
end;
try
    gunzip(fname0);
catch errid
    cd(curr_dir);
    if ~silent
        add_msg_board(sprintf('ERROR: PDB file with ID "%s" could not be unzipped.',PDBID));
        add_msg_board(errid.message);
    end
    return;
end;
try
    delete(fname0);
catch errid
end;

% add_msg_board(sprintf('PDB file with identifier %s was downloaded.',PDBID));
fname=[general.tmp_files 'pdb' lower(PDBID) '.ent'];
cd(curr_dir);