function fname = get_pdb_file(PDBID,silent)
% function fname = get_pdb_file
%
% retrieve a PDB file with given PDB identifier PDBID directly from the PDB 
% output is the full local file name (in the temporary directory)
% if download is not successful, the filename fname is empty
%
% if flag silent exists and is true, error output is suppressed
%
% G. Jeschke, 2009-2024

global general

if ~exist('silent','var') || isempty(silent)
    silent = false;
end

% PDB_ftp=web_adr.PDB;
% basdir=queries.PDB_structures;

curr_dir=pwd;

cd(general.tmp_files);

fname='';


query = sprintf('https://files.rcsb.org/download/%s.pdb.gz',lower(PDBID));
fname0 = [lower(PDBID) '.pdb.gz'];
try
    websave(fname0,query);
catch errid
    if ~silent
        add_msg_board('ERROR: ftp connection to PDB could not be opened.');
        add_msg_board(errid.message);
    end
    cd(curr_dir);
    return;
end

try
    gunzip(fname0);
catch errid
    cd(curr_dir);
    if ~silent
        add_msg_board(sprintf('ERROR: PDB file with ID "%s" could not be unzipped.',PDBID));
        add_msg_board(errid.message);
    end
    return;
end
try
    delete(fname0);
catch
end

fname=[general.tmp_files lower(PDBID) '.pdb'];
cd(curr_dir);