function m = get_pdb_dir
% function fname = get_pdb_dir
%
% retrieve the current PDB directory
%
% ### This does not work. MATLAB ftp dir appears to have problems with too
% long file lists ###
%
% G. Jeschke, 2009-2012

global web_adr
global general
global queries

PDB_ftp=web_adr.PDB;
basdir=queries.PDB_structures;

curr_dir=pwd;

cd(general.tmp_files);


add_msg_board('Downloading current PDB directory...');
try
    ftp_obj=ftp(PDB_ftp,'anonymous','anonymous');
catch errid
    add_msg_board('ERROR: ftp connection to PDB could not be opened.');
    add_msg_board(errid.message);
    cd(curr_dir);
    return;
end;
try
    cd(ftp_obj,basdir);
catch errid
    cd(curr_dir);
    add_msg_board('ERROR: Change to PDB structure directory not successful.');
    add_msg_board(errid.message);
    return;
end;
try
    binary(ftp_obj);
catch errid
    cd(curr_dir);
    add_msg_board('ERROR: Change to ftp binary mode not successful.');
    add_msg_board(errid.message);
    return;
end;
try    
    m=dir(ftp_obj);
catch errid
    cd(curr_dir);
    add_msg_board('ERROR: PDB directory could not be retrieved.');
    add_msg_board(errid.message);
    return;
end;