function [dssp,s,w] = get_dssp(fname)
% dssp = get_dssp(fname)
%
% gets DSSP information on a PDB file, empty if no success
%
% fname name of input PDB file, with or without extension
%
% dssp  array of struct containing the DSSP information
% s     status, 0 if successful
% w     console output or error message
%
% G. Jeschke, 10.01.2020


dssp=[];
dospath=which('dssp-2.0.4-win32.exe');
if isempty(dospath)
    dospath=which('dssp.exe');
end
if isempty(dospath)
    dospath=which('dsspcmbi.exe');
end
if ~isempty(dospath) % suppress this if DSSP not known or MMM rotamer format or direct input
    add_msg_board('DSSP geometry analysis is performed');
    infile=which(fname);
    poi=strfind(infile,'.pdb');
    if isempty(poi)
        poi=strfind(infile,'.ent');
    end
    if isempty(poi) || poi(1)<2
        basname=infile;
    else
        basname=infile(1:poi(1)-1);
    end
    outfile=[basname '.dssp'];
    cmd=[dospath ' ' infile ' ' outfile];
    [s,w] = dos(cmd);
    if s~=0
        add_msg_board('Warning: DSSP did not run successfully on this PDB file.');    
        % disp(w);
    else
        add_msg_board('Result file:');
        add_msg_board(outfile);
        dssp_file=outfile;
        dssp=rd_dssp(dssp_file);
    end
else
    s = -100;
    w = 'ERROR: DSSP not found on Matlab path';
end
