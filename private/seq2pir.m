function [fname,message]=seq2pir(basname,tempname,targname,tempseq,ares,achain,eres,echain,targseq,arestarg)
% function [fname,message]=seq2pir(basname,tempname,targname,tempseq,targseq)
%
% generates a sequence alignment file in PIR format, as needed for Modeller
% input, given a template sequence and an optional target sequence
%
% basname   basis name for the output file (no directory, no extension),
%           file is saved in the tmp subdirectory
% tempname  name of the template PDB file
% targname  name of the target structure
% tempseq   sequence of the template structure in Modeller format, use
%           wr_pdb_paradigm to prepare the template PDB file and obtain
%           this sequence
% ares      first residue in template sequence
% achain    chain identifier for first residue in template sequence
% eres      last residue in template sequence
% echain    chain identifier for last residue in template sequence
% targseq   (optional) sequence of the target structure, if missing or 
%           empty, the target has 100% sequence identity, which makes sense 
%           only when providing restraints
% arestarg  (optional) first residue number for target chain, not used by
%           Modeller, but needed for consistent documentation, defaults to
%           ares
%
% fname     full file name of the output file, extension is .ali
% message   error message with fields .error (number) and .text, .error=0
%           for no eror
%
% Remarks:
% - missing residues in the template file are codes as '-' in the templae
%   sequence
% - the first residue of the template sequence is always residue 1
% - basname is also used as protein identifier for the target structure
%
% G. Jeschke, 2011

global general

message.error=0;
message.text='No error.';

if nargin<9 || isempty(targseq),
    targseq=tempseq;
end;

if nargin<10 || isempty(arestarg) || arestarg<1,
    arestarg=ares;
end;

if strcmpi(achain,echain),
    erestarg=eres+arestarg-ares;
else
    erestarg=eres;
end;

modeltype='structure';

fname=[general.tmp_files basname '.ali'];
fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='Alignment file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'C; Alignment file based on template structure %s\n\n',tempname);

molecule='undefined';
organism='undefined';
resolution=-1;
rfactor=-1;

% write template record
fprintf(fid,'>P1;%s\n',tempname);
fprintf(fid,'%s:%s:%5i:%s:%5i:%s:%s:%s:%4.2f:%4.2f\n',modeltype,tempname,ares,achain,eres,echain,molecule,organism,resolution,rfactor);
tsequence=tempseq;
while length(tsequence)>70,
    fprintf(fid,'%s\n',tsequence(1:70));
    tsequence=tsequence(71:end);
end;
fprintf(fid,'%s*\n\n',tsequence);

modeltype='sequence';
sequence=targseq;
% write target record
fprintf(fid,'>P1;%s\n',targname);
fprintf(fid,'%s:%s:%5i:%s:%5i:%s:%s:%s:%4.2f:%4.2f\n',modeltype,targname,arestarg,achain,erestarg,echain,molecule,organism,resolution,rfactor);
while length(sequence)>70,
    fprintf(fid,'%s\n',sequence(1:70));
    sequence=sequence(71:end);
end;
fprintf(fid,'%s*\n\n',sequence);
fclose(fid);
