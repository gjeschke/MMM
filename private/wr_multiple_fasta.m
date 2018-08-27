function message=wr_multiple_fasta(fname,alignment,selection)
% function wr_multiple_fasta(fname,alignment,selection)
%
% generates a multiple sequence file in FASTA format when given a multiple 
% sequence record or alignment and, optionally, a list of
% sequence numbers to be included, this can be used, among else, for
% creating input for the sequence alignment program MUSCLE by R. C. Edgar
%
% alignments can be read with 
% get_multiple_fasta.m
% get_multiple_pir.m
% get_multiple_clustal.m
%
% fname     full file name, including extension
% alignment array of alignment struct, as created by the read routines mentioned above
%           as a minimum, each item must have fields
%           .db
%           .sequence
%           .name
% selection optional list of sequence numbers (items in alignment) that are
%           included into output
%           
% message   error message record (message.error and message.text) 
%
% Remarks:
% - assumed format for the so-called FASTA defline is based on the 
%   "Sequence identifiers" section of
%   http://en.wikipedia.org/wiki/FASTA_format (as of 21.03.2011)
%
%   generally, the FASTA defline is ill-defined, but major databases use
%   well-defined subformats, for pdb we use the NCBI Blast flavor of the
%   defline, rather than the very non-standard format provided by the PDB
%   homepage
%   sister routine get_multiple_fasta.m accepts both PDB defline formats
%
% G. Jeschke, 2011


message.error=0;
message.text='No error.';

if isempty(alignment),
    message.error=1;
    message.text='Empty alignment.';
    return
end;

fid=fopen(fname,'wt');
if fid==-1,
    message.error=2;
    message.text='File could not be written';
    return;
end;

for cseq=1:length(alignment),
    alg=alignment(cseq);
    % defline preamble, which is independent of database (replacing the
    % messy PDB FASTA format by a sensible one, as NCBI does)
    fprintf(fid,'>%s|',alg.db);
    % rest of defline depends on data base
    switch alg.db,
        case 'gi'
            fprintf(fid,'%s|%s|%s| %s\n',alg.gi_number,alg.type,alg.accession,alg.name);
        case 'pir'
            fprintf(fid,'|%s\n',alg.name);
        case 'prf'
            fprintf(fid,'|%s\n',alg.name);
        case 'sp'
            fprintf(fid,'%s|%s\n',alg.accession,alg.name);
        case 'pdb'
            fprintf(fid,'%s|%s\n',alg.name,alg.chain);
        case 'pat'
            fprintf(fid,'%s|%s\n',alg.country,alg.name);
        case 'bbs'
            fprintf(fid,'%s\n',alg.name);
        case 'gnl'
            fprintf(fid,'%s|%s\n',alg.type,alg.name);
        case 'ref'
            fprintf(fid,'%s|%s\n',alg.accession,alg.name);
        case 'tr'
            fprintf(fid,'%s|%s\n',alg.accession,alg.name);
        case 'lcl'
            fprintf(fid,'%s\n',alg.name);
        otherwise
            fprintf(fid,'%s\n',alg.name);
    end;
    % write sequence, up to 70 residues per line
    tsequence=alg.sequence;
    while length(tsequence)>70,
        fprintf(fid,'%s\n',tsequence(1:70));
        tsequence=tsequence(71:end);
    end;
    fprintf(fid,'%s\n\n',tsequence);
end;

fclose(fid);