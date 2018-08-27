function [message,outfile] = align_sequences(seqs,sindices,is_silent)
% function [message,outfile]=align_sequences(seqs,sindices,is_silent)
%
% aligns two peptide sequences by calling MUSCLE by Robert C. Edgar
% the alignment in Clustal output format is displayed in the report editor
% if MUSCLE is not found on the Matlab path, the routine returns with an
% error message
%
% seqs      cell array of peptide sequences (at least two) 
% sindices  (optional) array of MMM indices for the chains, use empty
%           array if you wish to provide flag 'is_silent' without providing
%           indices
% is_silent (optional) flag that suppresses report editor output, defaults
%           to false (alignment is shown in report editor)
%
% message   error message
% outfile   name of alignment output file (Clustal format)
%
% G. Jeschke, 2011

global model
global general
global help_files
global third_party
global hMain

message.error=0;
message.text='No error.';

outfile='';

if nargin<3,
    is_silent=false;
end;

for k=1:length(seqs), % replace unknown residue identifier to conform to common usage
    seq=seqs{k};
    for kk=1:length(seq),
        if char(seq(kk))=='?',
            seq(kk)='X';
        end;
    end;
    seqs{k}=seq;
end;

if length(seqs)<2,
    message.error=1;
    message.text='ERROR: At least two peptide sequences required for alignment.';
    return
end;

is_indexed=false;
if nargin>1 && ~isempty(sindices),
    [m,n]=size(sindices);
    if m==length(seqs),
        is_indexed=true;
    end;
end;

for k=1:length(seqs),
    alignment(k).db='MMM';
    if is_indexed
        alignment(k).name=mk_address(sindices(k,:));
    else
        alignment(k).name=sprintf('seq%i',k);
    end;
    alignment(k).sequence=seqs{k};
end;

dospath=which('muscle.exe');
entry=strcat(help_files,'third_party.html#MUSCLE');
if isempty(dospath),
    add_msg_board('ERROR: MUSCLE could not be found on the Matlab path.');
    add_msg_board('Please check whether MUSCLE is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;
% add the reference, if it does not yet exist
muscle_ref=true;
id=tag2id('Edgar:2004_muscle',third_party.tags,[],'|');
if ~isempty(id),
    if isfield(model,'auto_references'),
        if ~isempty(find(id==model.auto_references, 1)),
            muscle_ref=false;
        end;
    else
        model.auto_references=[];
    end;
    if muscle_ref,
        if ~isfield(model,'references'),
            model.references(1)=third_party.references(id);
        elseif isempty(model.references)
            model=rmfield(model,'references');
            model.references(1)=third_party.references(id);
        else
            model.references(end+1)=third_party.references(id);
        end;
        model.auto_references(end+1)=id;
    end;
end;
infile=fullfile(general.tmp_files,'doalign.fa');
wr_multiple_fasta(infile,alignment);
outfile=[general.tmp_files 'muscle_alignment.aln'];
comd=[dospath ' -in ' infile ' -clwout ' outfile];
if ispc
    [s, w] = dos(comd);
elseif isunix
    [s, w] = unix(comd);
end
if s~=0,
    message.error=3;
    message.text='MUSCLE error.';
    add_msg_board('ERROR: MUSCLE did not run successfully.');
    disp(w);
    return
else
    if ~is_silent,
        rem=w;
        while ~isempty(rem),
            [token,rem]=strtok(rem,char(10));
            if ~isempty(token),
                add_msg_board(token);
            end;
        end;    
        hMain.report_file=outfile;
        report_editor;
    end;
end;
