function [fname,gaplist,message]=mk_pir(basname,tempname,template,sequence,cid,res1)
% function fname=mk_pir(snum,sequence)
%
% generates a sequence alignment file in PIR format, as needed for Modeller
% input
%
% basname   basis name for the output file (no directory, no extension),
%           file is saved in the tmp subdirectory
% tempname  name of the template PDB file
% template  indices of the template chain in the MMM model
% sequence  (optional) sequence for which Modeller should create a homology
%           model, default is the full sequence of the template chain,
%           including missing residues
% cid       (optional) chain identifier in the modelled structure, defaults
%           to the chain identifier in the template, must be an uppercase
%           letter or a space
% res1      (optional) first residue number in the modelled structure,
%           defaults to 1
%
% fname     full file name of the output file, extension is .ali
% gaplist   list of gaps in the template, empty if there are no gaps
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

global model
global general

extend=true;

message.error=0;
message.text='No error.';

cleanup=false;

gaplist=[];

snum=template(1);
cnum=template(2);

tsequence=model.structures{snum}(cnum).sequence;
if nargin<4,
    sequence=tsequence;
end;

if nargin<6,
    res1=1;
end;

if isfield(model.info{snum},'cryst'),
    modeltype='structureX';
else
    modeltype='structure';
end;

fname=[general.tmp_files basname '.ali'];
fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='Alignment file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

adr=mk_address(template);

fprintf(fid,'C; Alignment file based on template structure %s\n\n',adr);

[stag,ctag]=mk_address_parts(template);

if nargin<5,
    cid=ctag;
end;

% prepare sequence file of template with missing residue information
if ~isfield(model.structures{snum}(cnum),'seqexist') || length(model.structures{snum}(cnum).seqexist)~=length(tsequence),
    seqexist=zeros(1,length(tsequence));
    mres=length(model.structures{snum}(cnum).residues{1}.info);
    numbers=zeros(1,mres);
    for k=1:mres,
        numbers(k)=model.structures{snum}(cnum).residues{1}.info(k).number;
    end;
    for k=1:length(tsequence),
        if ~isempty(find(numbers==k, 1)),
            seqexist(k)=1;
        end;
    end;
    model.structures{snum}(cnum).seqexist=seqexist;
else
    seqexist=model.structures{snum}(cnum).seqexist;
end;
for k=1:length(seqexist),
    if ~seqexist(k),
        tsequence(k)='-';
    end;
    if strcmp(tsequence(k),'?'),
        tsequence(k)='.'; % '-'
    end;
end;
for k=1:length(sequence),
    if strcmp(sequence(k),'?'),
        sequence(k)='-'; % '-'
    end;
end;
tres1=1;
trese=length(tsequence);
rese=length(sequence);
tseq0=tsequence;
seq0=sequence;

% remove leading missing residues in template sequence
while strcmp(tsequence(1),'-'),
    tres1=tres1+1;
    res1=res1+1;
    tsequence=tsequence(2:end);
    sequence=sequence(2:end);
    cleanup=true;
end;
if ~extend,
    % remove trailing missing residues in template sequence
    while strcmp(tsequence(end),'-'),
        trese=trese-1;
        rese=rese-1;
        tsequence=tsequence(1:end-1);
        sequence=sequence(1:end-1);
        cleanup=true;
    end;    
end;

if ~cleanup,
    tsequence=tseq0;
    sequence=seq0;
    rshift=0;
else
    rshift=tres1-1;
end;

% prepare gap list
gaplist=zeros(100,2);
poi=0;
gap=false;
for k=1:length(tsequence),
    if ~gap && strcmp(tsequence(k),'-'),
        gap=true;
        poi=poi+1;
        gaplist(poi,1)=rshift+k;
    end;
    if gap && ~strcmp(tsequence(k),'-'),
        gap=false;
        gaplist(poi,2)=rshift+k-1;
    end;
end;
if gap,
    gaplist(poi,2)=rshift+k;
end;
gaplist=gaplist(1:poi,:);

if isfield(model.info{snum},'molecule'),
    molecule=model.info{snum}.molecule;
    if ~isempty(strfind(molecule,':')),
        molecule='undefined';
    end;
else
    molecule='undefined';
end;

if isfield(model.info{snum},'organism'),
    organism=model.info{snum}.organism;
    if ~isempty(strfind(organism,':')),
        organism='undefined';
    end;
else
    organism='undefined';
end;

if isfield(model.info{snum},'resolution') && ~isempty(model.info{snum}.resolution),
    resolution=model.info{snum}.resolution;
else
    resolution=-1;
end;
rfactor=-1;

% write template record
fprintf(fid,'>P1;tmp\n');
fprintf(fid,'%s:%s:%5i:%s:%5i:%s:%s:%s:%4.2f:%4.2f\n',modeltype,tempname,tres1,ctag,trese,ctag,molecule,organism,resolution,rfactor);
while length(tsequence)>70,
    fprintf(fid,'%s\n',tsequence(1:70));
    tsequence=tsequence(71:end);
end;
fprintf(fid,'%s*\n\n',tsequence);

modeltype='sequence';
if nargin>=4, % the sequence is not the same as in the template
    molecule='undefined';
    organism='undefined';
end;
resolution=-1;
rfactor=-1;
% write target record
fprintf(fid,'>P1;fill\n');
fprintf(fid,'%s:%s:%5i:%s:%5i:%s:%s:%s:%4.2f:%4.2f\n',modeltype,basname,res1,cid,rese,cid,molecule,organism,resolution,rfactor);
while length(sequence)>70,
    fprintf(fid,'%s\n',sequence(1:70));
    sequence=sequence(71:end);
end;
fprintf(fid,'%s*\n\n',sequence);
fclose(fid);
