function message=wr_multiple_pir(fname,alignment,selection)
% function message=wr_multiple_pir(fname,alignment,selection)
%
% generates a sequence alignment file in PIR format, as needed for Modeller
% input when given a multiple sequence alignment and, optionally, a list of
% sequence numbers to be included
%
% alignments can be read with 
% get_multiple_fasta.m
% get_multiple_pir.m
% get_multiple_clustal.m
%
% fname     full file name, including extension
% alignment array of alignment struct, as created by the read routines mentioned above
%           as a minimum, each item must have fields
%           .db (data base, determines existing information)
%           .sequence
%           .name
%           for pdb inputs, field .chain is also used
%           for original pir inputs (.db=='pir_native'), fields .seqtype
%           and .description are also used
%           for Modeller-type input (.db=='Modeller'), all of the following
%           fields are used
%           .seqtype
%           .type
%           .pdb
%           .first_res
%           .first_chain
%           .last_res
%           .protein
%           .source
%           .resolution
%           .Rfactor
% selection optional list of sequence numbers (items in alignment) that are
%           included into output
%           
% message   error message record (message.error and message.text) 
%
% Remarks:
% - the length of all sequences (.sequence of all items in alignment)
%   should match, otherwise the routine tests, whethe the third-party
%   program MUSCLE by Robert C. Edgar is present, if so, sequences are
%   aligned with that program, if not, an error message is returned

%
% G. Jeschke, 2011

global model
global general
global help_files
global third_party

message.error=0;
message.text='No error.';

if isempty(alignment),
    message.error=1;
    message.text='Empty alignment.';
    return
end;

if nargin<3,
    selection=1:length(alignment);
end;

type='imported';
mismatch=0;
seqlength=length(alignment(selection(1)).sequence);
for k=2:length(selection),
    cseq=selection(k);
    if length(alignment(cseq).sequence)~=seqlength,
        mismatch=true;
    end;
end;
if mismatch,
    type='MUSCLE';
    add_msg_board('Warning: Sequence length mismatch.');
    add_msg_board('Alignment by MUSCLE will be attempted.');
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
    wr_multiple_fasta(infile,alignment,selection);
    outfile=[general.tmp_files 'isaligned.afa'];
    comd=[dospath ' -in ' infile ' -out ' outfile];
    [s, w] = dos(comd);
    if s~=0,
        message.error=3;
        message.text='MUSCLE error.';
        add_msg_board('ERROR: MUSCLE did not run successfully.');
        return
    else
        rem=w;
        while ~isempty(rem),
            [token,rem]=strtok(rem,char(10));
            if ~isempty(token),
                add_msg_board(token);
            end;
        end;
        add_msg_board('Now importing aligned sequences.');
        alignment=get_multiple_fasta(outfile);
    end;
end;


fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='Alignment file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'C; Alignment file written by MMM based on %s alignment\n\n',type);

names=':';
codes=':';
for k=1:length(selection),
    alg=alignment(selection(k));
    if strcmpi(alg.db,'Modeller'),
       seqtype=alg.seqtype;
       name=alg.name;
       type=alg.type;
       idCode=alg.pdb;
       first_res=str2double(alg.first_res);
       first_chain=alg.first_chain;
       last_res=str2double(alg.last_res);
       last_chain=alg.last_chain;
       protein=alg.protein;
       source=alg.source;
       resolution=alg.resolution;
       Rfactor=alg.Rfactor;
    else
        seqtype='P1';
        name=strtok(alg.name,' :');
        if strcmpi(alg.db,'pdb'),
            type='structure';
            idCode=alg.name;
            first_chain=alg.chain;
            last_chain=alg.chain;
        else
            type='sequence';
            idCode=strtok(alg.name,' :'); % use only part until first space
            first_chain=' ';
            last_chain=' ';
        end;
        first_res=1;
        poi=0;
        last_res=0;
        while poi<length(alg.sequence),
            poi=poi+1;
            if char(alg.sequence(poi))~='-',
                last_res=last_res+1;
            end;
        end;
        protein='unknown';
        source='unknown';
        resolution='-1.00';
        Rfactor='-1.00';
    end;
    % make sure that codes and names are unique throughout file
    id=tag2id(name,names);
    while ~isempty(id),
        name=[name 'a'];
        id=tag2id(name,names);
    end;
    names=[names name ':'];
    id=tag2id(idCode,codes);
    while ~isempty(id),
        idCode=[idCode 'a'];
        id=tag2id(idCode,codes);
    end;
    codes=[codes idCode ':'];
    fprintf(fid,'>%s;%s\n',seqtype,name);
    fprintf(fid,'%s:%s:%i:%s:%i:%s:%s:%s:%s:%s\n',type,idCode,first_res,...
        first_chain,last_res,last_chain,protein,source,resolution,Rfactor);
    tsequence=alg.sequence;
    while length(tsequence)>70,
        fprintf(fid,'%s\n',tsequence(1:70));
        tsequence=tsequence(71:end);
    end;
    fprintf(fid,'%s*\n\n',tsequence);
end;
fclose(fid);
