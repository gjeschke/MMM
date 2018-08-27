function message=repack(snum,modnum,crystal,sequence)
% function message=repack(snum,modnum,crystal,sequence)
%    Repacks the sidechains of all native residues in a structure
%    using third-party software SCWRL4 (which must be installed)
%
% snum      number of the structure to be repacked
% modnum    optional number of the model in an ensemble, defaults to 1, in
%           ensemble mode, the original structure with number snum0 is
%           overwritten, use copy_structure before, if you wish a copy;
%           if modnum is not provided, a new structure snum is added to the
%           model, which is a copy of structure snum0 with repacked
%           sidechains
% crystal   optional flag, false (0) no inclusion of crystal environment,
%           true (1) crystal environment is included, defaults to false
%
% Limitations:
%
% in a multi-chain complex, invariably all chains are repacked
% residues of type UNK are treated as glycine
%
% G. Jeschke, 2010


global model
global general
global help_files
global third_party
global residue_defs
global hMain

if nargin<3,
    crystal=false;
end;

% if nargin>=4,
%     rem=sequence;
%     k=0;
%     while ~isempty(rem),
%         k=k+1;
%         [seq,rem]=strtok(rem,'/');
%         sequences{k}=seq;
%     end;
% else
%     sequences={};
% end;

standard='ARNDCEQGHILKMFPSTWYV'; % single letter codes of all standard amino acids

message.error=0;
message.text='No error.';

entry=strcat(help_files,'third_party.html#SCWRL4');

dospath=which('scwrl4.exe');
if isempty(dospath),
    message.error=5;
    message.text='SCWRL4 software not found on Matlab path.';
    add_msg_board('Sidechain repacking requires SCWRL4 by Krivov, Shapovalov, and Dunbrack');
    add_msg_board('ERROR: SCWRL4 could not be found on the Matlab path.');
    add_msg_board('Please check whether SCWRL4 is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;

chains=length(model.structures{snum}(:));
longseq='';
indices=zeros(10000,4);
repacked=zeros(1,10000);
ranges=zeros(chains,2);
poi=0;
for c=1:chains,
    residues=length(model.structures{snum}(c).residues{modnum}.info);
    newr=0;
    if c>1,
        ranges(c-1,2)=length(longseq);
    end;
    ranges(c,1)=length(longseq)+1;
    for r=1:residues,
        rtag=model.structures{snum}(c).residues{modnum}.info(r).name;
        atags=model.structures{snum}(c).residues{modnum}.info(r).atom_tags;
        backbone=true;
        id=tag2id('N',atags);
        if isempty(id), backbone=false; end;
        id=tag2id('O',atags);
        if isempty(id), backbone=false; end;
        id=tag2id('CA',atags);
        if isempty(id), backbone=false; end;
        id=tag2id('C',atags);
        if isempty(id), backbone=false; end;
        if backbone,
            newr=newr+1;
            cslc=tag2id(upper(rtag),upper(residue_defs.restags),residue_defs.single_letter_code);
            if isempty(cslc),
                cslc='g';
            else
                cslc=upper(cslc);
                poi=poi+1;
                indices(poi,1)=snum;
                indices(poi,2)=c;
                indices(poi,3)=modnum;
                indices(poi,4)=r;
                repacked(poi)=newr;
            end;
            longseq=[longseq cslc];
        end;
    end;
end;
ranges(chains,2)=length(longseq);

if poi==0,
    message.error=5;
    message.text='No residues found that could be repacked.';
    add_msg_board('ERROR: No repacking performed.');
    add_msg_board(message.text);
    return
end;

if nargin>3,
    if length(sequence)~=length(longseq),
        message.error=6;
        message.text='Target sequence length does not match sequence length in structure.';
        add_msg_board('ERROR: No repacking performed.');
        add_msg_board(message.text);
        return
    else
        longseq=sequence;
    end
end;

indices=indices(1:poi,:);
repacked=repacked(1:poi);

seqfile=[general.tmp_files 'rtmp.seq'];
fid=fopen(seqfile,'wt');
if fid==-1,
    message.error=7;
    message.text='Sequence file could not be written.';
    add_msg_board('ERROR: No repacking performed.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'%s\n',longseq);
fclose(fid);


add_msg_board('Saving PDB input file...');
set(gcf,'Pointer','watch');
drawnow;

outfile=[general.tmp_files 'rtmp.pdb'];
message=wr_pdb(outfile,'MMMM',modnum,':',true);

if message.error,
    add_msg_board(message.text);
end;

framefile=[general.tmp_files 'rframe.pdb'];
add_msg_board('Saving PDB frame file (heteroatoms)...');
message=wr_het_frame(framefile,'MHET',modnum,':');
if message.error,
    frame=false;
else
    frame=true;
end;

if message.error,
    add_msg_board(message.text);
end;

add_msg_board('Now packing residues...');
drawnow

infile=[general.tmp_files 'rscwrl4_out.pdb'];

comd=[dospath ' -i ' outfile ' -o ' infile];
if frame,
    comd=[comd ' -f ' framefile];
end;
comd=[comd ' -s ' seqfile ' -h -t'];
if crystal,
    comd=[comd ' -#'];
end;
[s, w] = dos(comd);
scwrl4_error=false;
if s~=0,
    message.error=8;
    message.text='SCWRL4 error.';
    add_msg_board('ERROR: SCWRL4 did not run successfully.');
    set(gcf,'Pointer','arrow');
    return
else
    rem=w;
    while ~isempty(rem),
        [token,rem]=strtok(rem,char(10));
        if ~isempty(token),
            add_msg_board(token);
            if strfind(token,'Error'),
                scwrl4_error=true;
            end;
        end;
    end;
    if scwrl4_error,
        message.error=8;
        message.text='SCWRL4 error.';
        add_msg_board('ERROR: SCWRL4 did not run successfully.');
        set(gcf,'Pointer','arrow');
        return
    end;
    add_msg_board('Now importing structure with repacked residues...');
end;

[mres,nres]=size(indices);
add_msg_board(sprintf('Replacing sidechains of %i residues.',mres));
[structure,info]=rd_pdb(infile,strcat(general.tmp_files,'rremarks.txt'));

for r=1:mres,
    adr=mk_address(indices(r,:),1);
    add_msg_board(sprintf('Replacing sidechain of residue %s',adr));
    cmd(hMain,sprintf('hide %s',adr));
    replace_sidechain(structure,indices(r,:),repacked(r));
end;

for c=1:chains,
    seq=longseq(ranges(c,1):ranges(c,2));
    seq0=model.structures{snum}(c).sequence;
    poi=1;
    while strcmp(seq0(poi),'?'),
        seq=['?' seq];
        poi=poi+1;
    end;
    model.structures{snum}(c).sequence=seq;
end;

set(gcf,'Pointer','arrow');
add_msg_board('Repacking is complete.');

drawnow;

% add the reference, if it does not yet exist
scwrl4_ref=true;
id=tag2id('Krivov:2009_scwrl4',third_party.tags,[],'|');
if ~isempty(id),
    if isfield(model,'auto_references'),
        if ~isempty(find(id==model.auto_references, 1)),
            scwrl4_ref=false;
        end;
    else
        model.auto_references=[];
    end;
    if scwrl4_ref,
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

