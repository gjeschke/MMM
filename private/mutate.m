function message=mutate(adr,type)
% function message=mutate(adr,type)
%    Mutates a residue to another amino acid, can also be used for
%    sidechain repacking of single residues, SCWRL4 must be installed
%
% adr       residue address, must address a single residue, can also be a
%           1x4 index vector for a residue
% type      optional type of amino acid after mutation, can be a
%           three-letter or single-letter code, must denote one of the 20
%           standard amino acids, if argument is missing, the sidechain is
%           not changed, but repacked
%           no mutation or repacking is performed if the residue after
%           this would not be a standrad amino acid
% message   error message with fields message.error (integer, 0 no error)
%           and message.text (error description)
%
% Limitations:
% in an ensemble model, adr addresses the residue only in one model
% (coordinate set), to ensure consistency, the residue must be mutated in
% all coordinate sets (they must have the same sequence), this is done here
% (and can take very long)
%
% in a multi-chain protein complex all chains must have the same number of
% models (coordinate sets)
%
% G. Jeschke, 2010


global model
global general
global help_files
global third_party
global residue_defs
global hMain

standard='ARNDCEQGHILKMFPSTWYV'; % single letter codes of all standard amino acids

mutnum=0;
message.error=0;
message.text='No error.';

if isa(adr,'char'),
    [indices,message]=resolve_address(adr);
    if message.error,
        return;
    end;
else
    indices=adr;
end;

[m,n]=size(indices);
if m==1,
    indices=indices(indices>0);
    [m,n]=size(indices);
end;
if m~=1 || n~=4,
    message.error=1;
    message.text='Addressed object is not a single residue.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

if ~isa(type,'char'),
    message.error=2;
    message.text='Requested residue type must be single-letter or three-letter code.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

if length(type)==1,
    slc=type;
    id=tag2id(upper(type),upper(residue_defs.single_letter_code));
    tlc=id2tag(id,residue_defs.restags);
elseif length(type)==3,
    slc=tag2id(upper(type),upper(residue_defs.restags),residue_defs.single_letter_code);
    tlc=type;
else
    message.error=2;
    message.text='Requested residue type must be single-letter or three-letter code.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

if length(slc)~=1,
    message.error=3;
    message.text='Three-letter code does not translate into amino acid single-letter code.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

k=strfind(standard,slc);
if isempty(k),
    message.error=4;
    message.text='Single-letter code does not correspond to a well defined standard amino acid.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

entry=strcat(help_files,'third_party.html#SCWRL4');

dospath=which('scwrl4.exe');
if isempty(dospath),
    message.error=5;
    message.text='SCWRL4 software not found on Matlab path.';
    add_msg_board('Mutation requires SCWRL4 by Krivov, Shapovalov, and Dunbrack');
    add_msg_board('ERROR: SCWRL4 could not be found on the Matlab path.');
    add_msg_board('Please check whether SCWRL4 is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;

snum0=model.current_structure;
snum=indices(1);
model.current_structure=snum;
cnum=indices(2);
rnum=indices(4);
modnum=length(model.structures{snum}(cnum).xyz);
if modnum>1,
    add_msg_board(sprintf('Ensemble structure with %i models. This can take long.',modnum));
end;

chains=length(model.structures{snum}(:));
longseq='';
mutated=0;
for c=1:chains,
    residues=length(model.structures{snum}(c).residues{1}.info);
    newr=0;
    for r=1:residues,
        rtag=model.structures{snum}(c).residues{1}.info(r).name;
        atags=model.structures{snum}(c).residues{1}.info(r).atom_tags;
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
                cslc='G';
            end;
            if c==cnum && r==rnum,
                cslc=upper(slc); % the residue to be mutated is uppercase
                mutated=newr;
            else
                cslc=lower(cslc); % othe residues are lowercase
            end;
            longseq=[longseq cslc];
        end;
    end;
end;

if mutated==0,
    message.error=5;
    message.text='Residue to be mutated not found or has incomplete backbone.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

seqfile=[general.tmp_files 'tmp.seq'];
fid=fopen(seqfile,'wt');
if fid==-1,
    message.error=6;
    message.text='Sequence file could not be written.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'%s\n',longseq);
fclose(fid);


for mnum=1:modnum,

    add_msg_board('Saving PDB input file...');
    set(gcf,'Pointer','watch');
    drawnow;

    outfile=[general.tmp_files 'tmp.pdb'];
    message=wr_pdb(outfile,'MMMM',mnum,':',true);

    if message.error,
        add_msg_board(message.text);
    end;

    framefile=[general.tmp_files 'frame.pdb'];
    add_msg_board('Saving PDB frame file (heteroatoms)...');
    message=wr_het_frame(framefile,'MHET',mnum,':');
    if message.error,
        frame=false;
    else
        frame=true;
    end;
    
    if message.error,
        add_msg_board(message.text);
    end;
    
    message.error=0;
    message.text='No error.';

    add_msg_board(sprintf('Mutating and packing residue in model number %i...',mnum));
    drawnow

    infile=[general.tmp_files 'scwrl4_out.pdb'];

    comd=[dospath ' -i ' outfile ' -o ' infile];
    if frame,
        comd=[comd ' -f ' framefile];
    end;
    comd=[comd ' -s ' seqfile ' -h -t'];
    [s, w] = dos(comd);
    save(fullfile(general.tmp_files,'scwrl4_diagnostics.mat'),'s','w','comd');
    if s~=0,
        message.error=7;
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
                if strfind(token,'Total minimal energy of the graph'),
                    [comment,energy]=strtok(token,'=');
                    en=str2num(energy(2:end));
                    conv_factor=(4.1868*1000); % conversion to CI, - if working per mol
                    add_msg_board(sprintf('Energy is %10.1f kJ/mol\n',conv_factor*en/1000));
                end;
            end;
        end;
        add_msg_board('Now importing structure with mutated residue...');
    end;

    cindices=indices;
    cindices(3)=mnum;
    cindices=cindices(cindices>0);
    adr=mk_address(indices,1);
    add_msg_board(sprintf('Replacing sidechain of residue %s',adr));
    cmd(hMain,sprintf('hide %s',adr));

    [structure,info]=rd_pdb(infile,strcat(general.tmp_files,'remarks.txt'));
    
    set_sequence(indices,slc,tlc);
    replace_sidechain(structure,indices,mutated,mnum);
    
    add_msg_board(sprintf('Mutation performed in data set (model) %i.',mnum));

end;

set(gcf,'Pointer','arrow');
add_msg_board('Mutation is complete.');

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

model.current_structure=snum0;

        
function set_sequence(indices,slc,tlc)
    % updates single-letter code and three-letter code in the sequence record
global model

rnum2=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
seq=model.structures{indices(1)}(indices(2)).sequence;
if indices(4)>=1 && indices(4)<=length(seq),
    seq(rnum2)=slc;
    poi=2+4*(rnum2-1);
    model.structures{indices(1)}(indices(2)).restags(poi:poi+2)=upper(tlc);
end;
model.structures{indices(1)}(indices(2)).sequence=seq;
if isfield(model.structures{indices(1)}(indices(2)),'mutations'),
    mutations=model.structures{indices(1)}(indices(2)).mutations;
    nm=length(mutations);
    mutations(nm+1).number=rnum2;
    mutations(nm+1).original=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
    mutations(nm+1).modified=upper(tlc);
    mutations(nm+1).comment=upper('in silico mutation');
else
    mutations(1).number=rnum2;
    mutations(1).original=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
    mutations(1).modified=upper(tlc);
    mutations(1).comment=upper('in silico mutation');
end;
model.structures{indices(1)}(indices(2)).mutations=mutations;
