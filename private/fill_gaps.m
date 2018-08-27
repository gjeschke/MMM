function message=fill_gaps(indices)
% function message=fill_gaps(indices)
%    Fills gaps in a chain structure using Modeller,
%    modeller must be installed
%
% indices   1x2 index vector specifying structure and chain
%
% message   error message with fields message.error (integer, 0 no error)
%           and message.text (error description)
%
% Limitations:
% - in an ensemble model, this can take very long
% - missing residues at the beginning and end of the chain are not
%   reconstructed
%
% G. Jeschke, 2011


global model
global general
global help_files
global third_party
global residue_defs
global hMain

message.error=0;
message.text='No error.';

[stag,ctag]=mk_address_parts(indices);
[fname,gaplist,message]=mk_pir('fill','tmp',indices);

if isempty(gaplist),
    add_msg_board(sprintf('No gaps of appropriate length in chain [%s](%s). Nothing modified.',stag,ctag));
    message.error=1;
    message.text='No gaps of appropriate length in structure.';
    return;
end;


entry=strcat(help_files,'third_party.html#Modeller');

dospath=which([third_party.modeller_version  '.exe']);
if isempty(dospath),
    message.error=2;
    message.text='Modeller software not found on Matlab path.';
    add_msg_board('Filling gaps requires Modeller from the Sali lab');
    add_msg_board('ERROR: Modeller could not be found on the Matlab path.');
    add_msg_board('Please check whether Modeller is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;
[modpath, modcmd] = fileparts(dospath);

snum0=model.current_structure;
snum=indices(1);
model.current_structure=snum;
cnum=indices(2);
modnum=length(model.structures{snum}(cnum).xyz);
if modnum>1,
    add_msg_board(sprintf('Ensemble structure with %i models. This can take long.',modnum));
end;

for mnum=1:modnum,

    offset = model.structures{snum}(cnum).residues{mnum}.info(1).number;

    add_msg_board('Saving PDB input file...');
    set(gcf,'Pointer','watch');
    drawnow;

    outfile=[general.tmp_files 'tmp.pdb'];
    message=wr_pdb(outfile,'MMMM',mnum,cnum);
    set(gcf,'Pointer','watch');
    drawnow;

    if message.error,
        add_msg_board(message.text);
    end;

    add_msg_board(sprintf('Modeling missing residues in model number %i...',mnum));
    adr=mk_address([indices mnum]);
    cmd(hMain,sprintf('hide %s',adr));
    set(gcf,'Pointer','watch');
    drawnow

    infile=[general.tmp_files 'run_modeller.py'];
    
    message=mk_modeller_input(infile,'fill','tmp',gaplist,offset);
    if message.error,
        continue;
    end;
    
    [batcmd,message]=mk_modeller_bat(modpath,modcmd,infile);
    if message.error,
        add_msg_board('ERROR: Modeller could not be initialized');
    end;
    my_dir=pwd;
    cd(modpath);
    [s, w] = dos(batcmd);
    if s~=0,
        rem=w;
        while ~isempty(rem),
            [token,rem]=strtok(rem,char(10));
            if ~isempty(token),
                add_msg_board(token);
            end;
        end;
        message.error=2;
        message.text='Modeller error.';
        add_msg_board('ERROR: Modeller did not run successfully.');
        set(gcf,'Pointer','arrow');
        cd(my_dir);
        return
    else
        analyze_modeller_log('run_modeller.log',true);
        add_msg_board('Now importing remodeled structure...');
    end;
    cd(my_dir);

    infile='fill.B99990001.pdb';
    structure=rd_pdb(infile,strcat(general.tmp_files,'remarks.txt'));
    set(gcf,'Pointer','arrow');
    [mg,ng]=size(gaplist);
    poi=0;
    for k=1:mg,
        if gaplist(k,2)-gaplist(k,1)<=7,
            poi=poi+1;
            gaplist(poi,:)=gaplist(k,:);
        elseif gaplist(k,2)-gaplist(k,1)<=9,
            add_msg_board(sprintf('Warning: Modeled loop from residue %i to %i may be unreliable.',gaplist(k,1),gaplist(k,2))); 
            poi=poi+1;
            gaplist(poi,:)=gaplist(k,:);
        else
            button = questdlg('Do you want to insert this loop to create an uncertain model?',sprintf('At length of %i loop model is possibly unreliable.',...
                gaplist(k,2)-gaplist(k,1)+1),'No','Yes','No');
            if strcmpi(button,'Yes'),
                poi=poi+1;
                gaplist(poi,:)=gaplist(k,:);
            end;
        end;                
    end;
    gaplist=gaplist(1:poi,:);
    set(gcf,'Pointer','watch');
    insert_residues(gaplist,structure,snum,cnum,mnum,offset);
    add_msg_board(sprintf('Gaps repaired in model %i.',mnum));

end;
message=consolidate_chain([snum,cnum]);
adr=mk_address([snum,cnum]);
cmd(hMain,sprintf('show %s ribbon',adr));

set(gcf,'Pointer','arrow');
add_msg_board('Gap repairing is complete.');

drawnow;

% add the reference, if it does not yet exist
modeller_ref=true;
id=tag2id('Eswar:2008_modeller',third_party.tags,[],'|');
if ~isempty(id),
    if isfield(model,'auto_references'),
        if ~isempty(find(id==model.auto_references, 1)),
            modeller_ref=false;
        end;
    else
        model.auto_references=[];
    end;
    if modeller_ref,
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

        
function message=mk_modeller_input(infile,basname,tempname,gaplist,offset)

message.error=0;
message.text='No error.';


if nargin < 5,
    offset = 0; 
end;

fid=fopen(infile,'wt');
if fid==-1,
    message.error=2;
    message.text='Modeller input file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

[m,n]=size(gaplist);

fprintf(fid,'from modeller import *\n');
fprintf(fid,'from modeller.automodel import *\n\n');
fprintf(fid,'log.verbose()\n');
fprintf(fid,'env = environ()\n\n');

fprintf(fid,'env = environ(rand_seed=-8123, restyp_lib_file=''$(LIB)/restyp_EPR.lib'', copy=None)\n\n');

fprintf(fid,'env.io.atom_files_directory = [''.'', ''../atomfiles'']\n\n');

fprintf(fid,'env.io.hetatm = True\n\n');

fprintf(fid,'class MyModel(automodel):\n');
fprintf(fid,'    def select_atoms(self):\n');
fprintf(fid,'        return selection(');
for k=1:m,
    fprintf(fid,'self.residue_range(''%i'', ''%i'')',gaplist(k,1) - offset,gaplist(k,2) - offset);
    if k<m,
        fprintf(fid,',\n        ');
    else
        fprintf(fid,')\n\n');
    end;
end;

fprintf(fid,'a = MyModel(env, alnfile = ''%s.ali'',\n',basname);
fprintf(fid,'            knowns = ''%s'', sequence = ''%s'')\n',tempname,basname);
fprintf(fid,'a.starting_model= 1\n');
fprintf(fid,'a.ending_model  = 1\n\n');

fprintf(fid,'a.make()\n');

fclose(fid);

function [batname,message]=mk_modeller_bat(modpath,modcmd,jobfile)

global general

batname='run_modeller.bat';
envfile=[modpath filesep 'modenv.bat'];
batfile=[modpath filesep  batname];

message.error=0;
message.text='No error.';

fid=fopen(batfile,'wt');
if fid==-1,
    message.error=2;
    message.text='Modeller batch file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

rfid=fopen(envfile,'rt');
if rfid==-1,
    message.error=3;
    message.text='Modeller environment file could not be read.';
    add_msg_board('ERROR: Modeller environment file could not be read.');
    add_msg_board(message.text);
    return;
end;

while 1
    tline = fgetl(rfid);
    if ~ischar(tline), break, end
    fprintf(fid,'%s\n',tline);
end;

fprintf(fid,'echo modeller installed\n');
k=strfind(general.tmp_files,':');
if ~isempty(k),
    fprintf(fid,'%s\n',general.tmp_files(1:k));
end;
fprintf(fid,'cd %s\n',general.tmp_files);
fprintf(fid,'%s %s\n',modcmd,jobfile);

fclose(fid);

function insert_residues(gaplist,structure,snum,cnum,mnum,offset)

% Bfactor is set to twice the maximum Bfactor of existing residues

global model
global residue_defs

if nargin < 6,
    offset = 0;
end;

[m,n]=size(gaplist);

if isempty(gaplist),
    return
end;

Bmax=max(model.info{snum}.B_range);

newresidues=gaplist(1,1):gaplist(1,2);
for k=2:m,
    newresidues=[newresidues gaplist(k,1):gaplist(k,2)];
end;
[ma,na]=size(model.structures{snum}(cnum).xyz{mnum});
[ma0,maxconn]=size(model.structures{snum}(cnum).conn);
newiso=[model.structures{snum}(cnum).isotopes; zeros(20*length(newresidues),2,'single')];
newxyz=[model.structures{snum}(cnum).xyz{mnum}; zeros(20*length(newresidues),3)];
newBfac=[model.structures{snum}(cnum).Bfactor{mnum}, zeros(1,20*length(newresidues))];
newBtens=[model.structures{snum}(cnum).Btensor{mnum}; zeros(20*length(newresidues),6,'int32')];
sequence=model.structures{snum}(cnum).sequence;
restags=model.structures{snum}(cnum).residues{mnum}.residue_tags;

newrnum=zeros(length(newresidues));
nr=length(model.structures{snum}(cnum).residues{mnum}.info);
rpoi=0;
for k=1:length(newresidues),
    rnum0=newresidues(k);
    rnum=0;
    for kk=1:length(structure(1).residues{1}.info),
        if structure(1).residues{1}.info(kk).number == rnum0 - offset + 1,
            rnum=kk;
            break;
        end;
    end;
    if rnum==0,
        add_msg_board(sprintf('Warning: Residue %i not found in Modeller result file',rnum0));
        continue;
    end;
    tag=structure(1).residues{1}.info(rnum).name;
    id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
    sequence(rnum0)=id;
    if mnum==1 && isfield(model.structures{snum}(cnum),'seqexist'),
        model.structures{snum}(cnum).seqexist(rnum0)=1;
    end;
    nr=nr+1;
    rpoi=rpoi+1;
    newrnum(rpoi)=nr;
    % model.structures{snum}(cnum).residues{mnum}.info(nr)=structure(1).residues{1}.info(rnum);
    restag=sprintf('%i:',rnum0);
    restags=strcat(restags,restag);
    pointers=structure(1).residues{1}.info(rnum).atom_numbers;
    for anum=1:length(pointers), % loop over atoms
        pointer=structure(1).residues{1}.info(rnum).atom_numbers{anum};
        [loc,n]=size(pointer);
        for lnum=1:loc, % loop over locations
            poi=pointer(lnum,1); % actual coordinate set number
            ma=ma+1;
            pointer(lnum,1)=ma;
            newiso(ma,:)=structure(1).isotopes(poi,:);
            newxyz(ma,:)=structure(1).xyz{1}(poi,:);
            newBfac(ma)=2*Bmax;
            newBtens(ma,:)=structure(1).Btensor{1}(poi,:);
        end;
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_numbers{anum}=pointer;
        model.structures{snum}(cnum).residues{mnum}.info(nr).name=tag;
        model.structures{snum}(cnum).residues{mnum}.info(nr).type=structure(1).residues{1}.info(rnum).type;
        model.structures{snum}(cnum).residues{mnum}.info(nr).secondary=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).hetflag=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).connected=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).number=structure(1).residues{1}.info(rnum).number + offset -1;
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_tags=structure(1).residues{1}.info(rnum).atom_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).elements=structure(1).residues{1}.info(rnum).elements;
        model.structures{snum}(cnum).residues{mnum}.info(nr).location_tags=structure(1).residues{1}.info(rnum).location_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).insertion_code=structure(1).residues{1}.info(rnum).insertion_code;
    end;
end;
newiso=newiso(1:ma,:);
newxyz=newxyz(1:ma,:);
newBfac=newBfac(1:ma);
newBtens=newBtens(1:ma,:);

if mnum==1,
    model.structures{snum}(cnum).sequence=sequence;
    model.structures{snum}(cnum).isotopes=newiso;
    model.structures{snum}(cnum).conn=[model.structures{snum}(cnum).conn; zeros(ma-ma0,maxconn)];
end;

model.structures{snum}(cnum).xyz{mnum}=newxyz;
model.structures{snum}(cnum).Bfactor{mnum}=newBfac;
model.structures{snum}(cnum).Btensor{mnum}=newBtens;
model.structures{snum}(cnum).residues{mnum}.residue_tags=restags;

% make internal bonds in new residues
for k=1:length(newrnum),
    model.structures{snum}(cnum)=mk_internal_bonds(model.structures{snum}(cnum),newrnum(k),residue_defs);
end;

% sort residues by number
info=model.structures{snum}(cnum).residues{mnum}.info;
numbers=length(info);
for k=1:length(info),
    numbers(k)=info(k).number;
end;
info0=info;
[sorted,oldnumbers]=sort(numbers);
tags=':';
for k=1:length(oldnumbers),
    info(k)=info0(oldnumbers(k));
    tag=id2tag(oldnumbers(k),model.structures{snum}(cnum).residues{mnum}.residue_tags);
    tags=[tags tag ':'];
end;
model.structures{snum}(cnum).residues{mnum}.residue_tags=tags;
model.structures{snum}(cnum).residues{mnum}.info=info;

