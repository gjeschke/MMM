function [snum,infile]=repacked_copy(snum0,idCode,modnum,chain,crystal)
% Makes a copy of an existing structure with repacking of the sidechains by
% SCWRL4
%
% snum0     number of structure to be copied
% idCode    internal structure identifier in MMM, should have four
%           characters as a PDB identifier, but start with a character
%           optional, automatic identifier is created if needed
%           if modnum is given and idCode is empty, the idCode of the
%           existing structure is used
% modnum    optional number of the model in an ensemble, defaults to 1, in
%           ensemble mode, the original structure with number snum0 is
%           overwritten, use copy_structure before, if you wish a copy;
%           if modnum is not provided, a new structure snum is added to the
%           model, which is a copy of structure snum0 with repacked
%           sidechains
% chain     optional number of the chain, if present, only this chain is
%           repacked, if requested chain does not exist, output arguments
%           are empty and no copy is made, note that the original structure
%           is overwritten, use copy_structure before, if you wish a copy
%           all chains are repacked for chain=0 or chain=[]
% crystal   optional flag, false (0) no inclusion of crystal environmen,
%           true (1) crystal environment is included, defaults to false
%
% snum      number of copied structure with repacked sidechains
% infile    file name of PDB output file of SCWRL4 (in temporary folder)
%
% Note: the calling routine is responsible for the model with number modnum
%       to exist in structure number snum0 if modnum is provided
%
% G. Jeschke, 2010


global model
global general
global help_files
global third_party
global hMain

chains=length(model.structures{snum0}(:));

if nargin<5,
    crystal=false;
end;

if nargin>=4 && ~isempty(chain) && chain~=0,
    if chain<=chains && chain>=1,
        cnumbers=chain;
    else
        snum=[];
        infile='';
        return;
    end;
else
    cnumbers=1:chains;
end;

old_snum=model.current_structure;
model.current_structure=snum0;

entry=strcat(help_files,'third_party.html#SCWRL4');

dospath=which('scwrl4.exe');
if isempty(dospath),
    add_msg_board('Sidechain repacking requires SCWRL4 by Krivov, Shapovalov, and Dunbrack');
    add_msg_board('ERROR: SCWRL4 could not be found on the Matlab path.');
    add_msg_board('Please check whether SCWRL4 is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    snum=[];
    infile='';
    return
end;

add_msg_board('Saving PDB input file...');
set(gcf,'Pointer','watch');
drawnow;

if nargin<2 || isempty(idCode),
    idCode=model.info{snum0}.idCode;
    if isempty(idCode), idCode='AMMM'; end;
    idCode(1)=char(idCode(1)+16);
end;

if nargin<3,
    modnum=1;
end;

outfile=[general.tmp_files 'rctmp.pdb'];
message=wr_pdb(outfile,idCode,modnum);

if message.error,
    add_msg_board(message.text);
end;

add_msg_board('Repacking side groups...');
drawnow

infile=[general.tmp_files 'rcscwrl4_out.pdb'];

cmd=[dospath ' -i ' outfile ' -o ' infile ' -h -t'];
if crystal,
    cmd=[cmd ' -#'];
end;
[s, w] = dos(cmd);
if s~=0,
    add_msg_board('ERROR: SCWRL4 did not run successfully.');
    set(gcf,'Pointer','arrow');
    return
else
    rem=w;
    while ~isempty(rem),
        [token,rem]=strtok(rem,char(10));
        if ~isempty(token),
            add_msg_board(token);
        end;
    end;
    add_msg_board('Now importing structure with repacked sidegroups...');
end;

if nargin<3,
    [message,snum]=add_pdb(infile,idCode);
    if message.error,
        add_msg_board('ERROR: Reading of SCWRL4 PDB output file failed.');
        add_msg_board(message.text);
    else
        add_msg_board(sprintf('Structure %s added to model as structure number %i',idCode,snum));
    end;
else
    [structure,info]=rd_pdb(infile,strcat(general.tmp_files,'rcremarks.txt'));
    for cnum=cnumbers,
        model.structures{snum0}(cnum).residues{modnum}=structure(cnum).residues{1};
        model.structures{snum0}(cnum).xyz{modnum}=structure(cnum).xyz{1};
        model.structures{snum0}(cnum).atoms{modnum}=structure(cnum).atoms{1};
        model.structures{snum0}(cnum).conn=structure(cnum).conn;
        model.structures{snum0}(cnum).Bfactor{modnum}=structure(cnum).Bfactor{1};
        model.structures{snum0}(cnum).Btensor{modnum}=structure(cnum).Btensor{1};
    end;
end;
set(gcf,'Pointer','arrow');
add_msg_board('Sidegroup repacking is complete.');

model.current_structure=old_snum;
stag=id2tag(old_snum,model.structure_tags);
if isfield(model.info{snum},'resolution') && ~isempty(model.info{snum}.resolution),
    resstring=sprintf('%4.2f Å',model.info{snum}.resolution);
else
    resstring='not specified';
end;
    
set(hMain.MMM,'Name',sprintf('MMM - [%s](%s) Resolution %s',stag,model.current_chain,resstring));

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

