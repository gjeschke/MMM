function message=wr_het_frame(fname,idCode,modnum,chain,stub)
% function wr_het_frame(fname,idCode,modnum,chain)
%
% Writes a pseudo-PDB file in format V 3.20 of the current structure
% with only  
%   header, (HEADER)
%   title, (TITLE)
%   remarks on format and the originating program (REMARK 4, REMARK 5)
%   unit cell information (CRYST1), if present in the original file
%   atom coordinates of only hetero atoms (HETATM)
%   for hetero amino acids the backbone atoms are omitted
%   no ANISOU records are written
%   no CONNECT recors and MASTER are written
% this can be used as a "frame" file for SCWRL4 computations
%
% fname     file name for output, extension .pdb is appended, if no
%           extension is present
% idCode    four-letter PDB code or pseudo-PDB code
% modnum    optional argument, model (dataset) number for an ensemble
%           structure, if present, only this particular model (dataset) is
%           written
%           the argument should exist if the output is intended for SCWRL4
%           and the structure has more than one model
%           if the requested model does not exist, the whole structure is
%           written
%           if there are several chains, the requested model must exist in
%           all of them (or in none)
% chain     optional argument, chain number for a structure with multiple
%           chains, if present, only this particular chain is written
%           if the requested chain does not exist, no file is written
%           if chain=':' all chains are written
% stub      optional flag, if present and true, the finishing END is
%           omitted, defaults to false
%
% G. Jeschke, 2010

global model
global chemistry


message.error=0;
message.text='';

if nargin<5,
    stub=false;
end;

% initialize checksum values

numRemark=0; % number of REMARK records
numCoord=0;
numTer=0;


% append proper extension, if extension is missing
if isempty(strfind(fname,'.')), fname=strcat(fname,'.pdb'); end;

snum=model.current_structure;

% generate header line
header=sprintf('HEADER    %s',model.info{snum}.class);
header=fillstr(header,50);
today=date;
today=[today(1:7) today(10:11)];
header=sprintf('%s%s   %s',header,today,idCode);
% state supported format and originating program
format=sprintf('REMARK   4 %s COMPLIES WITH FORMAT V. 3.20, 01-DEC-08',idCode);
numRemark=numRemark+1;
origin=sprintf('REMARK   5 WRITTEN BY MMM (MODIFIED STRUCTURE)');
numRemark=numRemark+1;

fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='File could not be written';
    return;
end;

fprintf(fid,'%s\n',header);
[tl,n]=size(model.info{snum}.title);
for k=1:tl, % loop over all title lines
    if k==1,
        fprintf(fid,'TITLE     %s\n',model.info{snum}.title(k,:));
    else
        fprintf(fid,'TITLE   %2i%s\n',k,model.info{snum}.title(k,:));
    end;
end;
fprintf(fid,'%s\n',format);
fprintf(fid,'%s\n',origin);

chains=length(model.structures{snum}(:));

if nargin>=4 && ~strcmp(chain,':'),
    if chain<=chains && chain>=1,
        cnumbers=chain;
    else
        message.error=2;
        message.text='Chain number does not exist';
        return;
    end;
else
    cnumbers=1:chains;
end;

% write CRYST1 record, if unit cell information is present

if isfield(model.info{snum},'cryst'),
    cryst=model.info{snum}.cryst;
    fprintf(fid,'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s%4i\n',cryst.a,cryst.b,cryst.c,cryst.alpha,cryst.beta,cryst.gamma,cryst.sGroup,cryst.Z);
end;

% write ATOM and HETATM records

models=length(model.structures{snum}(1).residues);
modvec=1:models;
if nargin>2 && modnum<=models && modnum>=1,
    modvec=modnum;
end;

atoms=0;
for cnum=cnumbers,
    [catoms,n]=size(model.structures{snum}(cnum).xyz{1});
    atoms=atoms+catoms+1;
end;
atom_table=zeros(atoms,3); % preallocate atom number translation table
for mnum=modvec,
    atnum=0; % initialize atom serial number
    if models>1,
        fprintf(fid,'MODEL     %4i\n',mnum);
    end;
    for cnum=cnumbers,
        cid=model.structures{snum}(cnum).name;
        residues=length(model.structures{snum}(cnum).residues{1}.info);
        for rnum=1:residues,
            isaminoacid=model.structures{snum}(cnum).residues{1}.info(rnum).type==1;
            ishet=model.structures{snum}(cnum).residues{1}.info(rnum).hetflag;
            if isaminoacid && ~ishet, % don't write ATOM records for native amino acids
                continue;
            end;
            rtag=model.structures{snum}(cnum).residues{1}.info(rnum).name;
            backbone=false;
            if isaminoacid && ishet, % don't write backbone atoms for hetero amino acids
                atags=model.structures{snum}(cnum).residues{1}.info(rnum).atom_tags;
                id=tag2id('N',atags);
                if ~isempty(id), backbone=true; end;
                id=tag2id('O',atags);
                if ~isempty(id), backbone=true; end;
                id=tag2id('CA',atags);
                if ~isempty(id), backbone=true; end;
                id=tag2id('C',atags);
                if ~isempty(id), backbone=true; end;
            end;
            if backbone,
                continue;
            end;
            rid=sprintf('%4i ',model.structures{snum}(cnum).residues{1}.info(rnum).number);
            rid0=id2tag(rnum,model.structures{snum}(cnum).residues{1}.residue_tags);
            if double(rid0(end))<double('0') || double(rid0(end))>double('9'),
                rid(end)=rid0(end);
            end;
            tline='HETATM'; % everything here is a HETATM record, also nucleic acids (don't know how SCWRL4 treats them)
            hetflag=1;
            pointers=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers;
            if isfield(model.structures{snum}(cnum).residues{mnum}.info(rnum),'location_tags')
                ltags=model.structures{snum}(cnum).residues{mnum}.info(rnum).location_tags;
            else
                ltags=': :';
            end;
            for anum=1:length(pointers), % loop over atoms
                pointer=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum};
                atag=id2tag(anum,model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_tags);
                elnum=model.structures{snum}(cnum).residues{mnum}.info(rnum).elements(anum);
                element=upper(id2tag(elnum,chemistry.element_tags));
                if length(element)<2 && length(atag)<4,
                    atag=[' ' atag];
                end;
                if length(atag)<4,
                    atag=fillstr(atag,4);
                end;
                [loc,n]=size(pointer);
                if loc>26, loc=26; end;
                for lnum=1:loc, % loop over locations
                    poi=pointer(lnum,1); % actual coordinate set number
                    ltag='';
                    if n>1,
                        occupancy=pointer(lnum,2);
                        ltag=id2tag(pointer(lnum,3),ltags);
                    else
                        occupancy=1.00;
                    end;
                    if isempty(ltag), ltag=' '; end;
                    xyz=model.structures{snum}(cnum).xyz{mnum}(poi,:);
                    Bfactor=model.structures{snum}(cnum).Bfactor{mnum}(poi);
                    atnum=atnum+1;
                    atom_table(atnum,1)=poi;
                    atom_table(atnum,2)=cnum;
                    atom_table(atnum,3)=hetflag;
                    fprintf(fid,'%s%5i %4s%s%s %s%4s   %8.3f%8.3f%8.3f',tline,atnum,atag,ltag,rtag,cid,rid,xyz);
                    fprintf(fid,'%6.2f%6.2f          %2s\n',occupancy,Bfactor,element);
                    numCoord=numCoord+1;
                end;
            end;
            if isfield(model.structures{snum}(cnum).residues{mnum}.info(rnum),'terminal')
                if model.structures{snum}(cnum).residues{mnum}.info(rnum).terminal==1,
                    atnum=atnum+1;
                    fprintf(fid,'TER   %5i      %3s %s%4s\n',atnum,rtag,cid,rid);
                    numTer=numTer+1;
                end;
            end;
        end;
    end;
    if models>1,
        fprintf(fid,'ENDMDL\n');
    end;
end;

if ~stub
    % don't forget to close it politely
    fprintf(fid,'END   \n');
end;
fclose(fid);

if numCoord==0,
    message.error=2;
    message.text='No heteroatoms found. Frame file will be ommitted.';
end;

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');
% newstring = char(padarray(uint8(string)', newlength-length(string), 32,'post')');
