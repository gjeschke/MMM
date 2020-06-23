function [err_msg,atnum] = wr_residue_to_PDB(fid,adr,atnum,chain_id,resnum)
% write_residue_to_PDB(adr,chain_id,resnum)
%
% writes PDB file lines for a single residue
%
% fid       file identifier
% adr       MMM address of the residue
% atnum     atom number of last residue, defaults to 0
% chain_id  optional chain identifiere, defaults to the one specified by
%           adr
% resnum    optional residue number, defualts to the one specified by adr
%
% err_msg   error message, struct with fields .error and .msg 
%           .error  error code, 
%               0   no error
%               -1  file does not exist
%               1   address syntax error
%               2   no object with this address exists
%               3   adr does not specify a residue 
%           .text   clear text error message
%
% G. Jeschke, 23.6.2020

global model
global chemistry

if ~exist('atnum','var') || isempty(atnum)
    atnum = 0;
end

err_msg.error = -1;
err_msg.text = 'File does not exist';

if fid == -1
    return
end

[indices,err_msg]=resolve_address(adr);

if err_msg.error ~= 0
    return
end

if length(indices) ~= 4
    err_msg.error = 3;
    err_msg.tex = 'Address does not specify a residue';
    return
end

[~,ctag] = mk_address_parts(indices);

if ~exist('chain_id','var') || isempty(chain_id)
    chain_id = ctag;
end

err_msg.error = 0;
err_msg.text = 'No error';

snum = indices(1);
cnum = indices(2);
mnum = indices(3);
rnum = indices(4);
if ~exist('resnum','var') || isempty(resnum)
    resnum = model.structures{snum}(cnum).residues{1}.info(rnum).number;
end

rtag = model.structures{snum}(cnum).residues{1}.info(rnum).name;
rid = sprintf('%4i ',resnum);
rid0 = id2tag(rnum,model.structures{snum}(cnum).residues{1}.residue_tags);
if double(rid0(end)) < double('0') || double(rid0(end)) > double('9')
    rid(end)=rid0(end);
end
if model.structures{snum}(cnum).residues{1}.info(rnum).hetflag
    tline = 'HETATM';
else
    tline = 'ATOM  ';
end
pointers = model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers;
if isfield(model.structures{snum}(cnum).residues{mnum}.info(rnum),'location_tags')
    ltags = model.structures{snum}(cnum).residues{mnum}.info(rnum).location_tags;
else
    ltags = ': :';
end
for anum = 1:length(pointers) % loop over atoms
    pointer = model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum};
    atag = id2tag(anum,model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_tags);
    elnum = model.structures{snum}(cnum).residues{mnum}.info(rnum).elements(anum);
    element = upper(id2tag(elnum,chemistry.element_tags));
    if length(element)<2 && length(atag)<4
        atag = strcat(' ',atag);
    end
    if length(atag)<4
        atag = fillstr(atag,4);
    end
    [loc,n] = size(pointer);
    if loc>26 
        loc = 26; 
    end
    for lnum = 1:loc % loop over locations
        poi = pointer(lnum,1); % actual coordinate set number
        ltag='';
        if n > 1
            occupancy = pointer(lnum,2);
            ltag = id2tag(pointer(lnum,3),ltags);
        else
            occupancy = 1;
        end
        if isempty(ltag)
            ltag = ' '; 
        end
        xyz = model.structures{snum}(cnum).xyz{mnum}(poi,:);
        Bfactor = model.structures{snum}(cnum).Bfactor{mnum}(poi);
        atnum = atnum + 1;
        fprintf(fid,'%s%5i %4s%s%s %s%4s   %8.3f%8.3f%8.3f',tline,atnum,atag,ltag,rtag,chain_id,rid,xyz);
        fprintf(fid,'%6.2f%6.2f          %2s\n',occupancy,Bfactor,element);
        if sum(model.structures{snum}(cnum).Btensor{mnum}(poi,:)) >= 1 % write ANISOU record
            Btensor = model.structures{snum}(cnum).Btensor{mnum}(poi,:);
            fprintf(fid,'ANISOU%5i %4s%s%3s %s%4s %7i%7i%7i%7i%7i%7i',atnum,atag,ltag,rtag,cid,rid,Btensor);
            fprintf(fid,'      %2s\n',element);
        end
    end
end
if isfield(model.structures{snum}(cnum).residues{mnum}.info(rnum),'terminal')
    if model.structures{snum}(cnum).residues{mnum}.info(rnum).terminal == 1
        atnum = atnum+1;
        fprintf(fid,'TER   %5i      %3s %s%4s\n',atnum,rtag,chain_id,rid);
    end
end

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');