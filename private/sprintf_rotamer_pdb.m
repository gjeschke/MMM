function rotamer=sprintf_rotamer_pdb(rotamer_stats,rotnum,label,MD2PDB,backbone)
% function rotamer=sprintf_rotamer_pdb(rotamer_stats,rotnum,label,MD2PDB,backbone)
%
% writes a PDB file with all rotamers of a spin labeled residue 
% the file name is derived from the residue address adr, the label name 
% and the labeling temperature T
% 
% further input parameters:
% rotamer_stats     the rotamer statistics as computed by get_rotamers.m
% rotnum            number of the rotamer to be written out
% label             type of label
% MD2PDB            translation table for atom numbers
% backbone          optional flag to decide whether backbone atoms are
%                   written, defaults to true

global label_defs
global chemistry

if nargin<5,
    backbone=true;
end;

record='HETATM';
switch label
    case 'R1A'
        record='HETATM';
    case 'IA1'
        record='HETATM';
end;

rmsd=1; % standard atom r.m.s.d. in Å for B factor specification
Bfactor=8*pi^2*rmsd^2;

labnum=tag2id(label,label_defs.restags);
atags=label_defs.residues(labnum).atoms;

pop=rotamer_stats.all_potentials.pop_rot;
pop=ones(size(pop));

rotamer='';

atnum=0;
ecoor0=rotamer_stats.all_rotamers{rotnum};
ecoor=zeros(size(ecoor0));
[atoms,n]=size(ecoor);
for kk=1:length(MD2PDB),
    if MD2PDB(kk)>0 && MD2PDB(kk)<=atoms,
        ecoor(kk,:)=ecoor0(MD2PDB(kk),:);
    else
        ecoor(kk,:)=[NaN,NaN,NaN,NaN];
    end;
end;
rid=sprintf('%4i ',rotnum+9500);
pconn=0;
for kk=1:atoms,
    atag=id2tag(kk,atags);
    if ~strcmpi(atag,'H2') && ~strcmpi(atag,'OXT') && ~strcmpi(atag,'HXT') && ~sum(isnan(ecoor(kk,:))),
        atnum=atnum+1;
        pconn=pconn+1;
        xyz=ecoor(kk,2:4);
        elnum=ecoor(kk,1);
        element=upper(id2tag(elnum,chemistry.element_tags));
        if length(element)<2 && length(atag)<4,
            atag=[' ' atag];
        end;
        if length(atag)<4,
            atag=fillstr(atag,4);
        end;
        ltag=' ';
        cid='A';
        writeout=true;
        if ~backbone,
            if strcmpi(strtrim(atag),'N') || strcmpi(strtrim(atag),'C') || strcmpi(strtrim(atag),'O') || strcmpi(strtrim(atag),'CA'),
                writeout=false;
            end;
        end;
        if writeout,
            rotamer=sprintf('%s%s%5i %4s%s%s %s%4s   %8.3f%8.3f%8.3f',rotamer,record,atnum,atag,ltag,label,cid,rid,xyz);
            rotamer=sprintf('%s%6.2f%6.2f          %2s\n',rotamer,pop(rotnum),Bfactor,element);
        end;
    end;
end;
rotamer=sprintf('%sEND   \n',rotamer);

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');