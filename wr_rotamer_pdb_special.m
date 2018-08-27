function wr_rotamer_pdb_special
% function wr_rotamer_pdb(adr,rotamer_stats,T,label,MD2PDB,library)
%
% writes a PDB file with all rotamers of a spin labeled residue 
% the file name is derived from the residue address adr, the label name 
% and the labeling temperature T
% 
% further input parameters:
% rotamer_stats     the rotamer statistics as computed by get_rotamers.m
% MD2PDB            translation table for atom numbers
% library           rotamer library name

global label_defs
global chemistry

record='HETATM';

load R1A_298K_090619

label='R1A';

corr=load('rotlib_PDB_correspondence.dat');

rmsd=1; % standard atom r.m.s.d. in Å for B factor specification
Bfactor=8*pi^2*rmsd^2;

labnum=tag2id(label,label_defs.restags);
atags=label_defs.residues(labnum).atoms;
mylabel.tc=label_defs.residues(labnum).tc;
conn=label_defs.residues(labnum).conn;
[mconn,nconn]=size(conn);

fid=fopen('rotamer1_R1A_298K.pdb','w');

fprintf(fid,fillstr(sprintf('HEADER    SPIN LABELING WITH LIBRARY %s','R1A_298K_090619'),80)); 
fprintf(fid,'\n');
fprintf(fid,fillstr(sprintf('TITLE     LABEL R1A (MTSSL) AT 298 K'),80)); 
fprintf(fid,'\n');
fprintf(fid,fillstr('REMARK 4',80));
fprintf(fid,'\n');
fprintf(fid,fillstr('REMARK 4 MMM ROTAMERS (PSEUDO-PDB FORMAT)',80));
fprintf(fid,'\n');
atnum=0;

ecoor=rot_lib.library(1).ecoor;
ecoor0=ecoor;
ecoor0=[ecoor0;0,0,0,0];
ecoor=[ecoor;1,0,0,0];

for k=1:length(corr);
    ecoor0(k,1)=corr(k,2);
end;
[sorted,poi]=sort(ecoor0(:,1));
ecoor=ecoor(poi,:);
% ecoor=ecoor(3:end,:);

[atoms,n]=size(ecoor);

rid=sprintf('%4i ',1);
for kk=1:atoms,
    atag=id2tag(kk,atags);
    if ~strcmpi(atag,'OXT') && ~strcmpi(atag,'HXT') && ~sum(isnan(ecoor(kk,:))),
        atnum=atnum+1;
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
        fprintf(fid,'%s%5i %4s%s%s %s%4s   %8.3f%8.3f%8.3f',record,atnum,atag,ltag,label,cid,rid,xyz);
        fprintf(fid,'%6.2f%6.2f          %2s\n',1.0,Bfactor,element);
    end;
end;
fprintf(fid,'END                                                          \n');
fclose(fid);

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');