function wr_rotamer_pdb(adr,rotamer_stats,T,label,MD2PDB,library)
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

% load Rp_std_ORCA

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

mylabel.tc=label_defs.residues(labnum).tc;
conn=label_defs.residues(labnum).conn;
[mconn,nconn]=size(conn);

pop=rotamer_stats.all_potentials.pop_rot;
nrot=length(pop);

% fid_f=fopen(['frames_' adr '_' label '_' sprintf('%i',T) '.dat'],'w');
if nrot <= 1024,
    fid=fopen(['rotamers_' adr '_' label '_' sprintf('%i',T) '.pdb'],'w'); % for large libraries this can lead to a crash of Matlab 2016b
else
    fid=fopen(['rotamers_' adr '_' label '_' sprintf('%i',T) '.pdb'],'W');
end;

fprintf(fid,fillstr(sprintf('HEADER    SPIN LABELING WITH LIBRARY %s',library),80)); 
fprintf(fid,'\n');
fprintf(fid,fillstr(sprintf('TITLE     %s LABELED WITH %s AT %i K',adr,label,T),80)); 
fprintf(fid,'\n');
fprintf(fid,fillstr('REMARK 4',80));
fprintf(fid,'\n');
fprintf(fid,fillstr('REMARK 4 MMM ROTAMERS (PSEUDO-PDB FORMAT) RESIDUE NUMBERS ARE ROTAMER IDS',80));
fprintf(fid,'\n');
atnum=0;
basnum=zeros(1,nrot);

maxatags = 0;
found = true;
while found
    maxatags = maxatags + 1;
    atag=id2tag(maxatags,atags);
    if isempty(atag),
        found = false;
        maxatags = maxatags - 1;
    end;
end;

for k=1:nrot,
    ecoor0=rotamer_stats.all_rotamers{k};
    ecoor=zeros(size(ecoor0));
    [atoms,n]=size(ecoor);
    for kk=1:length(MD2PDB),
        if MD2PDB(kk)>0 && MD2PDB(kk)<=atoms,
            ecoor(kk,:)=ecoor0(MD2PDB(kk),:);
        else
            ecoor(kk,:)=[NaN,NaN,NaN,NaN];
        end;
    end;
    mypop=pop(k);
    
%     if pop(k) >= 0.001,
%         fprintf(fid_f,'%i  %6.4f\n',k,mypop);
%         coor_r = ecoor(:,2:4);
%         [m,~] = size(coor_r);
%         coor_r = coor_r - repmat(coor_r(24,:),m,1);
% 
%         xp = coor_r(23,:)-coor_r(49,:);
%         x = xp/norm(xp);
%         yp = coor_r(68,:)-coor_r(49,:);
%         yp = yp/norm(yp);
%         z = cross(x,yp);
%         z = z/norm(z);
%         y = cross(z,xp);
%         y = y/norm(y);
% 
%         dircos=[x;y;z];
%         Rpl=dircos; % rotation matrix for conversion to standard frame
% 
%         Rpx = Rpb*Rpl;
% 
%         for kk = 1:3,
%             fprintf(fid_f,'%8.4f%8.4f%8.4f\n',Rpx(kk,:));
%         end;
%         
%         coorx2 = coor_r;
%         for kk = 1:m,
%             newvec = Rpx*coor_r(kk,:)';
%             coorx2(kk,:) = newvec';
%         end;
%     end;
        
    rid=sprintf('%4i ',k);
    basnum(k)=atnum;
    translate=zeros(1,mconn);
    pconn=0;
    for kk=1:maxatags,
        atag=id2tag(kk,atags);
        if ~strcmpi(atag,'H2') && ~strcmpi(atag,'OXT') && ~strcmpi(atag,'HXT') && ~sum(isnan(ecoor(kk,:))),
            atnum=atnum+1;
            pconn=pconn+1;
            translate(kk)=pconn;
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
            fprintf(fid,'%6.2f%6.2f          %2s\n',pop(k),Bfactor,element);
        end;
    end;
end;
if strcmpi(record,'HETATM'), % this requires CONECT records
    [m,n]=size(conn);
    for k=1:nrot,
        for kk=1:m,
            if translate(kk)>0,
                bonds=conn(kk,:);
                for kkk=1:length(bonds), if bonds(kkk)~=0, bonds(kkk)=translate(bonds(kkk)); end; end;
                abonds=[];
                for kkk=1:length(bonds), if bonds(kkk)>0, abonds=[abonds bonds(kkk)]; end; end;
                poi=0;
                while poi<length(abonds),
                    tline=sprintf('CONECT%5i',basnum(k)+translate(kk));
                    apoi=0;
                    flush=0;
                    while apoi<4 && poi<length(abonds),
                        poi=poi+1;
                        apoi=apoi+1;
                        flush=1;
                        tline=sprintf('%s%5i',tline,abonds(poi)+basnum(k));
                    end;
                    if flush,
                        fprintf(fid,'%s\n',tline);
                    end;    
                end;
            end;
        end;
    end;
end;
fprintf(fid,'END                                                          \n');
fclose(fid);
% fclose(fid_f);

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');