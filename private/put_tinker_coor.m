function outname=put_tinker_coor(fname,coor,connect,types)

% make Tinker *.xyz input file from a coordinate set and known connection
% table
%
% fname - output file name basis (without extention)
% coor - xyz coordinates (row 2-4) with atomic numbers (row 1)
% connect - connections (get from the template file prepared beforehand)
% types - preassigned atom type indexes for the chosen forcefield

% keyboard

% Ye. Polyhach, 2009

pse=' HHeLiBe B C N O FNeNaMgAlSi P SClAr KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMdNoLr';

[m,n]=size(coor);
[mc,nc]=size(connect);
outname=strcat(fname,'.xyz');
wfile=fopen(outname,'w+');
line=sprintf(' %i',m);
fprintf(wfile,'%s\n',line);
for k=1:m,
    atn=k;
    label=sprintf('%6u',atn);
    element=2*coor(k,1)-1;
    symb=pse(element:element+1);
    symb=sprintf('%3s',symb');
    line=strcat(label,symb);
    xc=coor(k,2);
    yc=coor(k,3);
    zc=coor(k,4);
    xcoor=sprintf('%11.5f',xc);
    ycoor=sprintf('%11.5f',yc);
    zcoor=sprintf('%11.5f',zc);
    line=strcat(line,xcoor,ycoor,zcoor);
    type0=types(k);
    type0line=sprintf('%5.0f',type0);
    line=strcat(line,type0line);
    connect0=connect(k,:);
    connline=sprintf('%i');
    nonzero=find(connect0);
    for kc=1:length(nonzero)
        conn_entry=connect0(nonzero(kc));
        connline0=sprintf('%5.0f',conn_entry);
        connline=strcat(connline,connline0);
    end
    
%     for kc=1:nc
%         conn_entry=connect0(nc);
%         if conn_entry==0
%             return
%         end
%         connline0=sprintf('%5.0f',conn_entry);
%         connline=strcat(conn_start,connline0);
%     end
    line=strcat(line,connline);
    fprintf(wfile,'%s\n',line);
end;
fclose(wfile);
% 
% [m,n]=size(connect);
% for k=1:m,
%    line=sprintf('%s%5u','CONECT',k);
%    numcon=length(find(connect(k,:)));
%    for l=1:numcon,
%       bond=sprintf('%5u',connect(k,l));
%       line=[line bond];
%    end;
%    fprintf(wfile,'%s\n',line);
% end;
% fprintf(wfile,'%s\n','END');
% 
% fclose(wfile);
