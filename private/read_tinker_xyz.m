function [ecoor,connect,types]=read_tinker_xyz(fname,num)
%
% Reads Tinker xyz file
%
% fname   name of input file, if the file name has no extension, .xyz is
%         appended
% num     optional number of the processing step, if present '_(num)' is
%         appended to the file name
% ecoor   extended coordinates, 1st column is the atom number, columns 2:4
%         are Cartesian coordinates in Angstroem
% connect connection list, numbers of atoms that are connected to a given
%         atom
% types   force field atom types
%
% (c) G. Jeschke, 2008

pse=' HHeLiBe B C N O FNeNaMgAlSi P SClAr KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMdNoLr';

k=strfind(fname,'.');
if isempty(k)
    fname=strcat(fname,'.xyz');
end

if exist('num','var') && ~isempty(num)
    fname=strcat(fname,sprintf('_%i',num));
end

rfile=fopen(fname,'r');

line=fgets(rfile);
[num,rem]=strtok(line);
m=str2double(num);
ecoor=zeros(m,4);
types=zeros(1,m);
connect=zeros(m,8);
for kk=1:m,
    line=fgets(rfile);
    if line==-1,
        disp('### ERROR: Archive file shorter than expected. No output written. ###');
        fclose(rfile); return;
    end;
    if length(line)<53,
        disp('### ERROR: Archive file corrupted, coordinates missing. No output written. ###');
        fclose(rfile); return;
    end;
    [num,rem]=strtok(line);
    [symb,rem]=strtok(rem);
    symb=deblank(symb);
    symb = symb(1);
    if length(symb)<2, symb=[' ' symb]; end;
    element=(findstr(pse,symb)+1)/2;
    ecoor(kk,1)=element;
    [xstr,rem]=strtok(rem);
    [ystr,rem]=strtok(rem);
    [zstr,rem]=strtok(rem);
    aco=[str2double(xstr) str2double(ystr) str2double(zstr)];
    if sum(isnan(aco)),
        disp('### ERROR: Archive file corrupted, wrong coordinate format. No output written. ###');
        fclose(rfile); return;
    end;
    ecoor(kk,2:4)=aco;
    [typstr,rem]=strtok(rem);
    typ=str2double(typstr);
    types(kk)=typ;
    if ~isempty(rem)
     conn=str2num(rem);
    end;
    connect(kk,1:length(conn))=conn;
end;
fclose(rfile);
