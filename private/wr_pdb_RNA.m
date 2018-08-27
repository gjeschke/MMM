function fname = wr_pdb_RNA(fname,ecoor,atomtags,seq,seqnum,code)
% function fname = wr_pdb_RNA(fname,ecoor,atomtags,seq,seqnum,code)

% append proper extension, if extension is missing
if ~contains(fname,'.'), fname=strcat(fname,'.pdb'); end

if ~exist('seqnum','var') || isempty(seqnum)
    seqnum = ecoor(1,1):ecoor(1,1)+length(seq)-1;
end

if ~exist('code','var')
    code = ones(1,length(seq));
end

idCode = 'XXXX';
% generate header line
header=sprintf('HEADER    RNA');
header=fillstr(header,50);
today=date;
today=[today(1:7) today(10:11)];
header=sprintf('%s%s   %s',header,today,idCode);
% state supported format and originating program
format=sprintf('REMARK   4 %s COMPLIES WITH FORMAT V. 3.20, 01-DEC-08',idCode);
origin=sprintf('REMARK   5 WRITTEN BY MMM (MODELLED STRUCTURE)');

fid = fopen(fname,'wt');
if fid == -1
    return;
end

fprintf(fid,'%s\n',header);
fprintf(fid,'TITLE     RNA MODELLED by MMM                                                  \n');
fprintf(fid,'REMARK   4\n%s\n',format);
fprintf(fid,'REMARK   5\n%s\n',origin);

cid = 'A';
[m,~] = size(ecoor);
for k = 1:m
    seqpoi = find(seqnum == ecoor(k,1));
    rtag = seq(seqpoi);
    atag = atomtags{k};
    fprintf(fid,'ATOM  %5i %4s   %s %s%4i    %8.3f%8.3f%8.3f',k,atag,rtag,cid,ecoor(k,1),ecoor(k,2:4));
    fprintf(fid,'%6.2f%6.2f          %2s\n',1.0,code(seqpoi),atag(1));
end
fprintf(fid,'TER   %5i %4s   %s %s%4i   %24s',k+1,'   ',rtag,cid,ecoor(k,1),' ');

fclose(fid);

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');
% newstring = char(padarray(uint8(string)', newlength-length(string), 32,'post')');

function doesit = contains(mystring,pattern)
% Matlab 2017 has this, but it is defined here for compatibility

poi = strfind(mystring,pattern);
doesit = ~isempty(poi);
