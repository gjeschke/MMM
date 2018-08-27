function text=rd_txt(fname)

fid=fopen(fname);
if fid==-1,
    info.no_file=1;
    return;
end;

text='';
nl=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    text=[text tline '\n'];
end;