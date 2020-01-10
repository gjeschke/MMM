function [tags,colors] = mk_color_list(fname)

ifid=fopen(fname);
if ifid==-1
    warning('rd_CA_trace:Input file could not be opened');
    return;
end

tags = ':';
colors = ':';
while 1
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    myline = textscan(tline,'%s');
    args=myline{1};
    tags = [tags char(args(1)) ':'];
    colors = [colors char(args(2)) ':'];
end
fclose(ifid);