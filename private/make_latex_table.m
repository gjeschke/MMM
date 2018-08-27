function make_latex_table(inname,outname)

in_fid=fopen(inname,'r');

fid=fopen(outname,'w');

while 1
    tline = fgetl(in_fid);
    if ~ischar(tline), break, end
     myline = textscan(tline,'%s');
     args=myline{1};
     for k=1:length(args)-1,
         fprintf(fid,'%s &',char(args(k)));
     end;
     fprintf(fid,'     \\hline\n');
     fprintf(fid,'%s \\\\ \n',char(args(end)));
end;
fclose(fid);
fclose(in_fid);