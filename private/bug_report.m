function bug_report

global hMain

fid=fopen(hMain.logfile,'a+');
if ~get(hMain.checkbox_log,'Value'),
    fprintf(fid,'BUG> logging was switched off when bug happened.\n');
end;
fprintf(fid,'BUG> follows the output into the Matlab command window:\n');
fclose(fid);
diary(hMain.logfile);
diary off
fid=fopen(hMain.logfile,'a+');
fprintf(fid,'BUG>--- end of Matlab command window output ---\n');
fclose(fid);
hMain.report_file=hMain.logfile;
report_editor;
web('mailto:gjeschke@ethz.ch');
diary on
