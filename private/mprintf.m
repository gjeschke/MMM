function mprintf(fid,line)
% mprintf(fid,line)
%
% Writes output line to file or, if fid == 1, to the message board 
% behaves like fprintf with command window output being redirected to 
% message board and no further arguments, with line feed in file
%
% G. Jeschke, 2014

if fid == 1,
    add_msg_board(sprintf(line));
else
    fprintf(fid,[line '\n']);
end;