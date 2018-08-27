function ascii_res=distr2ASCII(ascii_name,X,Y,format)

% make *.dat ASCII file from X-Y dependent data
% X - (1,m) vector, used as X data
% Y - (1,m) vector, used as Y data
% ascii_ name - string, used as a name for the output file
% format - optional format string, must accept two numbers as output
% arguments, defaults to '%7.4f %7.4f'
%
% Ye. Polyhach, 2009
%--------------------------------------------------------------------------

if nargin<4, % added by G. Jeschke, 29.10.2020
    format='%7.4f %7.4f'; 
end;

length_check=length(X)==length(Y);
if length_check==0
    disp(sprintf('%s','lengths of vectors X and Y are not the same!'));
    return
end

record=zeros(length(X),2); % initialize array to be recorded as ASCII;
ascii_name=strcat(ascii_name,'.dat');
fid=fopen(ascii_name, 'wt');
record(:,1)=X';
record(:,2)=Y';
record=record';
fprintf(fid,format,record);
fclose(fid);

ascii_res=sprintf('%s%s','ASCII data ',ascii_name,' have been created.');
disp(ascii_res);