function [errors,warnings,variables,values]=analyze_modeller_log(logfile,echo)
% function [errors,warnings,variables,values]=analyze_modeller_log(logfile,echo)
%
% Extracts warning and error messages as well as output variables and their
% values from a Modeller log file
%
% lines with output variable must start with &, followed by the variable
% name, the variable name must be limited by one of the characters space,
% comma, semicolon, = or full stop, the delimiter must be followed by a
% value or expression that can be interpreted by Matlab function str2num
% (and nothing else), if the expression cannot be interpreted, an empty
% array is returned for the value
%
% logfile   filename of the Modeller log file
% echo      (optional) flag, true: warnings and errors are output to
%           message board, false: routine is silent, defaults to false
%
% errors    cell array of error message strings
% warnings  cell array of warning strings
% variables cell array of variable names
% values    cell array of corresponding variable values
%
% G. Jeschke, 2011

errors={};
warnings={};
variables={};
values={};

if nargin<2,
    echo=false;
end;

fid=fopen(logfile);
if fid==-1,
    errors{1}='GENERAL: No log file found';
    if echo,
        add_msg_board('### Modeller ERRORS ###');
        add_msg_board(errors{1});
    end;
    return;
end;

nl=0;
ne=0;
nw=0;
nv=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    if length(tline)>1,
        if strcmp(tline(1),'&'),
            [myvar,rem]=strtok(tline(2:end),' ,;=.');
            if ~isempty(rem),
                while strfind(' ,;=.',rem(1)),
                    if length(rem)<2,
                        rem=''; break;
                    else
                        rem=rem(2:end);
                    end;
                end;
            end;
            if ~isempty(rem),
                nv=nv+1;
                variables{nv}=myvar;
                values{nv}=str2num(rem);
            end;
        else
            wpoi=strfind(tline,'_W');
            if ~isempty(wpoi),
                routine=strtok(tline(1:wpoi),'_');
                nw=nw+1;
                warnstr=[routine ':' tline(wpoi+2:end)];
                wline=' ';
                while ~isempty(wline),
                    wline = fgetl(fid);
                    if ~ischar(wline), break, end
                    nl=nl+1;
                    warnstr=[warnstr '|' wline];
                end;
                warnings{nw}=warnstr;
                tline='';
            end;
            epoi=strfind(tline,'_E');
            if ~isempty(epoi),
                routine=strtok(tline(1:epoi),'_');
                ne=ne+1;
                errstr=[routine ':' tline(epoi+2:end)];
                eline=' ';
                while ~isempty(eline),
                    eline = fgetl(fid);
                    if ~ischar(eline), break, end
                    nl=nl+1;
                    errstr=[errstr '|' eline];
                end;
                errors{ne}=errstr;
            end;
        end;
    end;
end;

if echo,
    if ~isempty(warnings),
        add_msg_board('--- Modeller warnings ---');
        for k=1:length(warnings),
            list2message_board(warnings{k});
        end;
    end;
    if ~isempty(errors),
        virgin=true;
        for k=1:length(errors),
            nonsense=textscan(errors{k},'%s','Delimiter','|');
            items=nonsense{1};
            if ~strcmpi(char(items{1}),'openf:PR.lib'), % not a Modeller ERROR, but output reported as ERROR format
                if virgin,
                    add_msg_board('### Modeller ERRORS ###');
                    virgin=false;
                end;
                list2message_board(errors{k});
            end;
        end;
    end;        
end;

function list2message_board(strlist)

nonsense=textscan(strlist,'%s','Delimiter','|');
items=nonsense{1};
for k=1:length(items),
    add_msg_board(items{k});
end;