function add_msg_board(msg)
% adds a one-line message to the message board

global hMain

max_lines=25;
max_columns=60;

if ~isa(msg,'char'),
    return;
end;

if ~isfield(hMain,'text_message_board') || ~ishandle(hMain.text_message_board),
    % fprintf(1,'%s\n',msg);
    return
end;

handles=guidata(hMain.text_message_board);
msg_board=get(hMain.text_message_board,'String');
if isempty(msg_board),
    msg_board{1}='--- Welcome to MMM ---';
elseif ischar(msg_board),
    msg0=msg_board;
    clear msg_board;
    msg_board{1}=msg0;
end;

words=textscan(msg,'%s');
words=words{1};
poi=0;
if length(words)>0,
    % cut words that are too long
    for k=1:length(words),
        cword=char(words(k));
        while length(cword)>max_columns,
            poi=poi+1;
            nwords{poi}=cword(1:max_columns);
            cword=cword(max_columns+1:end);
        end;
        poi=poi+1;
        nwords{poi}=cword;
    end;
    % now display, considering max width of text line
    line='';
    n=length(msg_board);
    poi=n;
    for k=1:length(nwords),
        if length(line)+length(nwords{k})<=max_columns-1,
            line=[line nwords{k} ' '];
        else
            poi=poi+1;
            msg_board{poi}=line;
            line=[nwords{k} ' '];
        end;
    end;
    if ~isempty(line),
        poi=poi+1;
        msg_board{poi}=line;
    end;    
end;
% now restrict to max number of lines
if poi>max_lines,
    n=poi;
    msg_board0=msg_board;
    clear msg_board
    poi=0;
    for k=n-max_lines+1:n,
        poi=poi+1;
        msg_board{poi}=msg_board0{k};
    end;
end;
set(hMain.text_message_board,'String',msg_board);
guidata(hMain.text_message_board,handles);

if get(hMain.checkbox_log,'Value'),
    fid=fopen(hMain.logfile,'a+');
    fprintf(fid,'msg> %s\n',msg);
    fclose(fid);
end;
drawnow;