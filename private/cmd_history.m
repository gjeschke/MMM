function [handles,veto]=cmd_history(handles,command,undo_cmd)
% handles=cmd_history(handles,command,undo_cmd)
%
% Stores a command in the history list and the corresponding inverse
% command in the undo list,
% if no undo command or an empty undo command is supplied, the user is
% asked whether he wants to proceed with this command, his answer is given
% back in the flag 'veto'
%
% handles  GUI handles structure
% command  comand to execute
% undo_cmd inverse command for undo function, use 'noundo' to suppress
%          question to user
% veto     flag that informs on user veto, 0 command accepted despite of
%          missing undo, 1 command rejected
%
% (c) G. Jeschke, 2009
global hMain

veto=0;
if (nargin<3 || isempty(undo_cmd)) && hMain.store_undo % undo impossible, but requested
    answer=questdlg('Do you want to execute this command?','No undo available','Execute','Cancel','Execute');
    switch answer
        case 'Execute'
            hMain.history{hMain.history_poi}=command;
            hMain.edit_poi=hMain.history_poi;
            hMain.history_poi=hMain.history_poi+1;
            if hMain.history_poi>hMain.history_num,
                hMain.history_poi=1;
            end;
            % The undo list must be flushed to avoid undefined states of
            % the program
            hMain.undo_poi=1;
            for k=1:hMain.history_num, hMain.undo{k}=''; end;
        case 'Cancel'
            veto=1;
    end;
else
    if hMain.store_undo, % to allow for history silent commands, for instance in calls by 'undo' itself
        hMain.history{hMain.history_poi}=command;
        hMain.edit_poi=hMain.history_poi;
        hMain.history_poi=hMain.history_poi+1;
        if hMain.history_poi>hMain.history_num,
            hMain.history_poi=1;
        end;
        if ~strcmpi(undo_cmd,'noundo')
            hMain.undo{hMain.undo_poi}=undo_cmd;
            hMain.undo_poi=hMain.undo_poi+1;
        end;
        if hMain.undo_poi>hMain.history_num,
            hMain.undo_poi=1;
        end;
    end;
    if nargin>2 && strcmpi(undo_cmd,'noundo'),
        % flush undo list
         hMain.undo_poi=1;
         for k=1:hMain.history_num, hMain.undo{k}=''; end;
    end;
end;  

hMain.store_undo=1;


        