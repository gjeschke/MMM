function cartoon_object_clicked(cb,eventdata)
% cartoon_object_clicked(cb,eventdata)
%
% Function executed when a user clicks on a cartoon graphics object

% global model
global hMain

mouseclick=get(gcf,'SelectionType');

if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')
    color_selection;
    if isempty(hMain.color_selection),
        add_msg_board('Color selection cancelled by user.');
    else
        if isfloat(hMain.color_selection),
            set(cb,'FaceColor',hMain.color_selection);
        else
            add_msg_board('Color schemes are not implemented for cartoons.');
        end;
    end;
end;

if strcmpi(mouseclick,'extend'),
    transparency_selection;
    if ~isempty(hMain.alpha),
        set(cb,'FaceAlpha',hMain.alpha);
    else
        add_msg_board('Transparency selection cancelled by user.');
    end;
end;

% mouseclick=get(gcf,'SelectionType');

% if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')
% 
%     % Deselect old selection
%     if isfield(model,'selected') && ~isempty(model.selected), % was anything selected
%         hMain=cmd(hMain,'unselect *');
%         model.secondary_selected={};
%     end;
% end;
% 
% if strcmpi(mouseclick,'open')
%     add_msg_board('--- All deselected ---');
%     return
% end;



