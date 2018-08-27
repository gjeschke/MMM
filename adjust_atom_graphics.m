function adjust_atom_graphics(automatic)
% adjusts atom graphics display according to number of displayed objects
%
% automatic     flag, true: use automatic values, false: revert to values
%               requested by user
% G. Jeschke, 2009

global hMain
global graph_settings

hMain.text_atom_request.Position = [0.75,0.05,0.18,0.4];

if automatic==hMain.atom_graphics_auto, % no change
    return
else
    hMain.atom_graphics_auto=automatic;
end;

if automatic,
    hMain.atom_graphics_requested=hMain.atom_graphics_mode;
    switch hMain.atom_graphics_requested,
        case 0
            set(hMain.text_atom_request,'String','off requested','ForegroundColor','r');
        case 1
            set(hMain.text_atom_request,'String','on requested','ForegroundColor',[0,0.6,0]);
        case 2
            set(hMain.text_atom_request,'String','low requested','ForegroundColor',[0.75,0.6,0]);
    end;
    full=0;
    try
        full=numel(get(hMain.atom_graphics,'Children'));
    catch except
    end;
    all_axis=numel(get(hMain.axes_model,'Children'));
    if full<=graph_settings.max_full_objects && all_axis<=graph_settings.max_full_objects,
        mode=hMain.atom_graphics_mode;
    else
        reduced=0;
        try
            reduced=numel(get(hMain.atom_graphics_reduced,'Children'));
        catch except
        end;
        if reduced<=graph_settings.max_reduced_objects,
            mode=2;
            if hMain.atom_graphics_mode==0,
                mode=0;
            end;
        else
            mode=0;
        end;
    end;
    set(hMain.togglebutton_atom_graphics,'Enable','off');
else
    set(hMain.text_atom_request,'String','');
    hMain.atom_graphics_mode=hMain.atom_graphics_requested;
    mode=hMain.atom_graphics_mode;
    set(hMain.togglebutton_atom_graphics,'Enable','on');
end;
ison=mode>0;
set(hMain.togglebutton_atom_graphics,'Value',ison);

switch_it=true;
if ~isfield(hMain,'atom_graphics') ||  ~ishandle(hMain.atom_graphics),
    switch_it=false;
end;
if ~isfield(hMain,'atom_graphics_reduced') ||  ~ishandle(hMain.atom_graphics_reduced),
    switch_it=false;
end;

add_msg_board('Display update may take a while...');
switch mode
    case 0
        if switch_it,
            set(hMain.atom_graphics,'Visible','off');
            set(hMain.atom_graphics_reduced,'Visible','off');
        end;
        hMain.atom_graphics_mode=0;
        set(hMain.togglebutton_atom_graphics,'String','Atom graphics off','ForegroundColor','r');
    case 1
        if switch_it,
            set(hMain.atom_graphics_reduced,'Visible','off');
            set(hMain.atom_graphics,'Visible','on');
        end;
        hMain.atom_graphics_mode=1;
        set(hMain.togglebutton_atom_graphics,'String','Atom graphics on','ForegroundColor',[0,0.6,0]);
    case 2
        if switch_it,
            set(hMain.atom_graphics,'Visible','off');
            set(hMain.atom_graphics_reduced,'Visible','on');
        end;
        hMain.atom_graphics_mode=2;
        set(hMain.togglebutton_atom_graphics,'String','Atom graphics low','ForegroundColor',[0.75,0.6,0]);
end;
