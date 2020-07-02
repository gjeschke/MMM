function varargout = sites_window(varargin)
% SITES_WINDOW M-file for sites_window.fig
%      SITES_WINDOW, by itself, creates a new SITES_WINDOW or raises the existing
%      singleton*.
%
%      H = SITES_WINDOW returns the handle to a new SITES_WINDOW or the handle to
%      the existing singleton*.
%
%      SITES_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SITES_WINDOW.M with the given input arguments.
%
%      SITES_WINDOW('Property','Value',...) creates a new SITES_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sites_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sites_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sites_window

% Last Modified by GUIDE v2.5 12-Jan-2018 12:58:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sites_window_OpeningFcn, ...
                   'gui_OutputFcn',  @sites_window_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sites_window is made visible.
function sites_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sites_window (see VARARGIN)

% Choose default command line output for sites_window
handles.output = hObject;

global model
% global MMM_icon
global hMain

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.output = hObject;

stag=id2tag(model.current_structure,model.structure_tags,model.structure_ids);
handles.preamble=['[' stag ']'];

set(handles.my_figure,'Name',sprintf('Inspect or define sites for structure %s',handles.preamble));

numsites=0;
if isfield(model.info{model.current_structure},'sites'),
    handles.sites=model.info{model.current_structure}.sites;
    numsites=length(handles.sites);
    first_tag=get_tag(handles.sites);
else
    first_tag='AC1';
end;

handles.numsites=numsites;
handles.edited=0;

handles.sites(numsites+1).tag=first_tag;
handles.sites(numsites+1).residues={};
handles.sites(numsites+1).evidence='AUTHOR';
handles.sites(numsites+1).declaration='<New site>';
handles.sites(numsites+1).evidence_type=2;

[selection,msg]=resolve_address('*');
[mm,nn]=size(selection);
poi=0;
for k=1:mm,
    cindices=selection(k,:);
    cindices=cindices(cindices>0);
    if length(cindices)==4 && cindices(1)==model.current_structure, % residue selection
        adr=mk_address(cindices);
        if ~isempty(adr),
            poi=poi+1;
            pp=strfind(adr,']');
            if pp<length(adr),
                handles.sites(numsites+1).residues{poi}=adr(pp+1:end);
            end;
        end;
    end;
end;

if numsites==0 && isempty(handles.sites(1).residues),
    add_msg_board('WARNING: No sites defined in this structure and no residues selected.');
    add_msg_board('Please select resdiues to define a new site.');    
    delete(hObject);
    return;
elseif isempty(handles.sites(numsites+1).residues),
    handles.sites=handles.sites(1:numsites);
    set(handles.pushbutton_new,'Enable','off');
end;

handles=update_sitelist(handles);

poi=get(handles.listbox_sites,'UserData');
site_id=poi(1);

handles.current=site_id;

adr_list{1}='';
if ~isempty(handles.sites(site_id).residues),
    m=length(handles.sites(site_id).residues);
    for k=1:m,
        address=sprintf('%s%s',handles.preamble,handles.sites(site_id).residues{k});
        adr_list{k}=address;
    end;
end;
adr_list=sort(adr_list);

set(handles.listbox_residues,'Value',1);
set(handles.listbox_residues,'String',adr_list);
set(handles.listbox_residues,'Min',1);
set(handles.listbox_residues,'Max',m);

etype=handles.sites(site_id).evidence_type;
set(handles.popupmenu_evidence,'Value',etype);
set(handles.edit_description,'String',handles.sites(site_id).declaration);


load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sites_window wait for user response (see UIRESUME)
% uiwait(handles.my_figure);


% --- Outputs from this function are returned to the command line.
function varargout = sites_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on selection change in listbox_sites.
function listbox_sites_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_sites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_sites contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sites

global model

site_id=get(hObject,'Value');
poi=get(hObject,'UserData');
site_id=poi(site_id);

handles.current=site_id;

adr_list{1}='';
if ~isempty(handles.sites(site_id).residues),
    m=length(handles.sites(site_id).residues);
    for k=1:m,
        address=sprintf('%s%s',handles.preamble,handles.sites(site_id).residues{k});
        cindices=resolve_address(address);
        resname=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).name;
        address=[address '; ' resname];
        adr_list{k}=address;
    end;
end;
adr_list=sort(adr_list);

set(handles.listbox_residues,'Value',1);
set(handles.listbox_residues,'String',adr_list);
set(handles.listbox_residues,'Min',1);
set(handles.listbox_residues,'Max',m);

etype=handles.sites(site_id).evidence_type;
set(handles.popupmenu_evidence,'Value',etype);
set(handles.edit_description,'String',handles.sites(site_id).declaration);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_sites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_sites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_residues.
function listbox_residues_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_residues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_residues contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_residues


% --- Executes during object creation, after setting all properties.
function listbox_residues_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_residues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if handles.edited,
    ButtonName = questdlg('Close window and cancel edits of new site?', 'New site edited but not stored.', 'Yes', 'No', 'No');
    if strcmp(ButtonName,'Yes'),
        model.info{model.current_structure}.sites=handles.sites(1:handles.numsites);
        delete(handles.my_figure);
    end;
else
    model.info{model.current_structure}.sites=handles.sites(1:handles.numsites);
    delete(handles.my_figure);
end;



% --- Executes on button press in pushbutton_annotation.
function pushbutton_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

sel=get(handles.listbox_residues,'Value');

contents = get(handles.listbox_residues,'String');
address = contents{sel(1)};

indices=resolve_address('*');
button='Yes';
if ~isempty(indices),
    button = questdlg('Existing selections will be cancelled.','Are you sure?','No','Yes','No');
end;
set(handles.my_figure,'Pointer','watch');
drawnow;
if strcmp(button,'Yes'),
    hMain=cmd(hMain,'unselect *');
    hMain=cmd(hMain,sprintf('select %s',address));
    hMain.keyword_request='binding sites';
    annotation_window;
end;
set(handles.my_figure,'Pointer','arrow');
guidata(hObject,handles);

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global hMain

sel=get(handles.listbox_residues,'Value');
contents = get(handles.listbox_residues,'String');
for k=1:length(sel),
    address = contents{sel(k)};
    hMain=cmd(hMain,sprintf('select %s',address));
end;

% --- Executes on button press in pushbutton_select_all.
function pushbutton_select_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

contents = get(handles.listbox_residues,'String');
for k=1:length(contents),
    address = contents{k};
    hMain=cmd(hMain,sprintf('select %s silent',address));
end;
if hMain.hierarchy_display,
%     if hMain.large,
%         hwhandle=hMain.hierarchy_window_large;
%         hhandles=guidata(hwhandle);
%         hierarchy_window_large('sequence_display',hhandles);
%     else
%         hwhandle=hMain.hierarchy_window;
%         hhandles=guidata(hwhandle);
%         hierarchy_window('sequence_display',hhandles);
%     end;
    hwhandle=hMain.hierarchy_window;
    hhandles=guidata(hwhandle);
    hierarchy_window('sequence_display',hhandles);
end;

% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain=cmd(hMain,'unselect *');

% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'sites_window.html');
webcall(entry,'-helpbrowser');


% --- Executes on button press in pushbutton_new.
function pushbutton_new_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.numsites=handles.numsites+1;
set(hObject,'Enable','off');
handles.edited=0;
if ~isempty(handles.sites(handles.numsites).residues),
    text=[handles.sites(handles.numsites).tag ': ' handles.sites(handles.numsites).declaration];
    evidence_type='conjecture';
    switch handles.sites(handles.numsites).evidence_type
        case 1
            evidence_type='atomistic structure';
        case 2
            evidence_type='conjecture';
        case 3
            evidence_type='homology';
        case 4
            evidence_type='mutation/assay';
        case 5
            evidence_type='other biochemistry';
        case 6
            evidence_type='spectroscopy';
    end;
    text=strvcat(text,sprintf('found by %s',evidence_type));
    for k=1:length(handles.sites(handles.numsites).residues),
        indices=resolve_address(handles.sites(handles.numsites).residues{k});
        if ~isempty(indices),
            add_annotation(indices,'Binding',text,{'binding sites'});
        end;
    end;
end;
guidata(hObject,handles);


function edit_description_Callback(hObject, eventdata, handles)
% hObject    handle to edit_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_description as text
%        str2double(get(hObject,'String')) returns contents of edit_description as a double

handles.sites(handles.current).declaration=get(hObject,'String');
poi=get(handles.listbox_sites,'UserData');
for k=1:length(poi),
    if poi(k)==handles.current, break; end;
end;
if handles.current>handles.numsites,
    handles.edited=1;
end;
sd=get(handles.listbox_sites,'String');
sd{k}=sprintf('%s: %s',handles.sites(handles.current).tag,get(hObject,'String'));
set(handles.listbox_sites,'String',sd);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_description_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_evidence.
function popupmenu_evidence_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_evidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_evidence contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_evidence

handles.sites(handles.current).evidence_type=get(hObject,'Value');
if handles.current>handles.numsites, handles.edited=1; end;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_evidence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_evidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=update_sitelist(handles)

sitelist={};
for k=1:length(handles.sites),
    sitelist{k}=sprintf('%s: %s',handles.sites(k).tag,handles.sites(k).declaration);
end;
[sitelist,poi]=sort(sitelist);
set(handles.listbox_sites,'String',sitelist);
set(handles.listbox_sites,'UserData',poi);
set(handles.listbox_sites,'Value',1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ButtonName = questdlg('Close window and cancel all changes?', 'All sites changes will be lost.', 'Yes', 'No', 'Yes');
if strcmp(ButtonName,'Yes'),
    delete(handles.my_figure);
end;

function newtag=get_tag(sites)

n=length(sites);

c1=1;
c2=3;
c3=1;
for k=1:n,
    if double(sites(k).tag(1))-64>c1,
        c1=double(sites(k).tag(1))-64;
    end;
end;
for k=1:n,
    if double(sites(k).tag(1))-64==c1,
        if double(sites(k).tag(2))-64>c2,
            c2=double(sites(k).tag(2))-64;
            if c2>26, c2=3; c1=c1+1; end;
        end;
    end;
end;
for k=1:n,
    if double(sites(k).tag(1))-64==c1,
        if double(sites(k).tag(2))-64==c2,
            if double(sites(k).tag(3))-48>c3,
                c3=double(sites(k).tag(2))-48;
                if c3>9, c3=1; c2=c2+1; end;
                if c2>26, c2=3; c1=c1+1; end;
            end;
        end;
    end;
end;
if c1>26,
    c1=1;
    c2=1;
    c3=1;
end;

found=1;
while found && c1<27,
    found=0;
    newtag=[char(c1+64) char(c2+64) char(c3+48)];
    for k=1:n,
        if strcmp(newtag,sites(k).tag),
            found=1;
            c3=c3+1;
            if c3>9, 
                c3=1;
                c2=c2+1;
                if c2>26,
                    c2=3;
                    c1=c1+1;
                end;
            end;
        end;
    end;
end;
