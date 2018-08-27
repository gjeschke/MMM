function varargout = find_by_key(varargin)
% FIND_BY_KEY M-file for find_by_key.fig
%      FIND_BY_KEY, by itself, creates a new FIND_BY_KEY or raises the existing
%      singleton*.
%
%      H = FIND_BY_KEY returns the handle to a new FIND_BY_KEY or the handle to
%      the existing singleton*.
%
%      FIND_BY_KEY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIND_BY_KEY.M with the given input arguments.
%
%      FIND_BY_KEY('Property','Value',...) creates a new FIND_BY_KEY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before find_by_key_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to find_by_key_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help find_by_key

% Last Modified by GUIDE v2.5 12-Jan-2018 12:45:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @find_by_key_OpeningFcn, ...
                   'gui_OutputFcn',  @find_by_key_OutputFcn, ...
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


% --- Executes just before find_by_key is made visible.
function find_by_key_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to find_by_key (see VARARGIN)

% Choose default command line output for find_by_key

global model
global MMM_icon
global hMain

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.output = hObject;

handles=update_keylist(handles);

poi=get(handles.listbox_keywords,'UserData');
indices=model.keys(poi(1)).indices;

adr_list{1}='';
if ~isempty(indices),
    [m,n]=size(indices);
    for k=1:m,
        cindices=indices(k,:);
        if ~isempty(cindices),
            address=mk_address(cindices,1);
            adr_list{k}=address;
        end;
    end;
end;

adr_list=sort(adr_list);
set(handles.listbox_objects,'String',adr_list);

load helpicon
set(handles.pushbutton_help_find_by_key,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes find_by_key wait for user response (see UIRESUME)
% uiwait(handles.my_figure);


% --- Outputs from this function are returned to the command line.
function varargout = find_by_key_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_keywords.
function listbox_keywords_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_keywords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_keywords contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_keywords

global model

key_id=get(hObject,'Value');
poi=get(hObject,'UserData');
key_id=poi(key_id);

indices=model.keys(key_id).indices;

adr_list{1}='';
if ~isempty(indices),
    [m,n]=size(indices);
    for k=1:m,
        cindices=indices(k,:);
        address=mk_address(cindices,1);
        adr_list{k}=address;
    end;
end;

adr_list=sort(adr_list);
set(handles.listbox_objects,'Value',1);
set(handles.listbox_objects,'String',adr_list);
set(handles.listbox_objects,'Min',1);
set(handles.listbox_objects,'Max',m);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listbox_keywords_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_keywords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_objects.
function listbox_objects_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_objects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_objects contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_objects


% --- Executes during object creation, after setting all properties.
function listbox_objects_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_objects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_okay.
function pushbutton_okay_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_okay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.my_figure);

% --- Executes on button press in pushbutton_annotation.
function pushbutton_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

sel=get(handles.listbox_objects,'Value');

contents = get(handles.listbox_objects,'String');
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
    key_id=get(handles.listbox_keywords,'Value');
    keywords=get(handles.listbox_keywords,'String');
    hMain.keyword_request=keywords{key_id};
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

sel=get(handles.listbox_objects,'Value');
contents = get(handles.listbox_objects,'String');
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

contents = get(handles.listbox_objects,'String');
for k=1:length(contents),
    address = contents{k};
    hMain=cmd(hMain,sprintf('select %s',address));
end;


% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain=cmd(hMain,'unselect *');

function handles=update_keylist(handles)

global model

if isfield(model,'keywords'),
    for k=1:length(model.keys),
        keylist{k}=id2tag(k,model.keywords);
        ukeylist{k}=upper(id2tag(k,model.keywords));
    end;
else
    keylist{1}='';
end;
[discard,poi]=sort(ukeylist);
newkeys={};
for k=1:length(poi),
    newkeys{k}=keylist{poi(k)};
end;
set(handles.listbox_keywords,'String',newkeys);
set(handles.listbox_keywords,'UserData',poi);
set(handles.listbox_keywords,'Value',1);


% --- Executes on button press in pushbutton_help_find_by_key.
function pushbutton_help_find_by_key_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_find_by_key (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'find_by_key_window.html');
webcall(entry,'-helpbrowser');
