function varargout = select_reference(varargin)
% SELECT_REFERENCE M-file for select_reference.fig
%      SELECT_REFERENCE, by itself, creates a new SELECT_REFERENCE or raises the existing
%      singleton*.
%
%      H = SELECT_REFERENCE returns the handle to a new SELECT_REFERENCE or the handle to
%      the existing singleton*.
%
%      SELECT_REFERENCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_REFERENCE.M with the given input arguments.
%
%      SELECT_REFERENCE('Property','Value',...) creates a new SELECT_REFERENCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_reference_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_reference_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_reference

% Last Modified by GUIDE v2.5 12-Jan-2018 12:56:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_reference_OpeningFcn, ...
                   'gui_OutputFcn',  @select_reference_OutputFcn, ...
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


% --- Executes just before select_reference is made visible.
function select_reference_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_reference (see VARARGIN)

global hMain
global model

% Choose default command line output for select_reference
handles.output = hObject;

if isempty(hMain.current_reference),
    add_msg_board('ERROR: No references available for selection.');
    delete(hObject);
    return;
elseif length(hMain.current_reference)==1,
    add_msg_board('Warning: Reference is already selected.');
    delete(hObject);
    return;
end;

for k=1:length(hMain.current_reference),
    reference=model.references(hMain.current_reference(k));
    short=reference.short;
    pp=strfind(short,'_');
    if ~isempty(pp) && pp(1)>1,
        short=short(1:pp(1)-1);
    end;
    if isempty(short), short='???:yyyy'; end;
    if length(short)>11,
        short=short(1:11);
    end;
    list_item=sprintf('%s',right_pad([short ': '],13,1,1));
    author='';
    if ~isempty(reference.authors),
        authors=textscan(reference.authors,'%s','Delimiter',';');
        if ~isempty(authors{1}),
            author=strtrim(char(authors{1}(1)));
        end;
    end;
    if isempty(author), author='???'; end;
    if length(author)>11,
        author=author(1:11);
    end;
    list_item=sprintf('%s%s',list_item,right_pad([author '; '],13,1,1));
    journal='J. ???';
    if ~isempty(reference.journal),
        journal=reference.journal;
    end;
    if length(journal)>15,
        journal=journal(1:15);
    end;
    list_item=sprintf('%s%s',list_item,right_pad([journal '. '],17,1,1));
    year='YYYY';
    if ~isempty(reference.year),
        year=reference.year;
    end;
    if length(year)>6,
        year=year(1:6);
    end;
    list_item=sprintf('%s%s',list_item,right_pad([year ', '],8,1,1));
    vol='???';
    if ~isempty(reference.volume),
        vol=reference.volume;
    end;
    if length(vol)>4,
        vol=vol(1:4);
    end;
    list_item=sprintf('%s%s',list_item,right_pad([vol ', '],6,1,1));
    pages='???-???';
    if ~isempty(reference.pages),
        pages=reference.pages;
    end;
    if length(pages)>9,
        pages=pages(1:9);
    end;
    list_item=sprintf('%s%s',list_item,right_pad([pages ': '],11,1,1));
    title='???';
    if ~isempty(reference.title),
        title=reference.title;
    end;
    if length(title)>29,
        title=[title(1:29) '...'];
    end;
    list_item=sprintf('%s%s',list_item,title);
    list{k}=list_item;
end;

[list,poi]=sort(list);
set(handles.listbox_references,'String',list);
set(handles.listbox_references,'UserData',poi);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_reference wait for user response (see UIRESUME)
uiwait(handles.my_figure);


% --- Outputs from this function are returned to the command line.
function varargout = select_reference_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on selection change in listbox_references.
function listbox_references_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_references (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_references contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_references


% --- Executes during object creation, after setting all properties.
function listbox_references_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_references (see GCBO)
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

global hMain

sel=get(handles.listbox_references,'Value');
poi=get(handles.listbox_references,'UserData');
sel=poi(sel);

hMain.current_reference=hMain.current_reference(sel);

delete(handles.my_figure);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.current_reference=[];

delete(handles.my_figure);

% --- Executes when user attempts to close my_figure.
function my_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to my_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global hMain

hMain.current_reference=[];

delete(hObject);


function str=right_pad(str0,len,cut,lefttrim)
% right pads a string with white space to length len
% if the input string is longer than len, nothing happens unless argument
% cut is provided and is nonzero
% if argument lefftrim is given and nonzero, left spaces are removed

if nargin<3,
    cut=0;
end;

if nargin>3 && lefttrim~=0,
    str0=strtrim(str0);
end;

str=str0;
if length(str0)<len,
    for k=1+length(str0):len,
        str=[str ' '];
    end;
end;

if cut,
    str=str(1:len);
end;
