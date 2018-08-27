function varargout = select_target(varargin)
% SELECT_TARGET M-file for select_target.fig
%      SELECT_TARGET, by itself, creates a new SELECT_TARGET or raises the existing
%      singleton*.
%
%      H = SELECT_TARGET returns the handle to a new SELECT_TARGET or the handle to
%      the existing singleton*.
%
%      SELECT_TARGET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_TARGET.M with the given input arguments.
%
%      SELECT_TARGET('Property','Value',...) creates a new SELECT_TARGET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_target_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_target_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_target

% Last Modified by GUIDE v2.5 23-Mar-2011 10:52:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_target_OpeningFcn, ...
                   'gui_OutputFcn',  @select_target_OutputFcn, ...
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


% --- Executes just before select_target is made visible.
function select_target_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_target (see VARARGIN)

% Choose default command line output for select_target
handles.output = hObject;

global hMain

if ~isfield(hMain,'container'),
    add_msg_board('ERROR: No input data for target sequence selection.');
    guidata(hObject, handles);
    delete(hObject);
    return
elseif ~isfield(hMain.container,'alignment') || ~isfield(hMain.container,'template_id'),
    add_msg_board('ERROR: Wrong input data for target sequence selection.');
    guidata(hObject, handles);
    delete(hObject);
    return
else    
    handles.alignment=hMain.container.alignment;
    handles.template_id=hMain.container.template_id;
end;

% Treat the trivial cases
if length(handles.alignment)<2,
    add_msg_board('Warning: Only one sequence in alignment. Target and template sequence are the same.');
    if isfield(hMain,'container'),
        hMain=rmfield(hMain,'container');
    end;
    hMain.container.alg=1;
    hMain.container.identity=1;
    delete(hObject);
    return
end;

if length(handles.alignment)==2,
    add_msg_board('Only two sequences in alignment. Automatic choice of target sequence.');
    if isfield(hMain,'container'),
        hMain=rmfield(hMain,'container');
    end;
    if handles.template_id==1,
        hMain.container.alg=2;
    else
        hMain.container.alg=1;
    end;
    ali(1)=alignment(hMain.container.alg);
    ali(2)=alignment(handles.template_id);
    hMain.container.identity=get_identity(ali);
    delete(hObject);
    return
end;

set(handles.text_template,'String',sprintf('%s (%s)',...
    handles.alignment(handles.template_id).name,...
    handles.alignment(handles.template_id).db));
handles=set_listbox(handles);

% Update handles structure
guidata(hObject, handles);

%UIWAIT makes select_target wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_target_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on selection change in listbox_target.
function listbox_target_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_target contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_target

handles=update_sequence(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_target_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_target (see GCBO)
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

if isfield(hMain,'container'),
    hMain=rmfield(hMain,'container');
end;
hMain.container.alg=handles.target_info.alg;
hMain.container.identity=handles.target_info.identity;
delete(handles.figure1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if isfield(hMain,'container'),
    hMain=rmfield(hMain,'container');
end;
hMain.container.alg=[];
hMain.container.identity=[];
delete(handles.figure1);



function handles=set_listbox(handles)

pointers=zeros(1,length(handles.alignment)-1);
poi=0;
for k=1:length(handles.alignment),
    if handles.template_id~=k,
        poi=poi+1;
        pointers(poi)=k;
        list{poi}=sprintf('%s (%s)',handles.alignment(k).name,handles.alignment(k).db);
    end;
end;
set(handles.listbox_target,'String',list);
set(handles.listbox_target,'Value',1);
set(handles.listbox_target,'UserData',pointers);
handles=update_sequence(handles);

function handles=update_sequence(handles)

sel=get(handles.listbox_target,'Value');
pointers=get(handles.listbox_target,'UserData');
target_id=pointers(sel);
handles.target_info.alg=target_id;
ali(1)=handles.alignment(target_id);
ali(2)=handles.alignment(handles.template_id);
set(handles.text_sequence,'String',ali(1).sequence);
[ident,len]=get_identity(ali);
handles.target_info.identity=ident;
set(handles.text_identity,'String',sprintf('%4.1f%%',100*ident));
if ident>0.4,
    col=[0,1,0];
elseif ident<0.2,
    col=[1,0,0];
else
    twilight=(ident-0.2)/0.2;
    col=twilight*[0,1,0]+(1-twilight)*[1,0,0];
end;
set(handles.text_identity,'ForegroundColor',col);
set(handles.text_sequence_length,'String',sprintf('%i',len));
    
function [ident,len]=get_identity(ali)
% in an alignement array with exactly two sequences, get_identity
% determines relative sequence identity between the sequences normalized to
% the length of the first sequence
% for an alignment of more than two sequences, identity of zero is returned
% identity of zero is also returned if sequence lengths do not match (not
% aligned) or sequences are empty
% the length len of the first sequence is also returned

ident=0;
len=0;
if length(ali)~=2,
    return
end;
seq1=ali(1).sequence;
seq2=ali(2).sequence;
if length(seq1)~=length(seq2),
    return
end;
len=0;
matches=0;
for k=1:length(seq1),
    if char(seq1(k))~='-', % no gap in sequence 1
        len=len+1; 
        if char(seq1(k))==char(seq2(k)),
            matches=matches+1;
        end;
    end;
end;
if len>0,
    ident=matches/len;
end;
