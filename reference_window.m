function varargout = reference_window(varargin)
% REFERENCE_WINDOW M-file for reference_window.fig
%      REFERENCE_WINDOW, by itself, creates a new REFERENCE_WINDOW or raises the existing
%      singleton*.
%
%      H = REFERENCE_WINDOW returns the handle to a new REFERENCE_WINDOW or the handle to
%      the existing singleton*.
%
%      REFERENCE_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REFERENCE_WINDOW.M with the given input arguments.
%
%      REFERENCE_WINDOW('Property','Value',...) creates a new REFERENCE_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reference_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reference_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reference_window

% Last Modified by GUIDE v2.5 19-Dec-2014 09:48:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reference_window_OpeningFcn, ...
                   'gui_OutputFcn',  @reference_window_OutputFcn, ...
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


% --- Executes just before reference_window is made visible.
function reference_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reference_window (see VARARGIN)

% Choose default command line output for reference_window

global model
global hMain
global MMM_icon
global dummy_reference
global reference_formats

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.output = hObject;

if hMain.current_reference==0,
    if isfield(model,'references') && isempty(model.references),
        model = rmfield(model,'references');
    end;
    if isfield(model,'references'),
        handles.refnum=length(model.references)+1;
    else
        handles.refnum=1;
    end;
    handles.current.format=0;
    handles.current.type=1;
    handles.current.title='';
    handles.current.book_title='';
    handles.current.authors='';
    handles.current.editors='';
    handles.current.journal='';
    handles.current.year='';
    handles.current.volume='';
    handles.current.issue='';
    handles.current.chapter='';
    handles.current.pages='';
    handles.current.DOI='';
    handles.current.publisher='';
    handles.current.city='';
    handles.current.short=sprintf('[%i]',handles.refnum);
    handles.current.URL='';
    handles.current.PMID=0;
    model.references(handles.refnum)=handles.current;
else
    handles.refnum=hMain.current_reference;
    handles.current=model.references(handles.refnum);
end;

handles.old=model.references;

handles=update_page(handles);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.reference_open=hObject;

hMain.auxiliary=[hMain.auxiliary hObject];

if ~isempty(hMain.annotation_open), 
    set(hObject,'WindowStyle','modal');
else
    set(hObject,'WindowStyle','normal');
end;

axes(handles.axes_example);
handles.example=text(0,0.5,'.');
set(handles.example,'Interpreter','Latex','FontSize',10,'FontName','Bookman','Color',[0,0,0.5]);
set(handles.example,'String','A. First, B. Second, ... and C. Last \textit{J. Format. Test.} \textbf{2009}, \textit{123}, 456--78. \textit{A test title...}');

def_reference_formats;

handles=setup_journals(handles);

sel=get(handles.listbox_journal,'Value');
formats=get(handles.listbox_journal,'UserData');
format=reference_formats(formats(sel));

output=format_reference(dummy_reference,format,'Latex');

axes(handles.axes_example);
cla;
handles.example=text(0,0.5,output);
set(handles.example,'Interpreter','Latex','FontSize',10,'Color',[0,0,0.5]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes reference_window wait for user response (see UIRESUME)
if ~isempty(hMain.annotation_open), 
    uiwait(handles.references);
end;

% --- Outputs from this function are returned to the command line.
function varargout = reference_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on selection change in popupmenu_type.
function popupmenu_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_type

handles.current.type=get(hObject,'Value');
handles=update_page(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_import.
function pushbutton_import_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general

my_path=pwd;
cd(general.biblio_files);

contents = get(handles.popupmenu_format,'String');
mode=contents{get(handles.popupmenu_format,'Value')};

switch mode
    case 'ISI CGI (EndNote)'
        ext='.ciw';
        title='Load references in ISI CGI format';
        import=1;
    case 'SciFinder tagged format'
        ext='.txt';
        title='Load references in SciFinder tagged format';
        import=2;
    case 'EndNote text export' 
        ext='.txt';
        title='Load references in EndNote text export format';
        import=3;
    case 'BibTeX'
        title='Load references in BibTex format';
        ext='.bib';
        import=4;
    case 'MEDLINE'
        title='Load references in MEDLINE (PubMed) format';
        ext='.txt';
        import=5;
    case 'MMM'
        title='Load references in MMM format';
        ext='.mat';
        import=6;
end;

[filename,pathname] = uigetfile(ext,title);
if pathname~=0,
    reset_user_paths(pathname);
    general.biblio_files=pathname;
end;

if isempty(model.references),
    poi=0;
else
    poi=length(model.references);
end;

if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Loading of references canceled by user');
else
    references=[];
    switch import
        case 1
            references=rd_ISI_CGI(fullfile(pathname,filename),poi);
        case 2
            references=rd_SciFinder(fullfile(pathname,filename),poi);
        case 3
            references=rd_EndNote(fullfile(pathname,filename),poi);
        case 4
            references=rd_BibTex(fullfile(pathname,filename),poi);
        case 5
            references=rd_MEDLINE(fullfile(pathname,filename),poi);
        case 6
            if exist('references','var'), 
                clear references 
            end;
            load(fullfile(pathname,filename));
            if ~exist('references','var'), 
                references=[];
            end;
    end;
    m=length(references);
    if m>0,
        if isempty(model.references),
            model.references=references;
        else
            poi=length(model.references);
            newrefs=0;
            oldrefs=0;
            for k=1:m, % add only new references, as judged by PubMed ID (if available)
                found=0;
                if references(k).PMID~=0,                   
                    for kk=1:poi,
                        if model.references(kk).PMID==references(k).PMID,
                            found=1; break;
                        end;
                    end;
                end;
                if ~found,
                    poi=poi+1;
                    newrefs=newrefs+1;
                    model.references(poi)=references(k);
                else
                    oldrefs=oldrefs+1;
                    references(k).short=model.references(kk).short;
                    if isempty(references(k).DOI),
                        references(k).DOI=model.references(kk).DOI;
                    end;
                    model.references(kk)=references(k);
                end;
            end;
        end;
        add_msg_board(sprintf('Imported %i reference(s) with new or without any PubMed IDs',newrefs)); 
        add_msg_board(sprintf('Updated %i reference(s) with already existing PubMed IDs',oldrefs)); 
    else
        add_msg_board('Empty reference file or wrong format.');
    end;
    if ~isempty(model.references),
        handles.refnum=length(model.references);
        handles.current=model.references(handles.refnum);
        handles=update_page(handles);
    end;
end;
cd(my_path);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_previous_reference.
function pushbutton_previous_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

model.references(handles.refnum)=handles.current;

if handles.refnum>1,
    handles.refnum=handles.refnum-1;
    handles.current=model.references(handles.refnum);
end;

handles=update_page(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_next_reference.
function pushbutton_next_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

model.references(handles.refnum)=handles.current;

n=length(model.references);
if handles.refnum<n,
    handles.current=model.references(handles.refnum+1);
else
    handles.current.format=0;
    handles.current.type=1;
    handles.current.title='';
    handles.current.book_title='';
    handles.current.authors='';
    handles.current.editors='';                
    handles.current.journal='';                
    handles.current.year='';
    handles.current.volume='';
    handles.current.issue='';
    handles.current.chapter='';
    handles.current.pages='';
    handles.current.DOI='';
    handles.current.publisher='';
    handles.current.city='';
    handles.current.short=sprintf('[%i]',handles.refnum);
    handles.current.URL='';
    handles.current.PMID=0;
    model.references(handles.refnum+1)=handles.current;
end;
handles.refnum=handles.refnum+1;
handles=update_page(handles);
guidata(hObject,handles);


function edit_authors_Callback(hObject, eventdata, handles)
% hObject    handle to edit_authors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_authors as text
%        str2double(get(hObject,'String')) returns contents of edit_authors as a double

handles.current.authors=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_authors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_authors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_title_Callback(hObject, eventdata, handles)
% hObject    handle to edit_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_title as text
%        str2double(get(hObject,'String')) returns contents of edit_title as a double

handles.current.title=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_year_Callback(hObject, eventdata, handles)
% hObject    handle to edit_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_year as text
%        str2double(get(hObject,'String')) returns contents of edit_year as a double

handles.current.year=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_volume_Callback(hObject, eventdata, handles)
% hObject    handle to edit_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_volume as text
%        str2double(get(hObject,'String')) returns contents of edit_volume as a double

handles.current.volume=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_issue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_issue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_issue as text
%        str2double(get(hObject,'String')) returns contents of edit_issue as a double

if handles.current.type==3,
    handles.current.chapter=get(hObject,'String');
else
    handles.current.issue=get(hObject,'String');
end;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_issue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_issue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pages_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pages as text
%        str2double(get(hObject,'String')) returns contents of edit_pages as a double

handles.current.pages=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_pages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_DOI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DOI as text
%        str2double(get(hObject,'String')) returns contents of edit_DOI as a double

handles.current.DOI=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_DOI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_editors_Callback(hObject, eventdata, handles)
% hObject    handle to edit_editors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_editors as text
%        str2double(get(hObject,'String')) returns contents of edit_editors as a double

if handles.current.type<4,
    handles.current.editors=get(hObject,'String');
else
    handles.current.URL=get(hObject,'String');
end;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_editors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_editors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_publisher_Callback(hObject, eventdata, handles)
% hObject    handle to edit_publisher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_publisher as text
%        str2double(get(hObject,'String')) returns contents of edit_publisher as a double

handles.current.publisher=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_publisher_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_publisher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_city_Callback(hObject, eventdata, handles)
% hObject    handle to edit_city (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_city as text
%        str2double(get(hObject,'String')) returns contents of edit_city as a double

handles.current.city=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_city_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_city (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
global hMain

model.references(handles.refnum)=handles.current;
hMain.current_reference=handles.refnum;
references_CloseRequestFcn(handles.references, eventdata, handles);
%delete(handles.references);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

button = questdlg('All edits and loaded references will be lost','Are you sure?','No','Yes','No');

if strcmp(button,'Yes'),
    model.references=handles.old;
    hMain.current_reference=0;
    references_CloseRequestFcn(handles.references, eventdata, handles);
%     delete(handles.references);
end;


function edit_short_Callback(hObject, eventdata, handles)
% hObject    handle to edit_short (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_short as text
%        str2double(get(hObject,'String')) returns contents of edit_short as a double

global model

new_short=get(hObject,'String');
exists=0;
if isfield(model,'references'),
    for k=1:length(model.references),
        if strcmp(model.references(k).short,new_short),
            exists=1; break;
        end;
    end;
end;

if exists,
    msgbox('Please use other short name','This short name is already in use.','error');
    set(hObject,'String',handles.current.short);
else
    handles.current.short=new_short;
end;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_short_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_short (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_book_title_Callback(hObject, eventdata, handles)
% hObject    handle to edit_book_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_book_title as text
%        str2double(get(hObject,'String')) returns contents of edit_book_title as a double

handles.current.book_title=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_book_title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_book_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=update_page(handles)

set(handles.edit_reference_number,'String',sprintf('%i',handles.refnum));
type=handles.current.type;
if type==0,
    type=1;
end;
set(handles.popupmenu_type,'Value',type);
set(handles.edit_title,'String',handles.current.title);
set(handles.edit_book_title,'String',handles.current.book_title);
set(handles.edit_authors,'String',handles.current.authors);
if type~=4,
    set(handles.text_editors_or_URL,'String','Editors');
    set(handles.pushbutton_URL_web,'Enable','off');
    set(handles.edit_editors,'String',handles.current.editors);
else
    set(handles.text_editors_or_URL,'String','URL');
    set(handles.pushbutton_URL_web,'Enable','on');
    set(handles.edit_editors,'String',handles.current.URL);
end;
if isempty(handles.current.editors) && ~isempty(handles.current.URL),
    set(handles.text_editors_or_URL,'String','URL');
    set(handles.pushbutton_URL_web,'Enable','on');
    set(handles.edit_editors,'String',handles.current.URL);
end;
set(handles.edit_journal,'String',handles.current.journal);
set(handles.edit_year,'String',handles.current.year);
set(handles.edit_volume,'String',handles.current.volume);
set(handles.edit_issue,'String',handles.current.issue);
if type==3,
    set(handles.edit_issue,'String',handles.current.chapter);
end;
if type==8,
    set(handles.text_volume,'String','Country');
    set(handles.text_chapter,'String','Patent No.');
else
    set(handles.text_volume,'String','Volume');
    set(handles.text_chapter,'String','Chapter/Issue');
end;
set(handles.edit_pages,'String',handles.current.pages);
set(handles.edit_DOI,'String',handles.current.DOI);
set(handles.edit_publisher,'String',handles.current.publisher);
set(handles.edit_city,'String',handles.current.city);
set(handles.edit_short,'String',handles.current.short);
% if isfield(handles.current,'PMID') && handles.current.PMID>0,
%     set(handles.pushbutton_abstract,'Enable','on');
% else
%     set(handles.pushbutton_abstract,'Enable','off');
% end;
set(handles.pushbutton_abstract,'Enable','on');

% --- Executes on button press in pushbutton_DOI_web.
function pushbutton_DOI_web_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DOI_web (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global queries

preamble=web_adr.DOI_resolver;

if ~isempty(handles.current.DOI),
    [rem,suffix]=strtok(handles.current.DOI,':');
    if strcmpi(strtrim(rem),'doi');
        URL=strcat(preamble,strtrim(suffix(2:end)));
        webcall(URL);
    else
        URL=strcat(preamble,strtrim(handles.current.DOI));
        webcall(URL);
    end;
else
    add_msg_board('No DOI available. Trying SFX.');
    query=[web_adr.SFX queries.SFX_preamble];
    length_empty=length(query);
    if ~isempty(strtrim(handles.current.title)),
        query=sprintf('%s%s%s',query,queries.SFX_title,SFX_string(strtrim(handles.current.title)));
    end;
    if ~isempty(handles.current.journal),
        query=sprintf('%s%s%s',query,queries.SFX_journal,SFX_string(strtrim(upper(handles.current.journal))));
        % the following fixes SFX's inability to find out whether
        % Biochemistry (Moscow) or Biochemistry of ACS has this reference,
        % it does so, by checking only ACS (sorry, I am not anti-Russian, but ACS is needed much more frequently)
        if strcmpi(handles.current.journal,'Biochemistry'),
            query=sprintf('%s&rft.issn=0006-2960',query);
        end;
    end;
    if ~isempty(handles.current.year),
        query=sprintf('%s%s%s',query,queries.SFX_year,strtrim(handles.current.year));
    end;
    if ~isempty(handles.current.volume),
        query=sprintf('%s%s%s',query,queries.SFX_volume,strtrim(handles.current.volume));
    end;
    if ~isempty(handles.current.pages)
        poi=strfind(handles.current.pages,'-');
        if isempty(poi) || length(poi)>1 || poi<=1 || poi==length(handles.current.pages),
            apage=strtrim(handles.current.pages);
            epage=[];
        else
            apage=strtrim(handles.current.pages(1:poi-1));
            epage=strtrim(handles.current.pages(poi+1:end));
        end;
        if ~isempty(apage),
            query=sprintf('%s%s%s',query,queries.SFX_first_page,apage);
        end;
    end;
    authors='';
    if ~isempty(handles.current.authors),
        authors=textscan(handles.current.authors,'%s','delimiter',';');
    end;
    if ~isempty(authors),
        if ~isempty(authors{1}),
            for k=2:length(authors{1}),
                author{k}=strtrim(char(authors{1}(k)));
            end;
            [surname,initials]=strtok(author{1});
            surname=strtrim(surname);
            initials=strtrim(initials);
            if ~isempty(initials),
                query=sprintf('%s%s%s',query,queries.SFX_initial,initials(1));
            end;
            if ~isempty(surname),
                query=sprintf('%s%s%s',query,queries.SFX_surname,surname);
            end;
            if length(author)>1,
                for k=1:length(author),
                    [surname,initials]=strtok(author{k});
                    surname=strtrim(surname);
                    initials=strtrim(initials);
                    if ~isempty(initials),
                        aq=[surname queries.SFX_sep initials(1)];
                    else
                        aq=surname;
                    end;
                    if ~isempty(aq),
                        query=sprintf('%s%s%s',query,queries.SFX_author,aq);
                    end;
                end;
            end;
        end;
    end;
    % a working query for test purposes:
    % query='sfx.ethz.ch:9003/sfx_locater?&url_ver=Z39.88-2004&url_ctx_fmt=info:ofi/fmt:kev:mtx:ctx&rft_val_fmt=info:ofi/fmt:kev:mtx:journal&rft.atitle=Isolation%20and%20characterization%20of%20lamellar%20aggregates%20of%20LHCII%20and%20LHCII-lipid%20macro-assemblies%20with%20light-inducible%20structural%20transitions&rft.auinit=I&rft.aulast=Simidjiev&rft.btitle=METHODS%20MOL%20BIOL&rft.date=2004&rft.spage=105&rft.volume=274&rft.au=Varkonyi%2C%20Z&rft.au=Garab%2C%20G';
    if length(query)>length_empty,
        webcall(query);
    else
        add_msg_board('Not enough information for SFX query.');
    end;
end;


% --- Executes on button press in pushbutton_URL_web.
function pushbutton_URL_web_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_URL_web (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.current.URL),
    webcall(handles.current.URL);
end;

% --- Executes during object creation, after setting all properties.
function pushbutton_DOI_web_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_DOI_web (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

load web_button
set(hObject,'CData',icon_data);


% --- Executes during object creation, after setting all properties.
function pushbutton_URL_web_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_URL_web (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

load web_button
set(hObject,'CData',icon_data);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'reference_window.html');
webcall(entry,'-helpbrowser');


% --- Executes on selection change in popupmenu_copy_format.
function popupmenu_format_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_copy_format contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_format


% --- Executes on selection change in popupmenu_copy_format.
function popupmenu_copy_format_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_copy_format contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_format


% --- Executes during object creation, after setting all properties.
function popupmenu_copy_format_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

linked=0;
if isfield(model,'annotations'),
    for k=1:length(model.annotations),
        refs=model.annotations(k).info.references;
        if ~isempty(refs),
            for kk=1:length(refs),
                if refs(kk)==handles.refnum,
                    linked=1;
                    adr=mk_address(model.annotations(k).indices);
                end;
            end;
        end;
    end;
end;
if isfield(hMain,'annotation_references'),
    if ~isempty(hMain.annotation_references),
        for k=1:length(hMain.annotation_references),
            if hMain.annotation_references(k)==handles.refnum,
                linked=1;
                adr='current object';
            end;
        end;
    end;
end;

n=length(model.references);
if linked,
    msgbox('Linked reference cannot be deleted. Please remove link first.',sprintf('This reference is used by %s.',adr),'warn');
    return;
elseif n==1,
    msgbox('At least one reference page must exist.','Cannot delete the only reference page','warn');
    return;
else
    answer=questdlg('Do you want to delete reference?','No undo available','Delete','Cancel','Delete');
    if strcmp(answer,'Delete'),
        if n>handles.refnum,
            for k=handles.refnum:n-1,
                model.references(k)=model.references(k+1);
            end;
            if isfield(model,'annotations'),
                for k=1:length(model.annotations), % renumber references in links
                    refs=model.annotations(k).info.references;
                    if ~isempty(refs),
                        for kk=1:length(refs),
                            if refs(kk)>handles.refnum,
                                refs(kk)=refs(kk)-1;
                            end;
                        end;
                    end;
                    model.annotations(k).info.references=refs;
                end;
            end;
            if isfield(hMain,'annotation_references'),
                if ~isempty(hMain.annotation_references),
                    for k=1:length(hMain.annotation_references),
                        if hMain.annotation_references(k)>handles.refnum,
                            hMain.annotation_references(k)=hMain.annotation_references(k)-1;
                        end;
                    end;
                end;
            end;
        end;
        model.references=model.references(1:n-1);
        if handles.refnum>length(model.references),
            handles.refnum=length(model.references);
        end;
        handles.current=model.references(handles.refnum);
    end;
end;

handles=update_page(handles);
guidata(hObject,handles);

function edit_journal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_journal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_journal as text
%        str2double(get(hObject,'String')) returns contents of edit_journal as a double
handles.current.journal=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_journal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_journal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_search_short.
function pushbutton_search_short_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_search_short (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

old_short=handles.current.short;
new_short=inputdlg('Enter short name:','Reference search by short name',1,{old_short});

found=0;
mynum=handles.refnum;
if isfield(model,'references'),
    for k=1:length(model.references),
        if strcmp(model.references(k).short,new_short),
            found=1;
            mynum=k;
            break;
        end;
    end;
end;
if found,
    handles.refnum=mynum;
    handles.current=model.references(handles.refnum);
else
    msgbox('Staying with old reference.','Short name not found.','error');
end;
handles=update_page(handles);
guidata(hObject,handles);

function edit_reference_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_reference_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_reference_number as text
%        str2double(get(hObject,'String')) returns contents of edit_reference_number as a double

global model

if isfield(model,'references'),
        maxnum=length(model.references);
else
    maxnum=1;
end;

[v,handles]=edit_update_MMM(handles,hObject,1,maxnum,handles.refnum,'%i',1);

model.references(handles.refnum)=handles.current;

handles.refnum=v;
handles.current=model.references(handles.refnum);
handles=update_page(handles);
guidata(hObject,handles);


% --- Executes when user attempts to close references.
function references_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to references (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global hMain

delete(hMain.reference_open);
hMain.reference_open=[];


% --- Executes on button press in pushbutton_search.
function pushbutton_search_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global queries
global model
global hMain

newrefs=0;
search_PubMed; % can be PubMed search or search in bibliography

if isempty(hMain.current_reference), % this was a PubMed search
    fname=strcat(general.tmp_files,queries.PubMed_file);
    if isempty(model.references),
        poi=0;
        oldnum=0;
    else
        poi=length(model.references);
        oldnum=poi;
    end;
    references=rd_MEDLINE(fname,poi);
    if isempty(references),
        add_msg_board('PubMed search does not return any references.');
    else
        m=length(references);
        add_msg_board(sprintf('PubMed search returned %i reference(s).',m));
        if m>0,
            if isempty(model.references),
                model.references=references;
                newrefs=length(references);
            else
                poi=length(model.references);
                newrefs=0;
                oldrefs=0;
                for k=1:m, % add only new references, as judged by PubMed ID (if available)
                    found=0;
                    if references(k).PMID~=0,                   
                        for kk=1:poi,
                            if model.references(kk).PMID==references(k).PMID,
                                found=1; break;
                            end;
                        end;
                    end;
                    if ~found,
                        poi=poi+1;
                        newrefs=newrefs+1;
                        model.references(poi)=references(k);
                    else
                        oldrefs=oldrefs+1;
                        references(k).short=model.references(kk).short;
                        if isempty(references(k).DOI),
                            references(k).DOI=model.references(kk).DOI;
                        end;
                        model.references(kk)=references(k);
                    end;
                end;
            end;
            add_msg_board(sprintf('Imported %i reference(s) with new or without any PubMed IDs',newrefs)); 
            add_msg_board(sprintf('Updated %i reference(s) with already existing PubMed IDs',oldrefs)); 
        end;
    end;

    if ~isempty(model.references),
        if newrefs>0
            handles.refnum=oldnum+1;
        end;
        handles.current=model.references(handles.refnum);
        handles=update_page(handles);
        hMain.current_reference=handles.refnum;
    end;
else % this was a (successful) search in the bibliography
    handles.refnum=hMain.current_reference;
    handles.current=model.references(handles.refnum);
    handles=update_page(handles);
end;

guidata(hObject,handles);
    
% --- Executes on button press in pushbutton_abstract.
function pushbutton_abstract_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_abstract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global web_adr
global queries

empty_ref=0;
if ~isempty(handles.current.PMID) && handles.current.PMID>0,
    query=sprintf('%s%i%s',web_adr.PubMed,handles.current.PMID,queries.PMID);
    fname=strcat(general.tmp_files,'PubMed.txt');
    [f,status]=urlwrite(query,fname);
else
    biblio='?term=';
    if ~isempty(handles.current.authors),
        authors=author_query(handles.current.authors);
        if ~isempty(authors),
            biblio=[biblio authors];
        end;
    end;
    if ~isempty(handles.current.year),
        years=strtrim(handles.current.year);
        for k=1:length(years),
            if strcmp(years(k),'-'),
                years(k)=':';
            end;
            if strcmp(years(k),' '),
                years(k)='+';
            end;
        end;
        if ~isempty(years),
            if length(biblio)>length('?terms='),
                biblio=[biblio '+AND+'];
            end;
            biblio=[biblio,years,'[dp]'];
        end;
    end;
    if ~isempty(handles.current.journal),
        journal=handles.current.journal;
        for k=1:length(journal),
            if strcmp(journal(k),'.'),
                journal(k)='+';
            end;
            if strcmp(journal(k),'-'),
                journal(k)=' ';
            end;
        end;
        journal=compact(journal);
        if ~isempty(journal),
            if length(biblio)>length('?terms='),
                biblio=[biblio '+AND+'];
            end;
            biblio=[biblio,journal,'[ta]'];
        end;
    end;
    if length(biblio)>length('?terms='),
        for k=1:length(biblio),
            if strcmp(biblio(k),' '),
                biblio(k)='+';
            end;
        end;
        query=sprintf('%s%s%s%s%s%s',web_adr.PubMed,biblio,'&dispmax=50',queries.search_PubMed,queries.tool,queries.email);
        fname=strcat(general.tmp_files,'PubMed.txt');
        [f,status]=urlwrite(query,fname);
    else
        add_msg_board('Not enough information to retrieve this reference.');
        status=0;
        empty_ref=1;
    end;
end;
if status,
    reference=rd_MEDLINE(f,1);
    if ~isempty(reference),
        page=abstract_MEDLINE(f,1);
        if ~isempty(page),
            webcall(strcat('file:///',page));
        else
            add_msg_board('ERROR: PubMed record did not contain abstract.');
        end;
        if length(reference)>1,
            add_msg_board('Warning: Initial search did not retrieve a unique reference.');
            Title='Several references match information.';
            Question='Accept the reference with the shown abstract for update?';
            for k=1:length(reference),
                ap0=strtok(handles.current.pages,'-');
                ap1=strtok(reference(k).pages,'-');
                if strcmp(ap0,ap1),
                    add_msg_board('Reference with matching first page found.');
                    v0=strtok(handles.current.volume);
                    v1=strtok(reference(k).volume);
                    if strcmp(v0,v1),
                        add_msg_board('Volume is also matching.');
                        if k~=1,
                            msgbox('Displayed abstract does not belong to matching reference. Please redisplay abstract.','Several references found.','warn');
                        end;
                        Title='Matching reference found.';
                        Question='Accept matching reference for update?';
                        reference=reference(k);
                        break
                    end;
                end;
            end;
            reference=reference(1);
        else
            Title='One reference matches information.';
            Question='Accept this reference for update?';
        end;
        ButtonName = questdlg(Question, Title, 'Yes', 'No', 'Yes');
        if strcmp(ButtonName,'Yes');
            handles.current.type=reference.type;
            handles.current.title=reference.title;
            handles.current.book_title=reference.book_title;
            handles.current.authors=reference.authors;
            handles.current.editors=reference.editors;
            handles.current.journal=reference.journal;
            handles.current.year=reference.year;
            handles.current.volume=reference.volume;
            handles.current.issue=reference.issue;
            handles.current.chapter=reference.chapter;
            handles.current.pages=reference.pages;
            if ~isempty(reference.DOI),
                handles.current.DOI=reference.DOI;
            end;
            handles.current.publisher=reference.publisher;
            handles.current.city=reference.city;
        end;
    else
        add_msg_board('No matching reference found in PubMed.');
    end;
else
    add_msg_board('ERROR: PubMed record could not be accessed.');
    if ~empty_ref,
        add_msg_board('Please check internet connection.');
    end;
end;

handles=update_page(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_export.
function pushbutton_export_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general

my_path=pwd;

contents = get(handles.popupmenu_format,'String');
mode=contents{get(handles.popupmenu_format,'Value')};

switch mode
    case 'ISI CGI (EndNote)'
        ext='.cgi';
        title='Save bibliography in ISI CGI format (for EndNote import)';
        success='Bibliography saved to ISI .cgi file:';
        export=1;
    case 'SciFinder tagged format'
        msgbox('Sorry. Please use another export format.','SciFinder tagged format not available for export','warn');
        return
    case 'EndNote text export' 
        msgbox('Sorry. Please use another export format.','EndNote text format not available for export','warn');
        return
    case 'BibTeX'
        title='Save bibliography in BibTex format';
        success='Bibliography saved to BibTex file:';
        ext='.bib';
        export=4;
    case 'MEDLINE'
        title='Save bibliography in MEDLINE format';
        success='Bibliography saved to MEDLINE file:';
        ext='.txt';
        export=5;
    case 'MMM'
        title='Save bibliography in MMM format';
        success='Bibliography saved to MMM file:';
        ext='.mat';
        export=6;
end;

cd(general.biblio_files);
[filename, pathname] = uiputfile(ext, title);
cd(my_path);
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Saving of bibliography cancelled by user');
else
    reset_user_paths(pathname);
    general.biblio_files=pathname;
    fname=fullfile(pathname, filename);
    switch export
        case 1
            message=wr_ISI_CGI(fname,model.references);
        case 4
            message=wr_BibTex(fname,model.references);
        case 5
            message=wr_MEDLINE(fname,model.references);
        case 6
            message.error=0;
            message.text='No error.';
            references=model.references;
            if ~isempty(references),
                save(fname,'references');
            else
                message.error=1;
                message.text='No references to save.';
            end;
    end;
    if message.error,
        add_msg_board(sprintf('ERROR: %s',message.error));
    else
        add_msg_board(success);
        add_msg_board(fname);
    end;
end;

function all_authors=author_query(authors0)

authors=textscan(authors0,'%s','delimiter',';');
all_authors='';
if ~isempty(authors),
    if ~isempty(authors{1}),
        for k=1:length(authors{1}),
            author0=strtrim(char(authors{1}(k)));
            author='';
            for kk=1:length(author0),
                if author0(kk)~='.' && author0(kk)~=',',
                    author=[author author0(kk)];
                end;
            end;                                
            [surname,initials]=strtok(author);
            surname=strtrim(surname);
            initials=strtrim(initials);
            if length(initials)>1,
                initials=initials(1); % PubMed stores two initials, but program should also work with full name
            end;
            if isempty(all_authors),
                all_authors=surname;
                if ~isempty(initials),
                    all_authors=[all_authors '+' initials];
                end;
            else
                all_authors=[all_authors '+AND+' surname];
                if ~isempty(initials),
                    all_authors=[all_authors '+' initials];
                end;
            end;
        end;
        all_authors=sprintf('%s[au]',all_authors);
    end;
end;

function str=SFX_string(str0)

global queries

str='';
for k=1:length(str0),
    if strcmp(str0(k),' '),
        str=[str queries.SFX_space];
    else
        str=[str str0(k)];
    end;
end;


% --- Executes on selection change in listbox_journal.
function listbox_journal_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_journal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_journal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_journal

handles=update_example(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_journal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_journal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global reference_formats
global formats

contents = get(handles.popupmenu_copy_format,'String');
copy_format=contents{get(handles.popupmenu_copy_format,'Value')};
[target,rem]=strtok(copy_format);
target=strtrim(target);

sel=get(handles.listbox_journal,'Value');
journals=get(handles.listbox_journal,'UserData');
format=reference_formats(journals(sel));

formatted=format_reference(handles.current,format,target);

if get(handles.checkbox_short,'Value'),
    switch lower(target)
        case {'rtf','latex','plain'}
            formatted=sprintf('[%s]\t%s',handles.current.short,formatted);
        case 'html'
            formatted=sprintf('<pre>[%s]\t%s</pre>',handles.current.short,formatted);
    end;
end;

if strcmpi(target,'rtf'),
    formatted=[formats.rtf_preamble formats.rtf_origin formatted '}'];
end;

clipboard('copy',formatted);

% --- Executes on key press with focus on listbox_journal and none of its controls.
function listbox_journal_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_journal (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

list=get(hObject,'String');
initial=double(upper(eventdata.Character));
if strcmp(eventdata.Modifier,'shift')
    sel=get(hObject,'Value');
    initial=double(upper(list{sel}(1)));
    found=0;
    for k=1:length(list),
        if initial<=double(upper(list{k}(1))),
            found=k;
            break;
        end;
    end;
    if found==0,
        return;
    end;
    ini2=double(upper(eventdata.Character));
    found2=found;
    for k=found:length(list),
        sstring=unblank(list{k});
        if ini2<=double(upper(sstring(2))),
            found2=k;
            break;
        end;
        if initial~=double(upper(sstring(1))),
            found2=found;
            break;
        end;
    end;
    found=found2;
else
    found=0;
    for k=1:length(list),
        if initial<=double(upper(list{k}(1))),
            found=k;
            break;
        end;
    end;
end;
if found>0,
    set(hObject,'Value',found);
end;
handles=update_example(handles);
guidata(hObject,handles);

function handles=setup_journals(handles)

global journals

formats=zeros(1,length(journals));
for k=1:length(journals),
    list{k}=journals(k).name;
    formats(k)=journals(k).format;
end;
[list,poi]=sort(list);
set(handles.listbox_journal,'String',list);
set(handles.listbox_journal,'UserData',formats(poi));

function handles=update_example(handles)

global dummy_reference
global reference_formats

sel=get(handles.listbox_journal,'Value');
formats=get(handles.listbox_journal,'UserData');
format=reference_formats(formats(sel));

output=format_reference(dummy_reference,format,'Latex');

axes(handles.axes_example);
cla;
handles.example=text(0,0.5,output);
set(handles.example,'Interpreter','Latex','FontSize',10,'Color',[0,0,0.5]);



% --- Executes during object creation, after setting all properties.
function popupmenu_format_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_save_current.
function pushbutton_save_current_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global reference_formats
global formats
global general

ext='.txt';
contents = get(handles.popupmenu_copy_format,'String');
copy_format=contents{get(handles.popupmenu_copy_format,'Value')};
[target,rem]=strtok(copy_format);
target=strtrim(target);

sel=get(handles.listbox_journal,'Value');
journals=get(handles.listbox_journal,'UserData');
format=reference_formats(journals(sel));

formatted=format_reference(handles.current,format,target);

if strcmpi(target,'latex') || strcmpi(target,'bibitem')
    ext='.tex';
end;

if get(handles.checkbox_short,'Value'),
    switch lower(target)
        case {'rtf','latex','plain'}
            formatted=sprintf('[%s]\t%s',handles.current.short,formatted);
    end;
end;

if strcmpi(target,'rtf'),
    formatted=[formats.rtf_preamble formats.rtf_origin formatted '}'];
    ext='.rtf';
end;

fname=[general.tmp_files 'reference' ext];
fid=fopen(fname,'w');
fprintf(fid,'%s\n',formatted);
fclose(fid);

% --- Executes on button press in pushbutton_save_all.
function pushbutton_save_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global reference_formats
global formats
global general
global model

contents = get(handles.popupmenu_copy_format,'String');
copy_format=contents{get(handles.popupmenu_copy_format,'Value')};
[target,rem]=strtok(copy_format);
target=strtrim(target);

if strcmpi(target,'rtf'),
    ext='.rtf';
elseif strcmpi(target,'latex'),
    ext='.tex';
elseif strcmpi(target,'bibitem'),
    ext='.tex';
elseif strcmpi(target,'html'),
    ext='.html';
else
    ext='.txt';
end;

fname=[general.tmp_files 'bibliography' ext];

currdir=pwd;
cd(general.biblio_files);
[FileName,PathName,FilterIndex] = uiputfile(['*' ext],sprintf('Save bibliography in %s format',copy_format),fname);
reset_user_paths(PathName);
general.biblio_files=PathName;
fname=fullfile(PathName,FileName);

cd(currdir);

if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Saving of bibliography in a text format cancelled by user.');
    return;
end;


sel=get(handles.listbox_journal,'Value');
journals=get(handles.listbox_journal,'UserData');
format=reference_formats(journals(sel));

fid=fopen(fname,'w');
if fid==-1,
    add_msg_board('ERROR. File could not be opened for writing. File name:');
    add_msg_board(fname);
else
    if strcmpi(target,'rtf'),
        fprintf(fid,'%s%s\n',formats.rtf_preamble,formats.rtf_origin);
    end;

    for k=1:length(model.references),
        if get(handles.checkbox_short,'Value'),
            switch lower(target)
                case {'rtf','latex','plain','html'}
                    fprintf(fid,'[%s]\t',model.references(k).short);
            end;
        end;
        if get(handles.checkbox_number,'Value'),
            switch lower(target)
                case {'rtf','latex','plain'}
                    fprintf(fid,'[%i]\t',k);
                case 'html'
                    fprintf(fid,'<pre>[%i]\t</pre>',k);
            end;
        end;
        refline=format_reference(model.references(k),format,target);
        if strcmpi(target,'rtf'),
            refline=[refline '\par '];
        end;
        fprintf(fid,'%s\n',refline);
        if strcmpi(target,'latex'),
            fprintf(fid,'\n');
        end;
        if strcmpi(target,'html'),
            fprintf(fid,'<p>');
        end;
    end;
    if strcmpi(target,'rtf'),
        fprintf(fid,'}');
    end;
    fclose(fid);
end;
guidata(hObject,handles);

% --- Executes on button press in checkbox_short.
function checkbox_short_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_short (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_short

if get(hObject,'Value'),
    set(handles.checkbox_number,'Value',0);
end;

guidata(hObject,handles);

% --- Executes on button press in checkbox_number.
function checkbox_number_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_number


if get(hObject,'Value'),
    set(handles.checkbox_short,'Value',0);
end;

guidata(hObject,handles);
