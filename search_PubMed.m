function varargout = search_PubMed(varargin)
% SEARCH_PUBMED M-file for search_PubMed.fig
%      SEARCH_PUBMED, by itself, creates a new SEARCH_PUBMED or raises the existing
%      singleton*.
%
%      H = SEARCH_PUBMED returns the handle to a new SEARCH_PUBMED or the handle to
%      the existing singleton*.
%
%      SEARCH_PUBMED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEARCH_PUBMED.M with the given input arguments.
%
%      SEARCH_PUBMED('Property','Value',...) creates a new SEARCH_PUBMED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before search_PubMed_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to search_PubMed_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help search_PubMed

% Last Modified by GUIDE v2.5 12-Jan-2018 12:55:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @search_PubMed_OpeningFcn, ...
                   'gui_OutputFcn',  @search_PubMed_OutputFcn, ...
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


% --- Executes just before search_PubMed is made visible.
function search_PubMed_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to search_PubMed (see VARARGIN)

global MMM_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

% Choose default command line output for search_PubMed
handles.output = hObject;

load helpicon
set(handles.pushbutton_help,'CData',cdata);

handles=make_keyword_list(handles);

handles=setup_autosearch(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes search_PubMed wait for user response (see UIRESUME)
uiwait(handles.my_figure);


% --- Outputs from this function are returned to the command line.
function varargout = search_PubMed_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'search_PubMed.html');
webcall(entry,'-helpbrowser');

function edit_search_terms_Callback(hObject, eventdata, handles)
% hObject    handle to edit_search_terms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_search_terms as text
%        str2double(get(hObject,'String')) returns contents of edit_search_terms as a double


% --- Executes during object creation, after setting all properties.
function edit_search_terms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_search_terms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over link_PubMed.
function link_PubMed_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to link_PubMed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr

webcall(web_adr.PubMed_interactive);



function edit_authors_Callback(hObject, eventdata, handles)
% hObject    handle to edit_authors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_authors as text
%        str2double(get(hObject,'String')) returns contents of edit_authors as a double


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


% --- Executes on selection change in listbox_keywords.
function listbox_keywords_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_keywords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
handles.current_key=lower(contents{get(hObject,'Value')});
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


% --- Executes on button press in pushbutton_keyword_only.
function pushbutton_keyword_only_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_keyword_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.edit_search_terms,'String',handles.current_key);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_plus.
function pushbutton_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

old_terms=get(handles.edit_search_terms,'String');
if ~isempty(old_terms),
    old_terms=[old_terms ';'];
end;
set(handles.edit_search_terms,'String',[old_terms handles.current_key]);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_query.
function pushbutton_query_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_query (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu_records.
function popupmenu_records_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_records (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_records contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_records


% --- Executes during object creation, after setting all properties.
function popupmenu_records_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_records (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_years_Callback(hObject, eventdata, handles)
% hObject    handle to edit_years (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_years as text
%        str2double(get(hObject,'String')) returns contents of edit_years as a double

years=get(hObject,'String');
years=strtrim(years);
if ~isempty(years),
    [syear,eyear]=strtok(years,'-');
    eyear=strtrim(eyear);
    syear=strtrim(syear);
    if strcmp(eyear,'-'),
        today=date;
        eyear=today(end-3:end);
        set(hObject,'String',[syear '-' eyear]);
    end;
    if isempty(eyear) && strcmp(years(1),'-'),
        set(hObject,'String',['1900-' syear]);
    end;
end;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_years_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_years (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_journal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_journal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_journal as text
%        str2double(get(hObject,'String')) returns contents of edit_journal as a double


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


% --- Executes on button press in pushbutton_search.
function pushbutton_search_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

if get(handles.checkbox_bibliography,'Value'),
    ref=search_bibliography(handles);
    if isempty(ref),
        status=0;
        add_msg_board('No matching reference found in bibliography.');
    else
        if length(ref)==1,
            add_msg_board('One matching reference found in bibliography.');
            hMain.current_reference=ref;
            status=1;
        else
            add_msg_board('Several matching references found in bibliography.');
            add_msg_board('Please select one.');
            hMain.current_reference=ref;
            select_reference;
            if isempty(hMain.current_reference),
                status=0;
            else
                status=1;
            end;
        end;
    end;
else
    hMain.current_reference=[];

    [handles,status]=do_search(handles);
end;

if status>0,
    delete(handles.my_figure);
end;

function [handles,status]=do_search(handles)

global web_adr
global queries
global general

search_string='?term=';

term_mode='[all]';
if get(handles.radiobutton_TIAB,'Value'),
    term_mode='[tiab]';
end;
if get(handles.radiobutton_keywords,'Value'),
    term_mode='[mh]';
end;
if get(handles.radiobutton_major,'Value'),
    term_mode='[majr]';
end;
terms=get(handles.edit_search_terms,'String');
[ma,na]=size(terms); % handle multi-line input
if ma>1,
    nterms=terms(1,:);
    for ka=2:ma,
        nterms=[nterms terms(ka,:)];
    end;
    terms=nterms;
end;
terms=compact(terms);
poi=strfind(terms,';');
if ~isempty(poi),
    new_terms='';
    poi0=1;
    for k=1:length(poi),
        new_terms=[new_terms terms(poi0:poi(k)-1) term_mode '+AND+'];
        poi0=poi(k)+1;
    end;
    terms=[new_terms terms(poi0:end)];
end;
if ~isempty(terms),
    search_string=[search_string,terms,term_mode];
end;

author_mode='[au]';
if get(handles.radiobutton_first_author,'Value'),
    author_mode='[1au]';
end;
if get(handles.radiobutton_last_author,'Value'),
    author_mode='[lastau]';
end;
authors=get(handles.edit_authors,'String');
[ma,na]=size(authors); % handle multi-line input
if ma>1,
    nauthors=authors(1,:);
    for ka=2:ma,
        nauthors=[nauthors authors(ka,:)];
    end;
    authors=nauthors;
end;
for k=1:length(authors),
    if strcmp(authors(k),'.'),
        authors(k)=' ';
    end;
    if strcmp(authors(k),','),
        authors(k)=' ';
    end;
end;
poi=strfind(authors,';');
if ~isempty(poi),
    new_authors='';
    poi0=1;
    for k=1:length(poi),
        new_authors=[new_authors authors(poi0:poi(k)-1) '[au]+AND+'];
        poi0=poi(k)+1;
    end;
    authors=[new_authors authors(poi0:end)];
end;
authors=compact(authors);
if ~isempty(authors),
    if length(search_string)>length('?term='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,authors,author_mode];
end;

journal=get(handles.edit_journal,'String');
for k=1:length(journal),
    if strcmp(journal(k),'.'),
        journal(k)=' ';
    end;
end;
journal=compact(journal);
if ~isempty(journal),
    if length(search_string)>length('?term='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,journal,'[ta]'];
end;

years=get(handles.edit_years,'String');
for k=1:length(years),
    if strcmp(years(k),'-'),
        years(k)=':';
    end;
end;
years=unblank(years);
if ~isempty(years),
    if length(search_string)>length('?term='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,years,'[dp]'];
end;

if get(handles.checkbox_reviews,'Value'),
    if length(search_string)>length('?term='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,'Review[pt]'];
end;
contents = get(handles.popupmenu_records,'String');
records=contents{get(handles.popupmenu_records,'Value')};
dispmax=sprintf('&dispmax=%s',strtrim(records));

for k=1:length(search_string),
    if strcmp(search_string(k),' '),
        search_string(k)='+';
    end;
end;

query=sprintf('%s%s%s%s',web_adr.PubMed,search_string,queries.search_PubMed,dispmax,queries.tool,queries.email);
fname=strcat(general.tmp_files,queries.PubMed_file);
set(gcf,'Pointer','watch');
drawnow
[f,status] = urlwrite(query,fname);
set(gcf,'Pointer','arrow');
if status==0,
    add_msg_board('ERROR: PubMed access failed.');
    add_msg_board('Please check internet connection.');
    add_msg_board('Use "Cancel" button if connection cannot be repaired.');
    msgbox('PubMed could not be accessed.','Internet connection may be broken.','error');
    if exist(fname,'file'),
        delete(fname);
    end;
end;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global queries
global hMain

hMain.current_reference=[];

fname=strcat(general.tmp_files,queries.PubMed_file);
if exist(fname,'file'),
    delete(fname);
end;

delete(handles.my_figure);

function newstr=unblank(str)
% remove all blanks from an input string

newstr='';
for k=1:length(str),
    if ~strcmp(str(k),' '),
        newstr=[newstr str(k)];
    end;
end;


function handles=make_keyword_list(handles)

global model

found=0;
list={};
if isfield(model,'info'),
    if ~isempty(model.info),
        for k=1:length(model.info),
            if ~isempty(model.info{k}.keywords),
                keys=textscan(model.info{k}.keywords,'%s','Delimiter',',');
            else
                keys{1}='';
            end;
            if ~isempty(keys{1})>0,
                found=1;
                for kk=1:length(keys{1}),
                    list{kk}=strtrim(char(keys{1}(kk)));
                end;
                list=sort(list);
                set(handles.listbox_keywords,'String',list);
                handles.current_key=lower(list{1});
            end;
        end;
    end;
end;

if ~found,
    set(handles.listbox_keywords,'Enable','off');
    set(handles.pushbutton_plus,'Enable','off');
    set(handles.pushbutton_keyword_only,'Enable','off');
    handles.current_key='';
end;


% --- Executes on button press in checkbox_reviews.
function checkbox_reviews_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_reviews (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_reviews


% --- Executes on selection change in popupmenu_autosearch.
function popupmenu_autosearch_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_autosearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_autosearch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_autosearch

global model

sid=get(hObject,'Value');

num=length(model.autosearch);

if sid>num,
    set(handles.pushbutton_activate,'Enable','off');
    set(handles.pushbutton_deactivate,'Enable','off');
    set(handles.pushbutton_delete_auto,'Enable','off');
    return
end;

search=model.autosearch(sid);

set(handles.edit_search_terms,'String',search.search_terms);

switch search.mode
    case 1
        set(handles.radiobutton_all_terms,'Value',1);
    case 2
        set(handles.radiobutton_TIAB,'Value',1);
    case 3
        set(handles.radiobutton_keywords,'Value',1);
    case 4
        set(handles.radiobutton_major,'Value',1);
end;

set(handles.edit_authors,'String',search.authors);

switch search.author_mode
    case 1
        set(handles.radiobutton_all_authors,'Value',1);
    case 2
        set(handles.radiobutton_first_author,'Value',1);
    case 3
        set(handles.radiobutton_last_author,'Value',1);
end;

set(handles.checkbox_reviews,'Value',search.reviews);

set(handles.edit_years,'String',search.years);

set(handles.popupmenu_records,'Value',search.records);

set(handles.edit_journal,'String',search.journal);

if search.active,
    set(handles.pushbutton_activate,'Enable','off');
    set(handles.pushbutton_deactivate,'Enable','on');
else
    set(handles.pushbutton_activate,'Enable','on');
    set(handles.pushbutton_deactivate,'Enable','off');
end;
set(handles.pushbutton_delete_auto,'Enable','on');

guidata(hObject,handles);
 
% --- Executes during object creation, after setting all properties.
function popupmenu_autosearch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_autosearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_delete_auto.
function pushbutton_delete_auto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if isfield(model,'autosearch'),
    if ~isempty(model.autosearch),
        num=length(model.autosearch);
        sid=get(handles.popupmenu_autosearch,'Value');
        if sid<num,
            for k=sid+1:num,
                model.autosearch(k-1)=model.autosearch(k);
            end;
        end;
        if sid<=num,
            model.autosearch=model.autosearch(1:num-1);
        end;
    end;
end;

 handles = setup_autosearch(handles);
 
 guidata(hObject,handles);
 
% --- Executes on button press in pushbutton_autosearch.
function pushbutton_autosearch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_autosearch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

hMain.current_reference=[];

if isfield(model,'autosearch'),
    sid=length(model.autosearch)+1;
else
    sid=1;
end;

search.search_terms=get(handles.edit_search_terms,'String');

search.mode=1;
if get(handles.radiobutton_TIAB,'Value');
    search.mode=2;
end;
if get(handles.radiobutton_keywords,'Value');
    search.mode=3;
end;
if get(handles.radiobutton_major,'Value');
    search.mode=4;
end;

search.authors=get(handles.edit_authors,'String');

search.author_mode=1;
if get(handles.radiobutton_first_author,'Value');
    search.author_mode=2;
end;
if get(handles.radiobutton_last_author,'Value');
    search.author_mode=3;
end;

search.reviews=get(handles.checkbox_reviews,'Value');

search.years=get(handles.edit_years,'String');

search.records=get(handles.popupmenu_records,'Value');

search.journal=get(handles.edit_journal,'String');

name = char(inputdlg('Name:','Select name for this autosearch',1));

if isempty(name),
    msgbox('Autosearch cancelled','Empty name.','warn');
    return
end;

found=0;
if sid>1,
    for k=1:sid-1,
        if strcmp(model.autosearch(k).name,name),
            found=1;
            break
        end;
    end;
end;

if found,
    msgbox('Autosearch cancelled','Name already exists.','warn');
    return
end;

search.name=name;
search.active=1;
search.date=[];

model.autosearch(sid)=search;

handles=setup_autosearch(handles);

[handles,status]=do_search(handles);

if status>0,
    model.autosearch(sid).date=now;
    delete(handles.my_figure);
end;


% --- Executes on button press in pushbutton_activate.
function pushbutton_activate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_activate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

sid=get(handles.popupmenu_autosearch,'Value');

model.autosearch(sid).active=1;

set(hObject,'Enable','off');
set(handles.pushbutton_deactivate','Enable','on');

guidata(hObject,handles);
 
% --- Executes on button press in pushbutton_deactivate.
function pushbutton_deactivate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deactivate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

sid=get(handles.popupmenu_autosearch,'Value');

model.autosearch(sid).active=0;

set(hObject,'Enable','off');
set(handles.pushbutton_activate','Enable','on');

guidata(hObject,handles);


function handles = setup_autosearch(handles)

global model

if ~isfield(model,'autosearch')
    set(handles.popupmenu_autosearch,'Enable','off');
else
    if isempty(model.autosearch),
        set(handles.popupmenu_autosearch,'Enable','off');
    else
        num=length(model.autosearch);
        list={};
        for k=1:num,
            list{k}=model.autosearch(k).name;
        end;
        list{num+1}='<none>';
        set(handles.popupmenu_autosearch,'String',list);
        set(handles.popupmenu_autosearch,'Value',num+1);
    end;
end;
set(handles.pushbutton_activate,'Enable','off');
set(handles.pushbutton_deactivate,'Enable','off');
set(handles.pushbutton_delete_auto,'Enable','off');


% --- Executes on button press in checkbox_bibliography.
function checkbox_bibliography_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_bibliography (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_bibliography

if get(hObject,'Value')
    set(handles.my_figure,'Name','Search bibliography for reference');
    set(handles.text_terms,'String','Title words');
else
    set(handles.my_figure,'Name','Search PubmMed for literature');
    set(handles.text_terms,'String','Concepts');
end;

function search=get_search(handles)
% gathers information from the GUI into a search structure

search.search_terms=get(handles.edit_search_terms,'String');

search.mode=1;
if get(handles.radiobutton_TIAB,'Value');
    search.mode=2;
end;
if get(handles.radiobutton_keywords,'Value');
    search.mode=3;
end;
if get(handles.radiobutton_major,'Value');
    search.mode=4;
end;

search.authors=get(handles.edit_authors,'String');

search.author_mode=1;
if get(handles.radiobutton_first_author,'Value');
    search.author_mode=2;
end;
if get(handles.radiobutton_last_author,'Value');
    search.author_mode=3;
end;

search.reviews=get(handles.checkbox_reviews,'Value');

search.years=get(handles.edit_years,'String');

search.records=get(handles.popupmenu_records,'Value');

search.journal=get(handles.edit_journal,'String');

search.name='bibliography';
search.active=0;
search.date=[];

function ref=search_bibliography(handles)
% searches the existing bibliography for a reference
% instead of concepts, words in the title are searched (concepts mode
% radiobuttons are ignored)
% the number of records setting is ignored
%
% ref   array with all matching references, empty if none is matching

global model

ref=[];
search=get_search(handles);

% Check, whether the search is empty
empty=1;
if ~isempty(search.search_terms), empty=0; end;
if ~isempty(search.authors), empty=0; end;
if search.reviews, empty=0; end;
if ~isempty(search.years), empty=0; end;
if ~isempty(search.journal), empty=0; end;

if empty,
    add_msg_board('Empty search in bibliography was not performed.');
    return
end;

% preprocess search terms, journal, authors, and years
term={};
poi=0;
if ~isempty(search.search_terms),
    terms=textscan(search.search_terms,'%s','Delimiter',';');
    if ~isempty(terms{1}),
        for k=1:length(terms{1}),
            if ~isempty(strtrim(char(terms{1}(k)))),
                poi=poi+1;
                term{poi}=strtrim(upper(char(terms{1}(k))));
            end;
        end;
    end;
end;
% for journals, each word is checked separately to avoid incorrect
% mismatches as far as possible, findstr instead of strfind is used
journal={};
poi=0;
if ~isempty(search.journal),
    jterms=textscan(search.journal,'%s');
    if ~isempty(jterms{1}),
        for k=1:length(jterms{1}),
            if ~isempty(strtrim(char(jterms{1}(k)))),
                poi=poi+1;
                journal{poi}=strtrim(upper(char(jterms{1}(k))));
            end;
        end;
    end;
end;
author={};
poi=0;
if ~isempty(search.authors),
    authors=textscan(search.authors,'%s','Delimiter',';');
    if ~isempty(authors{1}),
        for k=1:length(authors{1}),
            if ~isempty(strtrim(char(authors{1}(k)))),
                poi=poi+1;
                author{poi}=strtrim(upper(char(authors{1}(k))));
            end;
        end;
    end;
end;
if ~isempty(author),
    if search.author_mode==2, % consider only first author
        old_author=author;
        clear author;
        author{1}=old_author{1};
    end;
    if search.author_mode==3, % consider only last author
        old_author=author;
        clear author
        author{1}=old_author{poi};
    end;
end;
syear=[];
eyear=[];
if ~isempty(search.years),
    [syear,eyear]=strtok(search.years,'-');
    syear=str2double(syear);
    if isnan(syear),
        add_msg_board('Year or first year is not a number. Ignoring years.');
        syear=[];
    end;
    if ~isempty(eyear),
        eyear=str2double(eyear(2:end));
        if isnan(syear),
            add_msg_board('Last year is not a number. Searching only for first year.');
            eyear=[];
        end;
    end;         
end;

if isfield(model,'references'),
    if ~isempty(model.references),
        for rec=1:length(model.references),
            matching=1; % by default the reference is matching
            % check reviews, first, as it is fastest
            if search.reviews,
                if model.references(rec).type~=9,
                    matching=0;
                end;
            end;
            % check years, third to last, as it requires processing
            if ~isempty(syear) && matching,
                if isempty(eyear), eyear=syear; end;
                ref_year=str2double(model.references(rec).year);
                if ref_year<syear || ref_year>eyear,
                    matching=0; break;
                end;
            end;
            % check title words, this is done later, as it
            % involves a loop
            if ~isempty(term) && matching, % don't do it, if the reference does not match anyway
                for k=1:length(term),
                    found=strfind(upper(model.references(rec).title),term{k});
                    if isempty(found), matching=0; break; end;
                end;
            end;
            % check journal, this is done very late, as it takes much time
            % (nested loops, one by one comparsion)
            if ~isempty(journal) && matching, % don't do it, if the reference does not match anyway
                match=zeros(1,length(journal));
                for k=1:length(journal),
                    % set up author journal of reference
                    ref_journal={};
                    poi=0;
                    if ~isempty(model.references(rec).journal),
                        jterms=textscan(model.references(rec).journal,'%s');
                        if ~isempty(jterms{1}),
                            for kk=1:length(jterms{1}),
                                if ~isempty(strtrim(char(jterms{1}(kk)))),
                                    poi=poi+1;
                                    ref_journal{poi}=strtrim(upper(char(jterms{1}(kk))));
                                end;
                            end;
                        end;
                    end;
                    if ~isempty(ref_journal),
                        for kk=1:length(ref_journal),
                            found=findstr(ref_journal{kk},journal{k}); % a rather lenient comparison by findstr
                            if ~isempty(found), match(k)=1; break; end;
                        end;
                    end;
                end;
                if sum(match==0), matching=0; end;
            end;
            % check authors, this is done last, as it takes most time
            % (nested loops, one by one comparsion)
            if ~isempty(author) && matching, % don't do it, if the reference does not match anyway
                for k=1:length(author),
                    match=0;
                    % set up author array of reference
                    ref_author={};
                    poi=0;
                    if ~isempty(model.references(rec).authors),
                        authors=textscan(model.references(rec).authors,'%s','Delimiter',';');
                        if ~isempty(authors{1}),
                            for kk=1:length(authors{1}),
                                if ~isempty(strtrim(char(authors{1}(kk)))),
                                    poi=poi+1;
                                    ref_author{poi}=strtrim(upper(char(authors{1}(kk))));
                                end;
                            end;
                        end;
                    end;
                    if ~isempty(ref_author),
                        if search.author_mode==2, % consider only first author
                            old_author=ref_author;
                            clear ref_author;
                            ref_author{1}=old_author{1};
                        end;
                        if search.author_mode==3, % consider only last author
                            old_author=ref_author;
                            clear ref_author;
                            ref_author{1}=old_author{poi};
                        end;
                        for kk=1:length(ref_author),
                            found=strfind(ref_author{kk},author{k});
                            if ~isempty(found), match=1; break; end;
                        end;
                        if ~match, matching=0; break; end;
                    end;
                end;
            end;
            if matching,
                ref=[ref rec];
            end;
        end;
    end;
end;
