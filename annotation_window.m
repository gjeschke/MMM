function varargout = annotation_window(varargin)
% ANNOTATION_WINDOW M-file for annotation_window.fig
%      ANNOTATION_WINDOW, by itself, creates a new ANNOTATION_WINDOW or raises the existing
%      singleton*.
%
%      H = ANNOTATION_WINDOW returns the handle to a new ANNOTATION_WINDOW or the handle to
%      the existing singleton*.
%
%      ANNOTATION_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANNOTATION_WINDOW.M with the given input arguments.
%
%      ANNOTATION_WINDOW('Property','Value',...) creates a new ANNOTATION_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before annotation_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to annotation_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help annotation_window

% Last Modified by GUIDE v2.5 12-Jan-2018 12:37:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @annotation_window_OpeningFcn, ...
                   'gui_OutputFcn',  @annotation_window_OutputFcn, ...
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


% --- Executes just before annotation_window is made visible.
function annotation_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to annotation_window (see VARARGIN)

% Choose default command line output for annotation_window

global geometry_settings
global MMM_icon
global model
global hMain
global annotation_keys

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.output = hObject;

indices=resolve_address('*');
if isempty(indices), add_msg_board('Nothing selected for annotation.'); delete(hObject); return; end;
indices=indices(indices>0);
[m,n]=size(indices);
if m>1, add_msg_board('More than one object selected. Annotation requires a single object.'); delete(hObject); return, end;

my_adr=mk_address(indices);
handles.oindices=indices;
set(hObject,'Name',sprintf('Annotation of object %s',my_adr));
[msg,info]=get_object(indices,'info');
n1=length(info);
[msg,mass]=get_object(indices,'mass','-nowater');
info{n1+1}=sprintf('Mass: %6.1f g/mol',mass);
[msg,xyz]=get_object(indices,'xyz');

info{n1+2}=sprintf('Mean coordinates: %4.1f, %4.1f, %4.1f Å',mean(xyz,1));

poi=n1+2;
if length(indices)==1,
    poi=poi+1;
    info{poi}='---------------------------------------------------';
    poi=poi+1;
    info{poi}='Information from PDB file:';
    poi0=poi;
    if isfield(model.info{indices},'title'),
        poi=poi+1;
        info{poi}=model.info{indices}.title;
    end;
    if isfield(model.info{indices},'class'),
        poi=poi+1;
        info{poi}=model.info{indices}.class;
    end;
    if isfield(model.info{indices},'authors'),
        poi=poi+1;
        info{poi}=model.info{indices}.authors;
    end;
    if isfield(model.info{indices},'depDate'),
        poi=poi+1;
        info{poi}=model.info{indices}.depDate;
    end;
    if isfield(model.info{indices},'keywords'),
        poi=poi+1;
        info{poi}='--- PDB keywords: ---';
        poi=poi+1;
        info{poi}=model.info{indices}.keywords;
    end;
    if isfield(model.info{indices},'atoms'),
        text=sprintf('%i atoms',model.info{indices}.atoms);
        if isfield(model.info{indices},'residues'),
            text=sprintf('%s and %i residues.',text,model.info{indices}.residues);
        end;
        poi=poi+1;
        info{poi}=text;
    end;
    if poi==poi0,
        poi=poi+1;
        info{poi}='### Warning ### No annotation information in this PDB file.';
    end;
end;

if n==4, % this is a residue and DSSP information should be displayed
    [message,dssp]=get_residue(indices,'dssp');
    if message.error==2,
        poi=poi+1;
        info{poi}=message.text;
    end;
    if ~isempty(dssp),
        secondary='';
        switch dssp.sec;
            case 'H'
                secondary='H (alpha helix)';
            case 'B'
                secondary='B (isolated beta bridge)';
            case 'E'
                secondary='E (extended strand)';
            case 'G'
                secondary='G (3/10 helix)';
            case 'I'
                secondary='I (pi helix)';
            case 'T'
                secondary='T (hydrogen bonded turn)';
            case 'S'
                secondary='S (bend)';
            otherwise
                secondary='not defined';
        end;
        poi=poi+1;
        info{poi}=' ';
        poi=poi+1;
        info{poi}=sprintf('DSSP secondary structure assignment: %s',secondary);
        poi=poi+1;
        info{poi}=sprintf('DSSP water accessibility: %i Å^2 (ignores cofactors)',dssp.acc);
        poi=poi+1;
        info{poi}=sprintf('Backbone dihedrals: phi = %6.1f°, psi = %6.1f°',dssp.phi,dssp.psi);
        poi=poi+1;
        info{poi}='see next page for hydrogen bonds of peptide carbon and amide nitrogen';
    end;
end;

handles.info=info;
handles.page=1;
set(handles.edit_annotation,'String',handles.info);

switch n
    case 1, % structure has no context
        radius=[];
        title='';
    case 2, % chain context
        radius=geometry_settings.chain_context;
        title='Chain context radius';
    case 3, % chain model context        
        radius=geometry_settings.chain_context;
        title='Chain context radius';
    case 4, % residue context        
        radius=geometry_settings.residue_context;
        title='Residue context radius';
    case 5, % atom context        
        radius=geometry_settings.atom_context;
        title='Atom context radius';
    case 6, % location context        
        radius=geometry_settings.atom_context;
        title='Location context radius';
end;

if ~isempty(radius),
    [cindices,info]=context(indices,radius);
else
    cindices=[];
    info='No context defined for a complete structure.';
end;

if n==4, % this is a residue and DSSP information should be displayed
    poi=numel(info)+1;
    info{poi}=' ';
    poi=poi+1;
    info{poi}='DSSP context information';
    if message.error==2,
        poi=poi+1;
        info{poi}=message.text;
    end;
    if ~isempty(dssp),
        full_dssp=model.info{indices(1)}.dssp;
        stag=mk_address(indices(1));
        if sum(dssp.bp)>0,
            poi=poi+1;
            info{poi}='Bridge partners for beta sheet';
            for k=1:2,
                if dssp.bp(k)>0 && dssp.bp(k)<=numel(full_dssp),
                    partner=sprintf('%s(%s)%s',stag,full_dssp(dssp.bp(k)).chain,full_dssp(dssp.bp(k)).tag);
                    poi=poi+1;
                    info{poi}=partner;
                end;
            end;
        end;
        poi=poi+1;
        info{poi}='Hydrogen bonding partners for amide nitrogen';
        if dssp.NHO(2).hp~=0,
            info{poi}=[info{poi} ' (bifurcated)'];
        end;
        for k=1:2,
            if dssp.NHO(k).hp~=0,
                num=dssp.my_num+dssp.NHO(k).hp;
                if num>0 && num<=numel(full_dssp),
                    partner=sprintf('%s(%s)%s',stag,full_dssp(num).chain,full_dssp(num).tag);
                    poi=poi+1;
                    info{poi}=sprintf('%s with energy %6.1f kcal/mol',partner,dssp.NHO(k).energy);
                end;
            end;
        end;
        poi=poi+1;
        info{poi}='Hydrogen bonding partners for carbonyl oxygen';
        if dssp.OHN(2).hp~=0,
            info{poi}=[info{poi} ' (bifurcated)'];
        end;
        for k=1:2,
            if dssp.OHN(k).hp~=0,
                num=dssp.my_num+dssp.OHN(k).hp;
                if num>0 && num<=numel(full_dssp),
                    partner=sprintf('%s(%s)%s',stag,full_dssp(num).chain,full_dssp(num).tag);
                    poi=poi+1;
                    info{poi}=sprintf('%s with energy %6.1f kcal/mol',partner,dssp.OHN(k).energy);
                end;
            end;
        end;
    end;
end;

handles.context=info;

[message,annotations]=get_annotations(indices);

if isempty(annotations),
    annotations.privacy=[];
    annotations.keywords=[];
    annotations.references=[];
    annotations.text{1}='';
end;

handles.current_annotations=annotations;

handles=update_keylist(handles);

handles=update_references(handles);

page=1;

if ~isempty(hMain.keyword_request),
    key=strtrim(hMain.keyword_request);
    types={};
    tpoi=0;
    for k=1:length(annotation_keys),
        if ~isempty(annotation_keys(k).keywords),
            for kk=1:length(annotation_keys(k).keywords),
                if strcmpi(key,strtrim(annotation_keys(k).keywords{kk})),
                    tpoi=tpoi+1;
                    types{tpoi}=annotation_keys(k).type;
                    break;
                end;
            end;
        end;
    end;
    if ~isempty(types),
        for k=1:length(types),
            test=strcat('*',types{k});
            for kk=1:length(annotations.text),
                old_text=char(annotations.text(kk));
                if ~isempty(old_text),
                    old_text=old_text(1,:);
                    token=strtok(old_text);
                    if strcmp(token,test),
                        page=kk+2;
                        break;
                    end;
                end;
            end;        
        end;
    end;
end;
hMain.keyword_request='';

handles.page=page;
handles=page_update(handles);

load helpicon
set(handles.pushbutton_help_annotation,'CData',cdata);

hMain.annotation_open=hObject;

hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes annotation_window wait for user response (see UIRESUME)
uiwait(handles.annotations);


% --- Outputs from this function are returned to the command line.
function varargout = annotation_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



% --- Executes on button press in pushbutton_previous_page.
function pushbutton_previous_page_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.page>1, 
    handles.page=handles.page-1;
    handles=page_update(handles);
end;

guidata(hObject,handles);

% --- Executes on button press in pushbutton_next_page.
function pushbutton_next_page_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.page=handles.page+1;
handles=page_update(handles);

guidata(hObject,handles);

function edit_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_annotation as text
%        str2double(get(hObject,'String')) returns contents of edit_annotation as a double

my_text=get(hObject,'String');
n=handles.page-2;
handles.current_annotations.text{n}=my_text;
handles.current_annotations.privacy(n)=0;
if get(handles.radiobutton_group,'Value'),
    handles.current_annotations.privacy(n)=1;
end;    
if get(handles.radiobutton_user,'Value'),
    handles.current_annotations.privacy(n)=2;
end;    
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_annotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_annotation (see GCBO)
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

% Hints: contents = get(hObject,'String') returns listbox_keywords contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_keywords


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

% --- Executes on button press in pushbutton_add_keyword.
function pushbutton_add_keyword_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_keyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

contents = get(handles.listbox_keywords,'String');
sel=get(handles.listbox_keywords,'Value');
keyword{1}=contents{sel};

register_keyword(keyword,handles.oindices);

if isfield(model,'keywords'),
   id=tag2id(keyword{1},model.keywords); % find keyword
   if ~isempty(id),
        n=length(handles.current_annotations.keywords);
        handles.current_annotations.keywords(n+1)=id;
   end;
end;

handles=page_update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_new_keyword.
function pushbutton_new_keyword_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new_keyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

key = inputdlg('Input new keyword','Keyword definition');
register_keyword(key,handles.oindices);
handles=update_keylist(handles);

if isfield(model,'keywords'),
   id=tag2id(key{1},model.keywords); % find keyword
   if ~isempty(id),
        n=length(handles.current_annotations.keywords);
        handles.current_annotations.keywords(n+1)=id;
   end;
end;

handles=page_update(handles);
guidata(hObject,handles);

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


% --- Executes on button press in pushbutton_add_reference.
function pushbutton_add_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

refnum = get(handles.listbox_references,'Value');
poi= get(handles.listbox_references,'UserData');
refnum=poi(refnum);

n=length(handles.current_annotations.references);
handles.current_annotations.references(n+1)=refnum;

handles=page_update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_new_reference.
function pushbutton_new_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.current_reference=0;
reference_window;
if hMain.current_reference>0,
    n=length(handles.current_annotations.references);
    handles.current_annotations.references(n+1)=hMain.current_reference;
end;
handles=update_references(handles);
handles=page_update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_okay.
function pushbutton_okay_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_okay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

annotations=handles.current_annotations;
n=handles.page-2;
if n>1,
    if get(handles.radiobutton_public,'Value'),
        annotations.privacy(n)=0;
    end;
    if get(handles.radiobutton_group,'Value'),
        annotations.privacy(n)=1;
    end;
    if get(handles.radiobutton_user,'Value'),
        annotations.privacy(n)=2;
    end;

    annotations.text{n}=get(handles.edit_annotation,'String');
end;
message=set_annotations(handles.oindices,annotations);
annotations_CloseRequestFcn(handles.annotations, eventdata, handles);
% delete(handles.annotations);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


button = questdlg('All annotation changes will be lost','Are you sure?','No','Yes','No');

if strcmp(button,'Yes'),
    annotations_CloseRequestFcn(handles.annotations, eventdata, handles);
    % delete(handles.annotations);
end;

% --- Executes on button press in pushbutton_remove_keyword.
function pushbutton_remove_keyword_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_keyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

contents = get(handles.listbox_keywords,'String');
sel=get(handles.listbox_keywords,'Value');
keyword{1}=contents{sel};

unregister_keyword(keyword,handles.oindices);

if isfield(model,'keywords'),
   id=tag2id(keyword{1},model.keywords); % find keyword
   keys=[];
   for k=1:length(handles.current_annotations.keywords),
       if handles.current_annotations.keywords(k)~=id, 
           keys=[keys handles.current_annotations.keywords(k)];
       end;
   end;
   handles.current_annotations.keywords=keys;
end;

handles=page_update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_remove_reference.
function pushbutton_remove_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

refnum = get(handles.listbox_references,'Value');
poi=get(handles.listbox_references,'UserData');
refnum=poi(refnum);

n=length(handles.current_annotations.references);
refs=[];
for k=1:n,
   if handles.current_annotations.references(k)~=refnum, 
       refs=[refs handles.current_annotations.references(k)];
   end;
end;
handles.current_annotations.references=refs;

handles=page_update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_edit_keyword.
function pushbutton_edit_keyword_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_edit_keyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

contents = get(handles.listbox_keywords,'String');
sel=get(handles.listbox_keywords,'Value');
keyword=contents{sel};


if isfield(model,'keywords'),
    poi=findstr(keyword,model.keywords);
    head=model.keywords(1:poi-1);
    tail=model.keywords(poi+length(keyword):end);
    key = inputdlg('Edit keyword','Keyword definition',1,{keyword});
    model.keywords=[head key{1} tail];
else
    msgbox('Please use New button to define a keyword','No keywords defined yet','error');
end;

handles=update_keylist(handles);
handles=page_update(handles);


% --- Executes on button press in pushbutton_edit_reference.
function pushbutton_edit_reference_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_edit_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

sel=get(handles.listbox_references,'Value');
poi=get(handles.listbox_references,'UserData');
hMain.current_reference = poi(sel);
if ~isfield(model,'references') || isempty(model.references),
    hMain.current_reference=0;
end;
set(handles.listbox_references,'Value',1);
hMain.annotation_references=handles.current_annotations.references;
reference_window;
handles.current_annotations.references=hMain.annotation_references;
handles=update_references(handles);
handles=page_update(handles);
guidata(hObject,handles);

function handles=page_update(handles)

global model

set(handles.text_page,'String',sprintf('%i',handles.page));

annotations=handles.current_annotations;

% keyboard

if ~isempty(annotations.keywords),
    keystring='';
    for k=1:length(annotations.keywords),
        keystring=sprintf('%s; %s',keystring,id2tag(annotations.keywords(k),model.keywords));
    end;
    keystring=keystring(2:end);
    set(handles.text_keywords,'String',keystring);
else
    set(handles.text_keywords,'String','');
end;
if ~isempty(annotations.references),
    refstring='';
    for k=1:length(annotations.references),
        new_ref=model.references(annotations.references(k)).short;
        refstring=sprintf('%s; %s',refstring,new_ref);
    end;
    refstring=refstring(2:end);
    set(handles.text_references,'String',refstring);
else
    set(handles.text_references,'String','');
end;

if handles.page==1,
    set(handles.edit_annotation,'String',handles.info);
    set(handles.edit_annotation,'Enable','inactive');
    return
elseif handles.page==2,
    set(handles.edit_annotation,'String',handles.context);
    set(handles.edit_annotation,'Enable','inactive');
    return
else
    set(handles.edit_annotation,'Enable','on');
    annotations=handles.current_annotations;
    existing=0;
    number=handles.page-2;
    if ~isempty(annotations),
        n=length(annotations.privacy);
        if n>=number,
            existing=1;
            switch annotations.privacy(number)
                case 0
                    set(handles.radiobutton_public,'Value',1);
                case 1
                    set(handles.radiobutton_group,'Value',1);
                case 2
                    set(handles.radiobutton_user,'Value',1);
            end;
            set(handles.edit_annotation,'String',annotations.text{number});
        end;
        if ~existing, % new annotation
            set(handles.radiobutton_public,'Value',1);
            set(handles.edit_annotation,'String','');
        end;
    end;
end;


function register_keyword(key,indices)

global model

if isempty(key),
    return;
end;

cindices=zeros(1,6);
cindices(1:length(indices))=indices;

if isfield(model,'keywords'),
   for k=1:length(key),
       id=tag2id(key{k},model.keywords); % find keyword
       if isempty(id), % keyword is new
           model.keywords=sprintf('%s%s:',model.keywords,key{k}); % extend keyword list
           poi=length(model.keys);
           model.keys(poi+1).indices=cindices;
       else
           found=0;
           [m,n]=size(model.keys(id).indices);
           for kk=1:m,
               if sum(abs(model.keys(id).indices(kk,:)-cindices))==0,
                   found=1;
                   break;
               end;
           end;
           if ~found,
               model.keys(id).indices=[model.keys(id).indices; cindices];
           end;
       end;
   end;
else
    model.keywords=sprintf(':%s:',key{1});
    model.keys(1).indices=cindices;
    if length(key)>1,
        for k=2:length(key),
           model.keywords=sprintf('%s%s:',model.keywords,key{k});
           model.keys(k).indices=cindices;
        end;
    end;
end;

function unregister_keyword(key,indices)

global model

cindices=zeros(1,6);
cindices(1:length(indices))=indices;

if isfield(model,'keywords'),
   for k=1:length(key),
       id=tag2id(key{k},model.keywords); % find keyword
       if ~isempty(id), % keyword is new
           oldindices=model.keys(id).indices;
           newindices=oldindices;
           [m,n]=size(oldindices);
           [m1,n1]=size(cindices);
           poi=0;
           for k=1:m,
               found=0;
               for k1=1:m1,
                   if sum(oldindices(k,:)-cindices(k1,:))==0,
                       found=1;
                       break
                   end;
               end;
               if ~found,
                   poi=poi+1;
                   newindices(poi,:)=oldindices(k,:);
               end;
           end;
           model.keys(id).indices=newindices(1:poi,:);
       end;
   end;
else
    for k=1:length(key),
       model.keywords=sprintf(':%s:',key{k});
       model.keys(k).indices=cindices;
    end;
end;

function handles=update_keylist(handles)

global model

if isfield(model,'keywords'),
    for k=1:length(model.keys),
        keylist{k}=id2tag(k,model.keywords);
    end;
else
    keylist{1}='';
end;
[keylist,poi]=sort(keylist);
set(handles.listbox_keywords,'String',keylist);
set(handles.listbox_keywords,'UserData',poi);

function handles=update_references(handles)

global model

reflist={''};
if isfield(model,'references'),
    for k=1:length(model.references),
        reflist{k}=model.references(k).short;
    end;
end;
[reflist,poi]=sort(reflist);
set(handles.listbox_references,'String',reflist);
set(handles.listbox_references,'UserData',poi);


% --- Executes on button press in radiobutton_public.
function radiobutton_public_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_public (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_public

set(handles.radiobutton_user,'Value',0);
set(handles.radiobutton_group,'Value',0);

page=handles.page-2;
if page>0,
    handles.current_annotations.privacy(page)=0;
end;

guidata(hObject,handles);

% --- Executes on button press in radiobutton_group.
function radiobutton_group_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_group

group=get(hObject,'Value');

page=handles.page-2;
if group,
    set(handles.radiobutton_public,'Value',0);
    set(handles.radiobutton_user,'Value',0);
    if page>0,
        handles.current_annotations.privacy(page)=1;
    end;
else
    set(handles.radiobutton_public,'Value',1);
    set(handles.radiobutton_user,'Value',0);
    if page>0,
        handles.current_annotations.privacy(page)=0;
    end;
end;

guidata(hObject,handles);


% --- Executes on button press in radiobutton_user.
function radiobutton_user_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_user

user=get(hObject,'Value');

page=handles.page-2;
if user,
    set(handles.radiobutton_public,'Value',0);
    set(handles.radiobutton_group,'Value',0);
    if page>0,
        handles.current_annotations.privacy(page)=2;
    end;
else
    set(handles.radiobutton_public,'Value',1);
    set(handles.radiobutton_user,'Value',0);
    if page>0,
        handles.current_annotations.privacy(page)=0;
    end;
end;

guidata(hObject,handles);


% --- Executes on button press in pushbutton_help_annotation.
function pushbutton_help_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'annotation_window.html');
webcall(entry,'-helpbrowser');


% --- Executes when user attempts to close annotations.
function annotations_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to annotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global hMain

delete(hMain.annotation_open);
hMain.annotation_open=[];
