function varargout = display_density(varargin)
% DISPLAY_DENSITY M-file for display_density.fig
%      DISPLAY_DENSITY, by itself, creates a new DISPLAY_DENSITY or raises the existing
%      singleton*.
%
%      H = DISPLAY_DENSITY returns the handle to a new DISPLAY_DENSITY or the handle to
%      the existing singleton*.
%
%      DISPLAY_DENSITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAY_DENSITY.M with the given input arguments.
%
%      DISPLAY_DENSITY('Property','Value',...) creates a new DISPLAY_DENSITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before display_density_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to display_density_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help display_density

% Last Modified by GUIDE v2.5 12-Jan-2018 12:42:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @display_density_OpeningFcn, ...
                   'gui_OutputFcn',  @display_density_OutputFcn, ...
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


% --- Executes just before display_density is made visible.
function display_density_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to display_density (see VARARGIN)

global model

% Choose default command line output for display_density
handles.output = hObject;

if ~isfield(model,'density_tags') || ~isfield(model,'densities') || isempty(model.densities),
    add_msg_board('ERROR: No density cubes defined. Use "File/Load density" first.');
    delete(handles.figure1);
    return;
else
    n=length(model.densities);
    for k=1:n,
        dlist{k}=model.densities{k}.tag;
    end;
    set(handles.popupmenu_tag,'String',dlist);
end;

handles.rgb=[0,0,1];
handles.falpha=0.5;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes display_density wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = display_density_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



% --- Executes on selection change in popupmenu_tag.
function popupmenu_tag_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_tag


% --- Executes during object creation, after setting all properties.
function popupmenu_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_color.
function pushbutton_color_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rgb=uisetcolor;
guidata(hObject,handles);

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

id=get(handles.popupmenu_tag,'Value');
contents=get(handles.popupmenu_tag,'String');
tag0=contents{id};
tag=['density:' strtrim(tag0)];

new_graph=[];
if isfield(model,'surfaces'),
    for k=1:length(model.surfaces),
        if strcmp(model.surfaces(k).tag,tag),
            if ishandle(model.surfaces(k).gobjects),
                delete(model.surfaces(k).gobjects);
            end;
        end;
    end;
end;

level=get(handles.slider_level,'Value');

if level>=1-1e-3,
    delete(handles.figure1);
    return
end;

cube=model.densities{id}.cube;
cube=cube/max(max(max(cube)));
x=model.densities{id}.x;
y=model.densities{id}.y;
z=model.densities{id}.z;

set(handles.figure1,'Pointer','watch');
drawnow;

if get(handles.checkbox_normalize,'Value'),
    sdens=sum(sum(sum(cube)));
    level0=level;
    level=0.9985;
    for k=1:99,
        mask=(cube>=k/100);
        test=sum(sum(sum(mask.*cube)));
        if test<=level0*sdens,
            level=k/100;
            break;
        end;
    end;
end;

cube_tag = handles.popupmenu_tag.String{handles.popupmenu_tag.Value};
poi = strfind(cube_tag,':');
if ~isempty(poi)
    cube_type = cube_tag(1:poi-1);
else
    cube_type = 'unknown';
end

axes(hMain.axes_model);
switch cube_type
    case 'MMMx'
        [xg,yg,zg]=meshgrid(y,x,z);
        p = patch(isosurface(yg,xg,zg,cube,level));
    otherwise
        [xg,yg,zg]=meshgrid(x,y,z);
        p = patch(isosurface(xg,yg,zg,cube,level));
end
    
set(p, 'FaceColor', handles.rgb, 'EdgeColor', 'none','FaceAlpha',handles.falpha,'FaceLighting','gouraud','Clipping','off');
set(p, 'CDataMapping','direct','AlphaDataMapping','none');
set(handles.figure1,'Pointer','arrow');
dg.gobjects=p;
dg.tag=tag;
dg.color=handles.rgb;
dg.level=level;
dg.transparency=handles.falpha;
dg.active=true;

if isfield(model,'surfaces'),
    template = model.surfaces(1);
    names = fieldnames(template);
    for k = 1:length(names),
        if ~isfield(dg,names{k}),
            dg.(names{k}) = [];
        end;
    end;
    model.surfaces=[model.surfaces dg];
else
    model.surfaces=dg;
end;
camlookat(hMain.axes_model);
delete(handles.figure1);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

add_msg_board('Display of density cancelled by user.');
delete(handles.figure1);

% --- Executes on button press in checkbox_normalize.
function checkbox_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_normalize



function edit_level_Callback(hObject, eventdata, handles)
% hObject    handle to edit_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_level as text
%        str2double(get(hObject,'String')) returns contents of edit_level as a double

[v,handles]=edit_update_MMM(handles,hObject,0,1,0.5,'%4.2f',0);
set(handles.slider_level,'Value',v);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_level_Callback(hObject, eventdata, handles)
% hObject    handle to slider_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=get(handles.slider_level,'Value');
set(handles.edit_level,'String',sprintf('%4.2f',val));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_transparency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_transparency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_transparency as text
%        str2double(get(hObject,'String')) returns contents of edit_transparency as a double

[v,handles]=edit_update_MMM(handles,hObject,0,1,1,'%4.2f',0);
handles.falpha=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_transparency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_transparency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
