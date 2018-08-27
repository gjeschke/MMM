function varargout = depth_cue_control(varargin)
% DEPTH_CUE_CONTROL M-file for depth_cue_control.fig
%      DEPTH_CUE_CONTROL, by itself, creates a new DEPTH_CUE_CONTROL or raises the existing
%      singleton*.
%
%      H = DEPTH_CUE_CONTROL returns the handle to a new DEPTH_CUE_CONTROL or the handle to
%      the existing singleton*.
%
%      DEPTH_CUE_CONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEPTH_CUE_CONTROL.M with the given input arguments.
%
%      DEPTH_CUE_CONTROL('Property','Value',...) creates a new DEPTH_CUE_CONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before depth_cue_control_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to depth_cue_control_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help depth_cue_control

% Last Modified by GUIDE v2.5 12-Jan-2018 12:41:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @depth_cue_control_OpeningFcn, ...
                   'gui_OutputFcn',  @depth_cue_control_OutputFcn, ...
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


% --- Executes just before depth_cue_control is made visible.
function depth_cue_control_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to depth_cue_control (see VARARGIN)

% Choose default command line output for depth_cue_control
handles.output = hObject;

global hMain
global model
global MMM_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.old_bckg_color=hMain.color;
hMain.depthcue_control=hObject;

command=sprintf('bckg white');
cmd(hMain,command);

if isfield(model,'selections')
    if model.selections>=1,
        set(handles.slider_focus,'Value',0);
        set(handles.slider_back,'Value',1.5);
        set(handles.edit_focus,'String','0.00');
        set(handles.edit_back,'String','1.50');
    else
        set(handles.slider_focus,'Value',0);
        set(handles.slider_back,'Value',1.5);
        set(handles.edit_focus,'String','0.00');
        set(handles.edit_back,'String','1.50');
    end;
else
    set(handles.slider_focus,'Value',0);
    set(handles.slider_back,'Value',1.5);
    set(handles.edit_focus,'String','0.00');
    set(handles.edit_back,'String','1.50');
end;
guidata(hObject,handles);
handles=depth_cueing(handles);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes depth_cue_control wait for user response (see UIRESUME)
% uiwait(handles.figure);


% --- Outputs from this function are returned to the command line.
function varargout = depth_cue_control_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_focus_Callback(hObject, eventdata, handles)
% hObject    handle to slider_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

newval=get(hObject,'Value');
set(handles.edit_focus,'String',sprintf('%6.2f',newval));

% --- Executes during object creation, after setting all properties.
function slider_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_focus_Callback(hObject, eventdata, handles)
% hObject    handle to edit_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_focus as text
%        str2double(get(hObject,'String')) returns contents of edit_focus as a double

oldval=get(handles.slider_focus,'Value');
minval=get(handles.slider_focus,'Min');
maxval=get(handles.slider_focus,'Max');

newval=str2double(get(hObject,'String'));
if isnan(newval),
    set(hObject,'String',sprintf('%6.2f',oldval));
    msgbox('Only numbers are allowed input. Focus plane position was reset to previous value.','Input error','warn');
else
    if newval<minval,
        set(handles.slider_focus,'Min',newval);
    end;
    if newval>maxval,
        set(handles.slider_focus,'Max',newval);
    end;
    set(handles.slider_focus,'Value',newval);
    set(hObject,'String',sprintf('%6.2f',newval));
end;

% --- Executes during object creation, after setting all properties.
function edit_focus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_back_Callback(hObject, eventdata, handles)
% hObject    handle to slider_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

newval=get(hObject,'Value');
set(handles.edit_back,'String',sprintf('%6.2f',newval));


% --- Executes during object creation, after setting all properties.
function slider_back_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_back_Callback(hObject, eventdata, handles)
% hObject    handle to edit_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_back as text
%        str2double(get(hObject,'String')) returns contents of edit_back as a double

oldval=get(handles.slider_back,'Value');
minval=get(handles.slider_back,'Min');
maxval=get(handles.slider_back,'Max');

newval=str2double(get(hObject,'String'));
if isnan(newval),
    set(hObject,'String',sprintf('%6.2f',oldval));
    msgbox('Only numbers are allowed input. Back plane position was reset to previous value.','Input error','warn');
else
    if newval<minval,
        set(handles.slider_back,'Min',newval);
    end;
    if newval>maxval,
        set(handles.slider_back,'Max',newval);
    end;
    set(handles.slider_back,'Value',newval);
    set(hObject,'String',sprintf('%6.2f',newval));
end;



% --- Executes during object creation, after setting all properties.
function edit_back_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=depth_cueing(handles);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

set(hMain.togglebutton_depth_cueing,'Value',0);
handles=cancel_depth_cueing(handles);
delete(handles.figure);

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

set(hMain.togglebutton_depth_cueing,'Value',1);
delete(handles.figure);

function handles=depth_cueing(handles)

global hMain
global model
global hModel

set(hMain.figure,'Pointer','watch');
if hMain.detached,
    set(hModel.figure,'Pointer','watch');
end;

switch hMain.color,
    case 'white'
        bckg=[1,1,1];
    case 'grey'
        bckg=[0.941,0.941,0.941];
    case 'black'
        bckg=[0,0,0];
end;


m=model.graphics_lookup_pointer;
campos=get(hMain.axes_model,'CameraPosition');
camtarget=get(hMain.axes_model,'CameraTarget');
camvec=camtarget-campos;
camvec=camvec/norm(camvec);
z=linspace(0,1,m);
xyzmat=zeros(m,3);
for k=1:m,
    xyzext=model.graphics_xyz(k,:);
    xyz=(xyzext(1:3)+xyzext(4:6))/2-campos;
    xyzmat(k,:)=xyz;
    z(k)=norm(xyz.*camvec);
end;

minz=min(z);
maxz=max(z);

selection_flag=0;
if isfield(model,'selections')
    if model.selections>=1,
        minz=1e6;
        maxz=-1e6;
        selection_flag=1;
        z2=zeros(1,model.selections);
        for k=1:model.selections,
            cindices=model.selected{k};
            cindices=cindices(cindices>0);
            [msg,xyz]=get_object(cindices,'xyz');
            ma=max(xyz);
            ma=norm((ma-campos).*camvec);
            mi=min(xyz);
            mi=norm((mi-campos).*camvec);
            if mi>ma,
                tmp=mi; mi=ma; ma=tmp;
            end;
            if ma>maxz, maxz=ma; end;
            if mi<minz, minz=mi; end;
            z2(k)=norm((mean(xyz)-campos).*camvec);
        end;
        maxz=maxz+2.5;
        minz=minz-2.5;
    end;
end;

if minz>=maxz,
    minz=maxz-2.5;
    maxz=minz+2.5;
end;

dz=(z-minz)/(maxz-minz);

backplane=get(handles.slider_back,'Value');
focusplane=get(handles.slider_focus,'Value');
dz=(dz-focusplane)/(backplane-focusplane);
trans=ones(size(dz));
trans(find(dz<focusplane))=0;
dz(find(dz<0))=0;
dz(find(dz>1))=1;
dz=sqrt(dz);
db=1-dz;
dz=dz.^2;
admix=0.3;
ambient=sqrt((1-admix)*dz+admix);
for kk=1:m,
    indices=model.graphics_lookup(kk,2:end);
    if indices(4)<0,
        [msg,allgraphics]=get_object(indices(1:3),'secondary_graphics');
    else
        indices=indices(indices>0);
        [msg,allgraphics]=get_object(indices,'graphics');
    end;
    if ~isempty(allgraphics),
        for kkk=1:length(allgraphics),
            graphics=allgraphics(kkk);
            if ~isempty(graphics),
                setcolor=dz(kk)*bckg+db(kk)*graphics.color(1,:);
                for k=1:length(graphics.objects),
                    if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                        if isprop(graphics.objects(k),'Color'),
                            set(graphics.objects(k),'Color',setcolor);
                        elseif isprop(graphics.objects(k),'FaceColor'),
                            set(graphics.objects(k),'FaceColor',setcolor,'AmbientStrength',ambient(kk),'FaceAlpha',graphics.opaque(1,:)*trans(kk));
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

set(hMain.figure,'Pointer','arrow');
if hMain.detached,
    set(hModel.figure,'Pointer','arrow');
end;

function handles=cancel_depth_cueing(handles)

global hMain
global model
global hModel

set(hMain.figure,'Pointer','watch');
if hMain.detached,
    set(hModel.figure,'Pointer','watch');
end;

m=model.graphics_lookup_pointer;

for kk=1:m,
    indices=model.graphics_lookup(kk,2:end);
    if indices(4)<0,
        [msg,allgraphics]=get_object(indices(1:3),'secondary_graphics');
    else
        indices=indices(indices>0);
        [msg,allgraphics]=get_object(indices,'graphics');
    end;
    if ~isempty(allgraphics),
        for kkk=1:length(allgraphics),
            graphics=allgraphics(kkk);
            if ~isempty(graphics),
                for k=1:length(graphics.objects),
                    if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                        if isprop(graphics.objects(k),'Color'),
                            set(graphics.objects(k),'Color',graphics.color(1,:));
                        elseif isprop(graphics.objects(k),'FaceColor'),
                            set(graphics.objects(k),'FaceColor',graphics.color(1,:),'FaceAlpha',graphics.opaque(1,:),'AmbientStrength',0.3);
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

highlight_selection;
camlight(hMain.camlight);

command=sprintf('bckg %s',handles.old_bckg_color);
cmd(hMain,command);

set(hMain.figure,'Pointer','arrow');
if hMain.detached,
    set(hModel.figure,'Pointer','arrow');
end;



% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global hMain

set(hMain.togglebutton_depth_cueing,'Value',0);
handles=cancel_depth_cueing(handles);

delete(hObject);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'depth_cueing.html');
webcall(entry,'-helpbrowser');
