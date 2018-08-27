function varargout = graphics_selection(varargin)
% GRAPHICS_SELECTION M-file for graphics_selection.fig
%      GRAPHICS_SELECTION, by itself, creates a new GRAPHICS_SELECTION or raises the existing
%      singleton*.
%
%      H = GRAPHICS_SELECTION returns the handle to a new GRAPHICS_SELECTION or the handle to
%      the existing singleton*.
%
%      GRAPHICS_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRAPHICS_SELECTION.M with the given input arguments.
%
%      GRAPHICS_SELECTION('Property','Value',...) creates a new GRAPHICS_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before graphics_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to graphics_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help graphics_selection

% Last Modified by GUIDE v2.5 12-Jan-2018 12:47:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @graphics_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @graphics_selection_OutputFcn, ...
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


% --- Executes just before graphics_selection is made visible.
function graphics_selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to graphics_selection (see VARARGIN)

global MMM_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

% Choose default command line output for graphics_selection
handles.output = hObject;

load helpicon
set(handles.pushbutton_help_graphics,'CData',cdata);

handles=setup_surfaces(handles);

handles=setup_motion_arrows(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes graphics_selection wait for user response (see UIRESUME)
% uiwait(handles.figure);


% --- Outputs from this function are returned to the command line.
function varargout = graphics_selection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_okay.
function pushbutton_okay_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_okay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

set(handles.figure,'Pointer','watch');

axes(hMain.axes_model);

atom_mode=[];
if get(handles.radiobutton_atoms_wire,'Value'),
    atom_mode='wire';
end;
if get(handles.radiobutton_atoms_stick,'Value'),
    atom_mode='stick';
end;
if get(handles.radiobutton_atoms_ballstick,'Value'),
    atom_mode='ball&stick';
end;
if get(handles.radiobutton_atoms_spacefilling,'Value'),
    atom_mode='space-filling';
end;

residue_mode=[];
if get(handles.radiobutton_residues_wire,'Value'),
    residue_mode='CaWire';
end;
if get(handles.radiobutton_residues_coil,'Value'),
    residue_mode='coil';
end;
if get(handles.radiobutton_residues_ribbon,'Value'),
    residue_mode='ribbon';
end;
if get(handles.radiobutton_residues_allatom,'Value'),
    residue_mode=atom_mode;
end;

coarse_mode=[];
if get(handles.radiobutton_coarse_cartoon,'Value'),
    coarse_mode='cartoon';
end;
if get(handles.radiobutton_coarse_density,'Value'),
    coarse_mode='density';
end;
if get(handles.radiobutton_coarse_coil,'Value'),
    coarse_mode='string';
end;

water=get(handles.checkbox_water,'Value');

labels=get(handles.checkbox_labels,'Value');

if isfield(model,'selections') && isfield(model,'selected')
    if ~isempty(model.selected),
        for k=1:length(model.selected),
            indices=model.selected{k};
            [m,n]=size(indices);
            for kk=1:m,
                currind=indices(kk,:);
                currind=currind(currind>0);
                switch length(currind)
                    case {1,2,3} % to be extended for coarse mode
                        if isempty(residue_mode) || isempty(coarse_mode),
                            set_object(currind,'hide');
                        end;
                        if ~isempty(coarse_mode),
                            set_object(currind,'show',{coarse_mode});
                        end;
                        if ~isempty(residue_mode)
                            set_object(currind,'show',{residue_mode});
                        end;
                        if water
                            set_object(currind,'show',{'water'});
                        end;
                        if labels
                            set_object(currind,'show',{'label'});
                        end;
                    case 4 
                        if isempty(residue_mode)
                            resname=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
                            reswater=strcmpi(resname,'HOH'); % check whether this is a water
                            if ~reswater
                                set_object(currind,'hide');
                            elseif water
                                set_object(currind,'show',{'water'});
                            else
                                set_object(currind,'hide');
                            end;
                        else
                            fres_mode=residue_mode;
                            if water, fres_mode=strcat(fres_mode,'!'); end;
                            set_object(currind,'show',{fres_mode});
                        end;
                        set_object(currind,'show',{'label_hide'});
                        if labels,
                            set_object(currind,'show',{'label'});
                        end;
                    case {5,6}
                        if isempty(atom_mode)
                            set_object(currind,'hide');
                        else
                            set_object(currind,'show',{atom_mode});
                        end;
                end;
            end;
        end;
    end;
end;

% now treat the surfaces
if isfield(model,'surfaces') && ~isempty(model.surfaces),
    active=get(handles.popupmenu_surfaces,'UserData');
    for k=1:length(model.surfaces),
        objects=model.surfaces(k).gobjects;
        if active(k), state='on'; else state='off'; end;
        if ~isempty(objects),
            for kk=1:length(objects),
                if ishandle(objects(kk)),
                    set(objects(kk),'Visible',state);
                end;
            end;
        end;
    end;
end;

% and the sets of motion arrows
if isfield(model,'motion') && ~isempty(model.motion),
    active=get(handles.popupmenu_motion,'UserData');
    for k=1:length(model.motion),
        objects=model.motion(k).gobjects;
        if active(k), state='on'; else state='off'; end;
        if ~isempty(objects),
            for kk=1:length(objects),
                if ishandle(objects(kk)),
                    set(objects(kk),'Visible',state);
                end;
            end;
        end;
    end;
end;

lighting gouraud
highlight_selection;
unrecord_objects;
set(handles.figure,'Pointer','arrow');
delete(handles.figure);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure);

% --- Executes on button press in checkbox_water.
function checkbox_water_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_water


% --- Executes on button press in checkbox_labels.
function checkbox_labels_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_labels



% --- Executes on button press in checkbox_labels_light.
function checkbox_labels_light_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_labels_light (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_labels_light

val=get(hObject,'Value');
if val,
    set(handles.checkbox_labels,'Value',0);
end;


% --- Executes on button press in pushbutton_help_graphics.
function pushbutton_help_graphics_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_graphics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'graphics_window.html');
webcall(entry,'-helpbrowser');


% --- Executes on selection change in popupmenu_surfaces.
function popupmenu_surfaces_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_surfaces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_surfaces contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_surfaces

myval=get(hObject,'Value');
active=get(hObject,'UserData');
set(handles.checkbox_surface_active,'Value',active(myval));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_surfaces_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_surfaces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_surface_active.
function checkbox_surface_active_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_surface_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_surface_active

global model

myval=get(hObject,'Value');
ind=get(handles.popupmenu_surfaces,'Value');
model.surfaces(ind).active=myval;
active=get(handles.popupmenu_surfaces,'UserData');
active(ind)=myval;
set(handles.popupmenu_surfaces,'UserData',active);

guidata(hObject,handles);

function handles=setup_surfaces(handles)

global model

if isfield(model,'surfaces') && ~isempty(model.surfaces),
    active=zeros(1,length(model.surfaces));
    for k=1:length(model.surfaces),
        list{k}=model.surfaces(k).tag;
        active(k)=model.surfaces(k).active;
    end;
    set(handles.popupmenu_surfaces,'String',list);
    set(handles.popupmenu_surfaces,'Value',1);
    set(handles.popupmenu_surfaces,'UserData',active);
    set(handles.checkbox_surface_active,'Value',active(1));
else
    set(handles.popupmenu_surfaces,'Enable','off');
    set(handles.checkbox_surface_active,'Enable','off');
end;

function handles=setup_motion_arrows(handles)

global model

if isfield(model,'motion') && ~isempty(model.motion),
    active=zeros(1,length(model.motion));
    for k=1:length(model.motion),
        list{k}=model.motion(k).tag;
        active(k)=model.motion(k).active;
    end;
    set(handles.popupmenu_motion,'String',list);
    set(handles.popupmenu_motion,'Value',1);
    set(handles.popupmenu_motion,'UserData',active);
    set(handles.checkbox_motion_active,'Value',active(1));
else
    set(handles.popupmenu_motion,'Enable','off');
    set(handles.checkbox_motion_active,'Enable','off');
end;


% --- Executes on selection change in popupmenu_motion.
function popupmenu_motion_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_motion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_motion

myval=get(hObject,'Value');
active=get(hObject,'UserData');
set(handles.checkbox_motion_active,'Value',active(myval));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_motion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_motion_active.
function checkbox_motion_active_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_motion_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_motion_active

global model

myval=get(hObject,'Value');
ind=get(handles.popupmenu_motion,'Value');
model.motion(ind).active=myval;
active=get(handles.popupmenu_motion,'UserData');
active(ind)=myval;
set(handles.popupmenu_motion,'UserData',active);

guidata(hObject,handles);
