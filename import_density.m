function varargout = import_density(varargin)
% IMPORT_DENSITY M-file for import_density.fig
%      IMPORT_DENSITY, by itself, creates a new IMPORT_DENSITY or raises the existing
%      singleton*.
%
%      H = IMPORT_DENSITY returns the handle to a new IMPORT_DENSITY or the handle to
%      the existing singleton*.
%
%      IMPORT_DENSITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPORT_DENSITY.M with the given input arguments.
%
%      IMPORT_DENSITY('Property','Value',...) creates a new IMPORT_DENSITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before import_density_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to import_density_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help import_density

% Last Modified by GUIDE v2.5 06-Jun-2016 14:37:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @import_density_OpeningFcn, ...
                   'gui_OutputFcn',  @import_density_OutputFcn, ...
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


% --- Executes just before import_density is made visible.
function import_density_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to import_density (see VARARGIN)

global model

% Choose default command line output for import_density
handles.output = hObject;

n=0;
if isfield(model,'densities'),
    n=length(model.densities);
end;
defname=sprintf('density_%i',n+1);

set(handles.edit_tag,'String',defname);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes import_density wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = import_density_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



function edit_tag_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tag as text
%        str2double(get(hObject,'String')) returns contents of edit_tag as a double


% --- Executes during object creation, after setting all properties.
function edit_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tag (see GCBO)
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

tag=get(handles.edit_tag,'String');
tag=strtrim(tag);

poi=0;
overwrite=0;

id='';
if isfield(model,'density_tags'),
    id=tag2id(tag,model.density_tags);
    if isfield(model,'densities'),
        poi=length(model.densities);
        if poi==0;
            model=rmfield(model,'densities');
             model.density_tags=':';
        end;
    end;
else
    model.density_tags=':';
end;

if ~isempty(id),
    title=sprintf('Density %s already exists',tag);
    button = questdlg('Overwrite this density?',title,'No');
    switch button
        case 'No'
            return;
        case 'Cancel'
            delete(handles.figure1);
            return;
        case 'Yes'
            overwrite=1;
            id=id(1);
    end;
end;
mode=0;
if get(handles.radiobutton_MMM,'Value'),
    mode=1;
    ext='*.mat';
elseif get(handles.radiobutton_altenbach,'Value'),
    mode=2;
    ext='*.cube';
elseif get(handles.radiobutton_Kuprov,'Value'),
    mode=3;
    ext='*.mat';
end;
if ~mode,
    add_msg_board('ERROR: No valid density format selected. Nothing loaded.');
end;
[filename,pathname]=uigetfile(ext,'Import density file');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Density loading cancelled by user');
else
    fname=fullfile(pathname,filename);
    switch mode
        case 1
            load(fname);
            if exist('x','var') && exist('y','var') && exist('z','var') && exist('cube','var'),
                if ~overwrite,
                    id=poi+1;
                    model.density_tags=sprintf('%s%s:',model.density_tags,tag);
                end;
                model.densities{id}.x=x;
                model.densities{id}.y=y;
                model.densities{id}.z=z;
                model.densities{id}.tag=tag;
                model.densities{id}.cube=cube;
                add_msg_board(sprintf('Density %s loaded from file %s',tag,filename));
            else
                add_msg_board('ERROR: Selected file does not conform to MMM density format.');
            end;
        case 2
            fid = fopen( fname, 'r','ieee-be.l64') ;
            nx = fread(fid,1,'float64');
            x0 = fread(fid,1,'float64');
            dx = fread(fid,1,'float64');
            ny = fread(fid,1,'float64');
            y0 = fread(fid,1,'float64');
            dy = fread(fid,1,'float64');
            nz = fread(fid,1,'float64');
            z0 = fread(fid,1,'float64');
            dz = fread(fid,1,'float64');
            cube = fread(fid,inf,'float64');
            fclose(fid);
            cube=reshape(cube,nx,ny,nz);
            cube=permute(cube,[2,1,3]);
            x=x0:dx:x0+(nx-1)*dx;
            y=y0:dy:y0+(ny-1)*dy;
            z=z0:dz:z0+(nz-1)*dz;
            if ~overwrite,
                id=poi+1;
                model.density_tags=sprintf('%s%s:',model.density_tags,tag);
            end;
            model.densities{id}.x=x;
            model.densities{id}.y=y;
            model.densities{id}.z=z;
            model.densities{id}.tag=tag;
            model.densities{id}.cube=cube;
            model.densities{id}.snum=model.current_structure;
        case 3
            load(fname);
            if exist('ranges','var') && exist('source_cube','var'),
                if ~overwrite,
                    id=poi+1;
                    model.density_tags=sprintf('%s%s:',model.density_tags,tag);
                end;
                [nx,ny,nz] = size(source_cube);
                dx = (ranges(2)-ranges(1))/(nx-1);
                dy = (ranges(4)-ranges(3))/(ny-1);
                dz = (ranges(6)-ranges(5))/(nz-1);
                x0 = ranges(1);
                y0 = ranges(3);
                z0 = ranges(5);
                x=x0:dx:x0+(nx-1)*dx;
                y=y0:dy:y0+(ny-1)*dy;
                z=z0:dz:z0+(nz-1)*dz;
                source_cube=permute(source_cube,[2,1,3]);
                model.densities{id}.x=x;
                model.densities{id}.y=y;
                model.densities{id}.z=z;
                model.densities{id}.tag=tag;
                model.densities{id}.cube=source_cube;
                add_msg_board(sprintf('Density %s loaded from file %s',tag,filename));
            else
                add_msg_board('ERROR: Selected file does not conform to Kuprov''s density format.');
            end;

    end;
end;

delete(handles.figure1);
