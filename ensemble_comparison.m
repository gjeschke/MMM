function varargout = ensemble_comparison(varargin)
% ENSEMBLE_COMPARISON MATLAB code for ensemble_comparison.fig
%      ENSEMBLE_COMPARISON, by itself, creates a new ENSEMBLE_COMPARISON or raises the existing
%      singleton*.
%
%      H = ENSEMBLE_COMPARISON returns the handle to a new ENSEMBLE_COMPARISON or the handle to
%      the existing singleton*.
%
%      ENSEMBLE_COMPARISON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENSEMBLE_COMPARISON.M with the given input arguments.
%
%      ENSEMBLE_COMPARISON('Property','Value',...) creates a new ENSEMBLE_COMPARISON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ensemble_comparison_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ensemble_comparison_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ensemble_comparison

% Last Modified by GUIDE v2.5 12-Dec-2019 10:21:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ensemble_comparison_OpeningFcn, ...
                   'gui_OutputFcn',  @ensemble_comparison_OutputFcn, ...
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


% --- Executes just before ensemble_comparison is made visible.
function ensemble_comparison_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ensemble_comparison (see VARARGIN)

% Choose default command line output for ensemble_comparison
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ensemble_comparison wait for user response (see UIRESUME)
% uiwait(handles.figure_ensemble_comparison);


% --- Outputs from this function are returned to the command line.
function varargout = ensemble_comparison_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
