function varargout = ENDOR_window(varargin)
% ENDOR_WINDOW M-file for ENDOR_window.fig
%      ENDOR_WINDOW, by itself, creates a new ENDOR_WINDOW or raises the existing
%      singleton*.
%
%      H = ENDOR_WINDOW returns the handle to a new ENDOR_WINDOW or the handle to
%      the existing singleton*.
%
%      ENDOR_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENDOR_WINDOW.M with the given input arguments.
%
%      ENDOR_WINDOW('Property','Value',...) creates a new ENDOR_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ENDOR_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ENDOR_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ENDOR_window

% Last Modified by GUIDE v2.5 14-Dec-2010 16:21:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ENDOR_window_OpeningFcn, ...
                   'gui_OutputFcn',  @ENDOR_window_OutputFcn, ...
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


% --- Executes just before ENDOR_window is made visible.
function ENDOR_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ENDOR_window (see VARARGIN)

% Choose default command line output for ENDOR_window
handles.output = hObject;

% global MMM_icon

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

handles.libration=0.5; % libration broadening of spin label position 

handles.updated=0;
handles.efficiency=1.0;
handles.center_frq=6000;
handles.center_frq_orig=6000;
handles.tau=600;
set(handles.edit_center_frq,'String',sprintf('%5i',handles.center_frq));
handles.range=[0.5,80];
handles.diamagnetic=[];
handles.labels={};
handles.expanded=false;
handles.structure=1;
handles.bsl_poly=0;

handles.trf=10;

set(handles.text_mean,'String','no distribution');
set(handles.text_stddev,'String',' ');
set(handles.text_coor,'String',' ');

handles.Pcoor=[];
handles.Pcoor_bilayer=[];
handles.Scoor=[];
handles.new_molecule=true;
handles.new_bilayer=true;
% Experimental data set defaults
handles.frq_orig=[];
handles.v_orig=[];
handles.dfrq=10;
handles.bas_name='';
handles.frq_corr=[];
handles.vexp=[];
handles.vb={};
handles.rexp=1.5:0.05:10;
handles.dsim=[];
handles.spcsim=[];
handles.frqsim=[];
handles.rmsd=[];
handles.rmean=[];
handles.stddev=[];

% kernel for fast Deer simulations
load('pake_base40_MMM.mat');
kernel=base-ones(size(base)); % kernel format for pcf2deer
handles.Pake_kernel=kernel;
handles.Pake_t=t;
handles.Pake_r=r;
handles.Pake_wd=wd;
handles.texp=t';
handles.tdip=t';

% copy control flags
handles.copy_ENDOR=0;
handles.copy_distr=0;

% experimental data names

handles.project_dir='';
handles.bas_name='';

[handles,success]=mk_label_list(handles);

load helpicon
set(handles.pushbutton_help,'CData',cdata);

if ~success,
   msgbox('Model must feature at least one spin label or selected atom','Mims ENDOR simulation impossible','error');
   ENDOR_CloseRequestFcn(handles.ENDOR, eventdata, handles); 
else
    set(handles.uipanel_diamagnetic,'SelectionChangeFcn',@radiobutton_SelectionChangeFcn);
    handles=update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;


% UIWAIT makes ENDOR_window wait for user response (see UIRESUME)
% uiwait(handles.ENDOR);


% --- Outputs from this function are returned to the command line.
function varargout = ENDOR_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton_Xepr.
function pushbutton_Xepr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Xepr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.exp_files);

ext='*.DTA;';
[fname,pname]=uigetfile(ext,'Load experimental dataset (Xepr)');
if isequal(fname,0)||isequal(pname,0),
    return; 
end;
reset_user_paths(pname);
general.exp_files=pname;
cd(pname);

% separate filename from extension
dots=findstr('.',fname);
ext_pos=length(fname);
if ~isempty(dots), ext_pos=dots(length(dots))-1; end;
bas_name=fname(1:ext_pos); % name without extension
figname=['P-31 Mims ENDOR simulation - ' fname]; % tell user, which file is current
set(handles.ENDOR,'Name',figname);
handles.project_dir=pname;
handles.bas_name=bas_name;
handles.source_file=[pname fname];

[x,y,z,vb]=get_elexsys_MMM(bas_name);

B0=get_vb2_MMM(vb,'B0VL'); % the magnetic field value is extracted from the parameter file
if ~isempty(B0),
    handles.field_known=true;
    mun=5.0507831720e-27; % nuclear magneton in SI units
    gnP31=2.26320; % nuclear g value of P-31
    h=6.6260687652e-34; % Planck's quantum of action in SI units
    Hz2kHz=1.0e-3; % conversion factor from Hz to kHz
    ny31P=Hz2kHz*B0*mun*gnP31/h; % P-31 nuclear Zeeman frequency at the given magnetic field
else
    handles.field_known=false;
    ny31P=1000*(x(1)+x(end))/2;
end;
handles.center_frq=ny31P;
handles.center_frq_orig=ny31P;
set(handles.edit_center_frq,'String',sprintf('%5.0f',ny31P));

handles.frq_orig=1000*x;
handles.v_orig=real(z);
handles.dfrq=1000*(x(2)-x(1));
handles.vb=vb;

cd(my_path);
handles=update(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_load_ASCII.
function pushbutton_load_ASCII_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_ASCII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;

cd(general.exp_files);
ext='*.dat;';
[fname,pname]=uigetfile(ext,'Load experimental dataset (ASCII)');
if isequal(fname,0)||isequal(pname,0),
    return; 
end;
reset_user_paths(pname);
general.exp_files=pname;
cd(pname);

% separate filename from extension
dots=findstr('.',fname);
ext_pos=length(fname);
if ~isempty(dots), ext_pos=dots(length(dots))-1; end;
bas_name=fname(1:ext_pos); % name without extension
figname=['P-31 Mims ENDOR simulation - ' fname]; % tell user, which file is current
set(handles.ENDOR,'Name',figname);
handles.project_dir=pname;
handles.bas_name=bas_name;
handles.source_file=[pname fname];

data=load(fname);
x=1000*data(:,1);
z=data(:,2);

ny31P=(x(1)+x(end))/2;
handles.field_known=false;
handles.center_frq=ny31P;
handles.center_frq_orig=ny31P;
set(handles.edit_center_frq,'String',sprintf('%5.0f',ny31P));

handles.frq_orig=x';
handles.v_orig=real(z)';
handles.dfrq=x(2)-x(1);

cd(my_path);

handles=update(handles);
guidata(hObject,handles);



% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

my_path=pwd;
cd(general.reports);

if get(handles.radiobutton_both,'Value'),
    mode='molecule_and_bilayer';
end;
if get(handles.radiobutton_molecule,'Value'),
    mode='molecule';
end;
if get(handles.radiobutton_bilayer,'Value'),
    mode='bilayer';
end;
name=get(handles.text_paramagnetic,'String');
labstring=sprintf('P31_Mims_%s_%s',name,mode);
suggestion=[labstring '_res.txt']; 
[fname,pname]=uiputfile('*.txt','Save results',suggestion);
if isequal(fname,0) || isequal(pname,0),
    add_msg_board('Saving of distance distribution and Mims ENDOR simulation canceled by user.');
    return;
end;
reset_user_paths(pname);
general.reports=pname;
% Remove (last) extension, if any
s=strfind(fname,'.');
if ~isempty(s),
    fname=fname(1:s(length(s))-1);
end;
% Remove suffix '_res', if present
s=strfind(fname,'_res');
if ~isempty(s),
    fname=fname(1:s(length(s))-1);
end;

fname_res=[fname '_res.txt'];

wfile=fopen(fullfile(pname,fname_res),'w+');
fprintf(wfile,'%s%s%s\n','>>> MMM 31P Mims ENDOR rotamer simulation of data set: ',handles.bas_name,' <<<');
fprintf(wfile,'%s%s%s\n\n','- ',mode,' -');
if ~isempty(handles.bas_name),
    fprintf(wfile,'%s\n','### Description of data set ###');
    fprintf(wfile,'%s\n%s\n','Source file:',handles.source_file);
    fprintf(wfile,'\n%s\n','### Pre-processing ###');
    fprintf(wfile,'%s%i\n','Background correction polynomial order: ',handles.bsl_poly);
    if get(handles.checkbox_fit,'Value'),
        fprintf(wfile,'%s\n','Center frequency fitted.');
    elseif handles.field_known && ~handles.center_edited,
        fprintf(wfile,'%s\n','Center frequency from known magnetic field.');
    elseif ~handles.center_edited
        fprintf(wfile,'%s\n','Original center frequency.');
    else
        fprintf(wfile,'%s\n','Center frequency manually edited.');
    end;
end;
fprintf(wfile,'\n%s\n','### Simulation ###');
fprintf(wfile,'%s%5.0f ns\n','Interpulse delay (manual input): ',handles.tau);
fprintf(wfile,'%s%5.0f us\n','Radiofrequency pulse length (manual input): ',handles.trf);
if ~isempty(handles.rmean),
    fprintf(wfile,'%s%4.1f ?\n','mean distance: ',handles.rmean);
end;
if ~isempty(handles.stddev),
    fprintf(wfile,'%s%4.1f ?\n','std. deviation: ',handles.stddev);
end;
if ~isempty(handles.rmsd),
    fprintf(wfile,'%s%9.6f\n','r.m.s. error of fit: ',handles.rmsd);
end;
fclose(wfile);


data=[handles.rsim',handles.dsim'];
fname_distr=[fname '_distr.dat'];
save(fullfile(pname,fname_distr),'data','-ascii');

data=[handles.frqsim',handles.spcsim'];
fname_spc=[fname '_spc.dat'];
save(fullfile(pname,fname_spc),'data','-ascii');

if ~isempty(handles.vexp),
    add_msg_board('Simulation and fit to experimental data saved');
    data=[handles.frq_corr',handles.vexp',handles.vsim'];
    fname_fit=[fname '_fit.dat'];
    save(fullfile(pname,fname_fit),'data','-ascii');
end;

cd(my_path);

guidata(hObject,handles);

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ENDOR_CloseRequestFcn(handles.ENDOR, eventdata, handles);


% --- Executes on selection change in listbox_label.
function listbox_label_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_label contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_label

global model

sel=get(hObject,'Value');
if isfield(model,'labels') && sel<=length(model.labels),
    NOpos=model.labels(sel).NOpos;
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    name=model.labels(sel).adr;
    indices=resolve_address(name);
    position=sprintf('at [%6.2f,%6.2f,%6.2f] ?',x,y,z);
    handles.Scoor=[NOpos(:,1:3),pop];
else
    if isfield(model,'labels'),
        asel=sel-length(model.labels);
    else
        asel=sel;
    end;
    indices=handles.atoms(asel,:);
    [msg0,xyz]=get_atom(indices,'xyz');
    [msg0,pop]=get_atom(indices,'populations');
    handles.Scoor=[xyz,pop];
    pop=pop/sum(pop);
    xyz=pop'*xyz;
    adr=mk_address(indices);
    name=adr;
    position=sprintf('at [%6.2f,%6.2f,%6.2f] ?',xyz(1),xyz(2),xyz(3));
end;
set(handles.text_coor,'String',position);
set(handles.text_paramagnetic,'String',name);
old_snum=handles.structure;
handles.structure=indices(1);
if indices(1)~=old_snum,
    bilayer_hole(indices(1));
end;
handles=get_phosphorous(handles);
handles.new_molecule=true;
handles.new_bilayer=true;
handles=update(handles);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function listbox_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'ENDOR_window.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_tau_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tau as text
%        str2double(get(hObject,'String')) returns contents of edit_tau as a double

[v,handles]=edit_update_MMM(handles,hObject,0,20000,600,'%4.0f',0);
handles.tau=v;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_fit.
function checkbox_fit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fit

if ~get(hObject,'Value'),
    handles.center_frq=handles.center_frq_orig;
    set(handles.edit_center_frq,'String',sprintf('%5.0f',handles.center_frq));
end;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_detach_ENDOR.
function pushbutton_detach_ENDOR_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_ENDOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy_ENDOR=1;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_detach.
function pushbutton_detach_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy_distr=1;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_lower_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower as text
%        str2double(get(hObject,'String')) returns contents of edit_lower as a double

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
else
    [v,handles]=edit_update_MMM(handles,hObject,min(handles.rsim),handles.range(2)-0.1,handles.range(1),'%4.1f',0);
    handles.range=[v,handles.range(2)];
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes during object creation, after setting all properties.
function edit_lower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_lower_left.
function pushbutton_lower_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lower_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(1)-0.1>=min(handles.rsim)
    handles.range(1)=handles.range(1)-0.1;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;


% --- Executes on button press in pushbutton_lower_right.
function pushbutton_lower_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lower_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(1)+0.1<=handles.range(2)-0.1,
    handles.range(1)=handles.range(1)+0.1;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

function edit_upper_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_upper as text
%        str2double(get(hObject,'String')) returns contents of edit_upper as a double

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
else
    [v,handles]=edit_update_MMM(handles,hObject,handles.range(1)+0.1,max(handles.rsim),handles.range(2),'%4.1f',0);
    handles.range=[handles.range(1),v];
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes during object creation, after setting all properties.
function edit_upper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_upper_left.
function pushbutton_upper_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_upper_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(2)-0.1>=handles.range(1)+0.1
    handles.range(2)=handles.range(2)-0.1;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes on button press in pushbutton_upper_right.
function pushbutton_upper_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_upper_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.rsim),
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(2)+0.1<=max(handles.rsim)
    handles.range(2)=handles.range(2)+0.1;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

function handles=range_update(handles)

if isfield(handles,'left_crsr'),
    set(handles.left_crsr,'XData',[handles.range(1) handles.range(1)]);
else
    axes(handles.axes_distr);
    handles.left_crsr=plot([handles.range(1) handles.range(1)],[min(handles.dsim),max(handles.dsim)],'b','LineWidth',0.5);
end;
set(handles.edit_lower,'String',sprintf('%4.1f',handles.range(1)));

if isfield(handles,'right_crsr'),
    set(handles.right_crsr,'XData',[handles.range(2) handles.range(2)]);
else
    axes(handles.axes_distr);
    handles.right_crsr=plot([handles.range(2) handles.range(2)],[min(handles.dsim),max(handles.dsim)],'m','LineWidth',0.5);
end;
set(handles.edit_upper,'String',sprintf('%4.1f',handles.range(2)));



% --- Executes on button press in togglebutton_expand.
function togglebutton_expand_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_expand

state=get(hObject,'Value');
handles.expanded=state;
if state,
    set(hObject,'TooltipString','Return to full view');
    set(hObject,'String','Full view');
else
    set(hObject,'TooltipString','Expand or shrink view to selected range');
    set(hObject,'String','Expand');
end;

handles=update(handles);
% Update handles structure
guidata(hObject, handles);

function edit_center_frq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_center_frq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_center_frq as text
%        str2double(get(hObject,'String')) returns contents of edit_center_frq as a double

[v,handles]=edit_update_MMM(handles,hObject,1000,100000,6000,'%5.0f',0);
handles.center_frq=v;
handles.center_frq_orig=v;
handles.center_edited=true;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_center_frq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_center_frq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close ENDOR.
function ENDOR_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to ENDOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

function [handles,success]=mk_label_list(handles)
% Fills listbox_label with the list of available labels
% returns a success flag success=1 if there are at least two labels,
% success=0 otherwise

global model

success=0;

indices=resolve_address('*');

if ~isfield(model,'labels')
    labels=0;
    if isempty(indices),
        return
    end;
end;

sites=0;

if isfield(model,'labels') && ~isempty(model.labels),
    sites=length(model.labels);
    labels=sites;
    for k=1:length(model.labels),
        label_list{k}=model.labels(k).adr;
    end;
end;


% make list of selected atoms
if ~isempty(indices),
    [m,n]=size(indices);
    poi=0;
    handles.atoms=zeros(m,5);
    for ko=1:m, % loop over all objects
        idepth=length(find(indices(ko,:)>0)); % determine type of current object
        if idepth==5,
            poi=poi+1;
            cindices=indices(ko,1:idepth);
            address=mk_address(cindices,true);
            label_list{poi+sites}=address;
            handles.atoms(poi,:)=cindices;
            handles.atom_adr{poi}=address;
        end;
    end;
    if poi>0,
        handles.atoms=handles.atoms(1:poi,:);
        sites=sites+poi;
    else
        handles.atoms=[];
    end;
end;

if sites>=1,
    success=1;
else
    set(handles.listbox_label,'String',' ');
    set(handles.listbox_label,'Value',1);
    return;
end;

set(handles.listbox_label,'String',label_list);
set(handles.listbox_label,'Value',1);

if labels>0,
    NOpos=model.labels(1).NOpos;
    pop=NOpos(:,4);
    pop=pop/sum(pop);
    handles.Scoor=[NOpos(:,1:3),pop];
    x=sum(NOpos(:,1).*pop);
    y=sum(NOpos(:,2).*pop);
    z=sum(NOpos(:,3).*pop);
    name=model.labels(1).adr;
    indices=resolve_address(name);
    position=sprintf('at [%6.2f,%6.2f,%6.2f] ?',x,y,z);
else
    indices=handles.atoms(1,:);
    [msg0,xyz]=get_atom(indices,'xyz');
    [msg0,pop]=get_atom(indices,'populations');
    pop=pop/sum(pop);
    handles.Scoor=[xyz,pop];
    xyz=pop'*xyz;
    adr=mk_address(indices);
    name=adr;
    position=sprintf('at [%6.2f,%6.2f,%6.2f] ?',xyz(1),xyz(2),xyz(3));    
end;
set(handles.text_coor,'String',position);
set(handles.text_paramagnetic,'String',name);
handles.structure=indices(1);
bilayer_hole(indices(1));
handles=get_phosphorous(handles);

function handles=get_phosphorous(handles)

global model

grid=0.5;

h=gcf;
set(h,'Pointer','watch');
drawnow;

% get structural phosphorous from molecule
adr=mk_address(handles.structure);
adr=sprintf('%s(:)."P"',adr);
indices=resolve_address(adr);
[m,n]=size(indices);
Pcoor=zeros(4*m,4);
poi=0;
for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    [mess,xyz]=get_object(cindices,'xyz');
    [ml,n]=size(xyz);
    [mess,pop]=get_object(cindices,'populations');
    Pcoor(poi+1:poi+ml,1:3)=xyz;
    Pcoor(poi+1:poi+ml,4)=pop;
    poi=poi+ml;
end;
handles.Pcoor=Pcoor(1:poi,:);
msg{1}=sprintf('%i phosphorous nuclei in this structure',m);

snum=handles.structure;
if ~isfield(model.info{snum},'bilayer') || isempty(model.info{snum}.bilayer),
    msg{2}='No lipid bilayer present.';
    handles.Pcoor_bilayer=[];
else
    xax=model.info{snum}.bilayer.xax;
    nxm=length(xax)-1;
    xmin=min(xax);
    xmax=max(xax);
    dx=nxm/(xmax-xmin);
    yax=model.info{snum}.bilayer.yax;
    nym=length(yax)-1;
    ymin=min(yax);
    ymax=max(yax);
    dy=nym/(ymax-ymin);
    upper_layer=model.info{snum}.bilayer.upper_layer;
    lower_layer=model.info{snum}.bilayer.lower_layer;
    Scoor=handles.Scoor;
    xyz=Scoor(:,1:3);
    pop=Scoor(:,4);
    center=pop'*xyz;
    xax1=-25:grid:25;
    yax1=-25:grid:25;
    zax1=-25:grid:25;
    poi=0;
    Pcoor_bilayer=zeros(250000,4);
    [zax,profile,A_per_lipid]=bilayer_P_profile;
    totpop=sum(profile); % normalization of z profile of 31P distribution
    pixel_area=grid^2; % pixel size (area of grid square)
    pixel_per_lipid=pixel_area/A_per_lipid; % number of pixels in area per lipid
    popscale=pixel_per_lipid/totpop; % normalization to one 31P per lipid area integrated over z
    thickness=model.info{snum}.bilayer.thickness;
    zax=zax*thickness/35;
    nzm=length(zax)-1;
    zmin=min(zax);
    zmax=max(zax);
    dz=nzm/(zmax-zmin);
    for kx=1:length(xax1),
        x=xax1(kx);
        for ky=1:length(yax1),
            y=yax1(ky);
            if norm([x,y])<=25,
                for kz=1:length(zax1),
                    z=zax1(kz);
                    if norm([x,y,z])<=25,
                        xyz=[x,y,z]+center;
                        if xyz(1)>xmin && xyz(1)<xmax && xyz(2)>ymin && xyz(2)<ymax,
                            poix=1+round(dx*(xyz(1)-xmin));
                            poiy=1+round(dy*(xyz(2)-ymin));
                            upop=upper_layer(poix,poiy);
                            lpop=lower_layer(poix,poiy);
                        else
                            upop=1;
                            lpop=1;
                        end;
                        if abs(xyz(3))<zmax,
                            poiz=1+round(dz*(abs(xyz(3))-zmin));
                            if xyz(3)>0,
                                pop=upop*profile(poiz);
                            else
                                pop=lpop*profile(poiz);
                            end;
                        else
                            pop=0;
                        end;
                        if pop>0.01,
                            poi=poi+1;
                            Pcoor_bilayer(poi,1:3)=xyz;
                            Pcoor_bilayer(poi,4)=pop*popscale;
                        end;
                    end;
                end;
            end;
        end;
    end;
    handles.Pcoor_bilayer=Pcoor_bilayer(1:poi,:);
    msg{2}=sprintf('%i phosphorous positions for lipid bilayer',poi);
end;

set(handles.text_phosphorous,'String',msg);
set(h,'Pointer','arrow');


function handles=update(handles)
% Display update

hgui=gcf;
set(hgui,'Pointer','watch');
drawnow;

if handles.copy_distr,
    figure(1); clf;
    set(gca,'FontSize',16);
else
    axes(handles.axes_distr);
    cla;
    if isfield(handles,'left_crsr'),
        handles=rmfield(handles,'left_crsr');
    end;
    if isfield(handles,'right_crsr'),
        handles=rmfield(handles,'right_crsr');
    end;
end;
hold on;

poi=0;
if isempty(handles.Pcoor),
    molecule=false;
else
    molecule=get(handles.radiobutton_molecule,'Value')+get(handles.radiobutton_both,'Value');
end;
if molecule,
    poi=poi+1;
end;
if isempty(handles.Pcoor_bilayer),
    bilayer=false;
else
    bilayer=get(handles.radiobutton_bilayer,'Value')+get(handles.radiobutton_both,'Value');
end;
if bilayer,
    poi=poi+1;
end;
rmin0=1e6;
rmax0=0;
% Display distance distributions
if molecule,
    if handles.new_molecule,
        [rax,act_distr]=get_distribution(10*handles.Scoor,10*handles.Pcoor,handles.libration,0,false,25);
        handles.molecule_rax=rax;
        handles.molecule_distr=act_distr;
        handles.new_molecule=false;
    else
        rax=handles.molecule_rax;
        act_distr=handles.molecule_distr;
    end;
    if sum(act_distr)<=1e-12*length(act_distr),
        act_distr=1e-12*ones(size(act_distr));
    else
        sc=sum(handles.Scoor(:,4));
        act_distr=act_distr/sc;
    end;
    h=plot(rax,act_distr,'b');
    full_distr=act_distr;
    simdistr=true;
end;
if bilayer,
    if handles.new_bilayer,
        [rax,act_distr]=get_distribution(10*handles.Scoor,10*handles.Pcoor_bilayer,handles.libration,0,false,25);
        handles.bilayer_rax=rax;
        handles.bilayer_distr=act_distr;
        handles.new_bilayer=false;
    else
        rax=handles.bilayer_rax;
        act_distr=handles.bilayer_distr;
    end;
    if sum(act_distr)<=1e-12*length(act_distr),
        act_distr=1e-12*ones(size(act_distr));
    else
        sc=sum(handles.Scoor(:,4));
        act_distr=act_distr/sc;
    end;
    h=plot(rax,act_distr,'g');
    if ~molecule,
        full_distr=act_distr;
    else
        full_distr=full_distr+act_distr;
    end;
    simdistr=true;
end;
if ~molecule && ~bilayer,
    add_msg_board('Warning: No phosphorous nuclei in this structure for this selection. Empty distance distribution.');
    rax=get_distribution([10,0,0,1],[0,0,0,1],0.05);
    rax=10*rax;
    full_distr=1e-12*ones(size(rax));
    simdistr=false;
end;

handles.rsim=rax;
handles.dsim=full_distr;
ndistr=full_distr/sum(full_distr);
rmean=sum(ndistr.*rax);
set(handles.text_mean,'String',sprintf('<r> = %4.1f ?',rmean));
handles.rmean=rmean;
dr=rax-rmean;
stddev=sqrt(sum(dr.^2.*ndistr));
set(handles.text_stddev,'String',sprintf('sr = %4.1f ?',stddev));
handles.stddev=stddev;
rfilled=rax(full_distr>=0.01*max(full_distr));
rmin=min([min(rfilled) rmin0]);
rmax=max([max(rfilled) rmax0]);
bsl=0.1*(rmax-rmin);
rmin=rmin-bsl;
rmax=rmax+bsl;
dmax=1.1*max(full_distr);
plot(rax,full_distr,'r','LineWidth',1);
dmin=-0.09*dmax;
if ~handles.expanded,
    if rmin> handles.range(1) || rmax< handles.range(2),
        handles.range(1)=rmin;
        handles.range(2)=rmax;
    end;
end;
handles=range_update(handles);
if handles.expanded,
    axis([handles.range(1),handles.range(2),dmin,dmax]);
else
    axis([rmin,rmax,dmin,dmax]);
end;

if handles.copy_ENDOR,
    figure(2); clf;
    set(gca,'FontSize',16);
    hold on;
else
    axes(handles.axes_ENDOR);
    cla;
    hold on;
end;

[frqax,spc]=mk_ENDOR_spc(handles,rax,full_distr);
max_effect=1-min(spc);
half_reduction=0.985;
poimin=find(spc<half_reduction,1);
poimax=find(spc<half_reduction,1,'last');
if isempty(poimin),
    set(handles.text_FWHR,'String','n.a.');
else
    FWHR=frqax(poimax)-frqax(poimin);
    set(handles.text_FWHR,'String',sprintf('%3.0f',FWHR));
end;
set(handles.text_max_effect,'String',sprintf('%5.2f',100*max_effect));
handles.frqsim=frqax;
handles.spcsim=spc;

plot(frqax,spc,'r');

frqfilled=frqax(1-spc>=0.002*(1-min(spc)));
frqmin=min(frqfilled);
frqmax=-frqmin;
bsl=0.1*(frqmax-frqmin);
if bsl<1,
    bsl=500;
end;
frqmin=frqmin-bsl;
frqmax=frqmax+bsl;
if max(spc)-min(spc)<1e-3,
    dspc=1e-3;
else
    dspc=max(spc)-min(spc);
end;
axis([frqmin,frqmax,1-1.1*dspc,1+0.1*dspc]);

if ~isempty(handles.v_orig),
    ndat=length(handles.v_orig);
    nbsl=round(ndat/10);
    spare=round(nbsl/4);
    xbsl=[handles.frq_orig(1+spare:nbsl) handles.frq_orig(end-nbsl+1:end)];
    bsl=[handles.v_orig(1+spare:nbsl) handles.v_orig(end-nbsl+1:end)];
    p=polyfit(xbsl,bsl,handles.bsl_poly);
    bslcorr=polyval(p,handles.frq_orig);
    handles.vexp=handles.v_orig./bslcorr;
    if get(handles.checkbox_fit,'Value'),
        v0=handles.center_frq;
        handles.center_frq=fminsearch(@rmsd_frq_corr,v0,[],handles.frq_orig,handles.vexp,frqax,spc);
        set(handles.edit_center_frq,'String',sprintf('%5.0f',handles.center_frq));
    end;
    handles.frq_corr=handles.frq_orig-handles.center_frq;
%     store=gca;
%     figure(13); clf;
%     plot(handles.frq_orig,handles.v_orig,'k');
%     figure(7); clf;
%     plot(handles.frq_corr,handles.vexp,'k');
%     axes(store);
    plot(handles.frq_corr,handles.vexp,'k');
    sim=interp1(frqax,spc,handles.frq_corr,'pchip',0);
    diff=handles.vexp-sim;
    rmsd=sqrt(sum(diff.^2)/length(diff));
    handles.rmsd=rmsd;
    set(handles.text_rmsd,'String',sprintf('%7.5f',rmsd));
    handles.vsim=sim;
    vaxis=axis(gca);
    dfrq=frqmax-frqmin;
    if min(handles.frq_corr)-0.1*frqmin<vaxis(1), 
        vaxis(1)=min(handles.frq_corr)-0.1*frqmin; 
    end;
    if max(handles.frq_corr)+0.1*frqmin>vaxis(2), 
        vaxis(2)=max(handles.frq_corr)+0.1*frqmin; 
    end;
    dexp=max(handles.vexp)-min(handles.vexp);
    if min(handles.vexp)-0.1*dexp<vaxis(3),
        vaxis(3)=min(handles.vexp)-0.1*dexp;
    end;
    if max(handles.vexp)+0.1*dexp>vaxis(4),
        vaxis(4)=max(handles.vexp)+0.1*dexp;
    end;
    axis(vaxis);
end;

handles.copy_distr=0;
handles.copy_ENDOR=0;

set(hgui,'Pointer','arrow');
drawnow;

function [frqax,spc]=mk_ENDOR_spc(handles,rax,distr)
% Simulates 31P-Mims-ENDOR spectrum from distance axis rax, distance
% distribution distr, and first interpulse delay tau (supplied in
% handles.tau)
% distance axis is supplied in Angstroem!
% this corresponds to having a frequency in kHz instead of MHz (time axis
% in ms)

reduce=-6.15575e-4; % ratio g_n mu_n/ g_e mu_B for 31P

pcf=get_std_distr_MMM(rax,distr,handles.Pake_r);
depth=sum(distr);
ff=get_formfactor_MMM(pcf,handles.Pake_kernel,handles.Pake_t);
if depth<1e-3,
    ff=zeros(size(ff));
else
    ff=depth*ff;
end;

tdip=handles.Pake_t;
fill=4;

tdmp=2*tdip/(1000*abs(reduce));
dmp=exp(-(tdmp/handles.trf).^2);
ffdmp=dmp.*ff;

dipevo2=[private_hamming(ones(size(ff)),1) zeros(1,(fill-1)*length(ff))];
dipevo2(1)=dipevo2(1)/2;
spc=real(fftshift(fft(dipevo2)));
corr=1/max(spc); % correction factor

dipevo2=[private_hamming(ffdmp,1) zeros(1,(fill-1)*length(ff))];
dipevo2(1)=dipevo2(1)/2;
spc=real(fftshift(fft(dipevo2)));
spc=corr*spc;
fmin=-length(dipevo2)/(2*fill*max(tdip));          % s. Schweiger/Jeschke book p. 106
fmax=(length(dipevo2)-1)/(2*fill*max(tdip));
% frequency axis divided by factor 2, as DEER detects splittings, 
% but ENDOR detects frequencies
% multiplied by factor 1000 since distance axis is in ? instead of nm
frqax=1000*abs(reduce)*linspace(fmin,fmax,length(spc))/2; 
% and now go to kHz frequency axis
frqax=1000*frqax;

% remove insignificant contributions
% spc(spc<0.002*max(spc))=0;

tau=handles.tau*1e-6; % going from ns to ms, corresponding to kHz frequency axis
% compute blindspot function for experimental tau value
blindspot=sin(2*pi*frqax*tau); % sin(2pi ny tau) 
blindspot=blindspot.^2; % sin^2(2pi ny tau)
spc=(spc+fliplr(spc))/2; % symmetrize spectrum
blindspot=(blindspot+fliplr(blindspot))/2; % symmetrize blindpsot pattern
spc=spc.*blindspot;
spc=1-spc;



function edit_trf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trf as text
%        str2double(get(hObject,'String')) returns contents of edit_trf as a double

[v,handles]=edit_update_MMM(handles,hObject,1,100,10,'%4.1f',0);
handles.trf=v;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_trf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rmsd=rmsd_frq_corr(v,frq_orig,vexp,frqax,spc)

frq_corr=frq_orig-v;
sim=interp1(frqax,spc,frq_corr,'pchip',0);
diff=vexp-sim;
rmsd=sqrt(sum(diff.^2)/length(diff));



function edit_bckg_poly_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bckg_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bckg_poly as text
%        str2double(get(hObject,'String')) returns contents of edit_bckg_poly as a double

[v,handles]=edit_update_MMM(handles,hObject,0,5,0,'%i',1);
handles.bsl_poly=v;
handles=update(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_bckg_poly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bckg_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_backg_order_left.
function pushbutton_backg_order_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_backg_order_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.bsl_poly>0
    handles.bsl_poly=handles.bsl_poly-1;
    set(handles.edit_bckg_poly,'String',sprintf('%i',handles.bsl_poly));
    handles=update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;

% --- Executes on button press in pushbutton_bckg_order_right.
function pushbutton_bckg_order_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_bckg_order_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.bsl_poly<6
    handles.bsl_poly=handles.bsl_poly+1;
    set(handles.edit_bckg_poly,'String',sprintf('%i',handles.bsl_poly));
    handles=update(handles);
    % Update handles structure
    guidata(hObject, handles);
end;


% --------------------------------------------------------------------
function uipanel_diamagnetic_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_diamagnetic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_both.
function radiobutton_SelectionChangeFcn(source, eventdata)
% hObject    handle to radiobutton_both (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(source);
handles=update(handles);
% Update handles structure
guidata(source, handles);

function [zax,profile,area]=bilayer_P_profile

load bilayer_definitions
lip=DPPC;
zax=linspace(0,35,351);
VPCN=lip.RPCN*lip.VHL;
cPCN=VPCN/(lip.A*lip.sigPCN);
profile=Pgauss(zax,cPCN,lip.zPCN,lip.sigPCN);
area=lip.A;

function Pi=Pgauss(z,ci,zi,sigi)
% double Gaussian distribution, Eq. (2) of Kucerka et al.
Pi=ci*(exp(-(z+zi).^2/(2*sigi^2))+exp(-(z-zi).^2/(2*sigi^2)))/sqrt(2*pi);
