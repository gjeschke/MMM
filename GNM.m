function varargout = GNM(varargin)
% GNM M-file for GNM.fig
%      GNM, by itself, creates a new GNM or raises the existing
%      singleton*.
%
%      H = GNM returns the handle to a new GNM or the handle to
%      the existing singleton*.
%
%      GNM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GNM.M with the given input arguments.
%
%      GNM('Property','Value',...) creates a new GNM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GNM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GNM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GNM

% Last Modified by GUIDE v2.5 20-Oct-2010 10:43:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GNM_OpeningFcn, ...
                   'gui_OutputFcn',  @GNM_OutputFcn, ...
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


% --- Executes just before GNM is made visible.
function GNM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GNM (see VARARGIN)

global MMM_icon
global hMain
global model
global ENM_param

% Choose default command line output for GNM
handles.output = hObject;

set(hMain.figure,'Pointer','watch');

drawnow;

h = msgbox('Please be patient. This can take several minutes.','Gaussian network model is computed');

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

hMain.auxiliary=[hMain.auxiliary hObject];

greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
greyscale=greyscale';
handles.greyscale=flipud(greyscale);
tempscale=[ones(1,51),linspace(1,0,50);linspace(0,1,50),1,linspace(1,0,50);linspace(0,1,50),ones(1,51)];
handles.tempscale=flipud(tempscale');
handles.saturation=8;

handles.crs1=[];
handles.crs2=[];

if ~exist('model','var') || ~isfield(model,'current_structure'),
    add_msg_board('Error: Cannot set up GNM, since no current structure exists.');
    guidata(hObject,handles);
    if ishandle(h),
        delete(h);
    end;
    set(hMain.figure,'Pointer','arrow');
    delete(hObject);
    return;
end;

snum=model.current_structure;
adr=mk_address(snum);
set(handles.figure1,'Name',sprintf('Gaussian network model for structure %s',adr));

if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
    [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues('!');
    if isempty(Ca_coor) || length(masses)<2,
        set(hMain.figure,'Pointer','arrow');
        add_msg_board('Error: Cannot set up GNM, since current structure has less than two amino acid residues.');
        if ishandle(h),
            delete(h);
        end;        
        guidata(hObject,handles);
        delete(hObject);
        return
    end;
    [m,n]=size(Ca_coor);
    resnums=zeros(1,m);
    for k=1:m,
        [stag,ctag,modelnum,resnum]=mk_address_parts(rindices(k,:));
        resnums(k)=resnum;
    end;
    disp(sprintf('Model has %i residues',m));
    model.coarse(snum).Ca_coor=Ca_coor;
    model.coarse(snum).indices=rindices;
    model.coarse(snum).masses=masses;
    model.coarse(snum).Bfactors=Bfactors;
    model.coarse(snum).restypes=restypes;
    model.coarse(snum).resnums=resnums;
    % make and store chain assignment
    [mr,nr]=size(rindices); 
    poi=1;
    chains=zeros(50,3);
    currchain=rindices(1,2);
    chains(1,1)=currchain;
    chains(1,2)=1;
    for k=2:mr,
        if rindices(k,2)~=currchain,
            currchain=rindices(k,2);
            chains(poi,3)=k-1;
            poi=poi+1;
            chains(poi,1)=currchain;
            chains(poi,2)=k;
        end;
    end;
    chains(poi,3)=mr;
    chains=chains(1:poi,:);
    model.coarse(snum).chains=chains;
else
    Ca_coor=model.coarse(snum).Ca_coor;
    restypes=model.coarse(snum).restypes;
    rindices=model.coarse(snum).indices;
    if isfield(model.coarse(snum),'resnums'),
        resnums=model.coarse(snum).resnums;
    else
        [m,n]=size(Ca_coor);
        resnums=zeros(1,m);
        for k=1:m,
            [stag,ctag,modelnum,resnum]=mk_address_parts(rindices(k,:));
            resnums(k)=resnum;
        end;
        model.coarse(snum).resnums=resnums;
    end;

end;


if ~isfield(model,'GNM') || length(model.GNM)<snum || isempty(model.GNM(snum).u),
   Gamma=setup_GNM_poly(Ca_coor);
    contacts=[];
%     [Gamma,contacts]=setup_GNM(Ca_coor);
    % [Gamma,contacts]=setup_GNM(Ca_coor,restypes,rindices); % residue-specific force constants, not recommended
    
    [u,D]=eig(Gamma);
    [m,n]=size(Gamma);

    lambda=zeros(1,m);
    for k=2:m,
        lambda(k)=D(k,k);
    end;
    model.GNM(snum).lambda=lambda;
    model.GNM(snum).u=u;
    model.GNM(snum).residues=m;
    model.GNM(snum).contacts=contacts;
    D(1,1)=0;
    fluc2=u*pinv(D)*u';
    model.GNM(snum).fluctuations=fluc2;
    msf=zeros(1,m);
    for k=1:m,
        msf(k)=3*fluc2(k,k)/ENM_param.gamma;
        fluc2(k,k)=0;
    end;
    ma=max(max(fluc2));
    for k=1:m,
        fluc2(k,k)=ma;
    end;
    model.GNM(snum).msf=msf;
   %  model.GNM(snum).fluctuations=fluc2;
else
    lambda=model.GNM(snum).lambda;
    m=model.GNM(snum).residues;
    fluc2=model.GNM(snum).fluctuations;
    for k=1:m,
        fluc2(k,k)=0;
    end;
    ma=max(max(fluc2));
    for k=1:m,
        fluc2(k,k)=ma;
    end;
end;

handles.xa=1;
handles.ya=1;
handles.xe=m;
handles.ye=m;
handles.mode=1;
handles.crsr_right=[];
handles.crsr_left=[];
handles.range=[];

axes(handles.axes_plot);
plot(lambda(2:end));
xlabel('Mode number');
ylabel('Normalized square frequency');

fluc2=255*fluc2/max(max(fluc2));

axes(handles.axes_fluctuations);
handles.corr_plot=image(fluc2);
set(gca,'YDir','normal');
set(handles.corr_plot,'ButtonDownFcn',@corr_click);

axis equal
axis tight
colormap(handles.greyscale.^handles.saturation);
hold on;

if ishandle(h),
    delete(h);
end;

set(hMain.figure,'Pointer','arrow');

set(gcf,'WindowButtonUpFcn',@my_button_up)
set(gcf,'WindowButtonMotionFcn',@my_motion_fct)
 
handles.crsr_mode=0; % selection mode

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GNM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GNM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'GNM.html');
webcall(entry,'-helpbrowser');

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);


% --- Executes on selection change in popupmenu_plot.
function popupmenu_plot_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot

handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);

function handles=update_auxiliary_plot(hObject,handles)

global model

contents = get(handles.popupmenu_plot,'String');
type=contents{get(handles.popupmenu_plot,'Value')};

axes(handles.axes_plot);
cla;

switch type
    case 'Mode distribution'
        plot(model.GNM(model.current_structure).lambda);
        xlabel('Mode number');
        ylabel('Normalized square frequency');
    case 'B factors'
        Bfac=model.coarse(model.current_structure).Bfactors;
        h1=plot(Bfac,'k');
        hold on;
        [msf,p,rmsd,r]=lin_fit(8*pi^2*model.GNM(model.current_structure).msf,model.coarse(model.current_structure).Bfactors);
        h2=plot(msf,'r');
        xlabel('Residue number');
        ylabel('B factor');
        handles.ya_aux=min([min(Bfac) min(msf)]);
        handles.ye_aux=max([max(Bfac) max(msf)]);
        set(handles.axes_plot,'ButtonDownFcn',@(hObject,eventdata)linked_cursor);
        set(h1,'ButtonDownFcn',@(hObject,eventdata)linked_cursor);
        set(h2,'ButtonDownFcn',@(hObject,eventdata)linked_cursor);
        msg=sprintf('rmsd: %5.2f, corr. coeff. %4.2f, offset: %5.2f, scale: %4.2f',rmsd,r,p(2),p(1));
        set(handles.text_selected_2,'String',msg);
    case 'Single mode'
        mode=model.GNM(model.current_structure).u(:,handles.mode+1);
        h1=plot(mode,'k');
        hold on
        h2=plot([1,length(mode)],[0,0],'g');
        xlabel('Residue number');
        ylabel('Mode eigenvector');
        handles.ya_aux=min(mode);
        handles.ye_aux=max(mode);
        set(handles.axes_plot,'ButtonDownFcn',@(hObject,eventdata)linked_cursor);
        set(h1,'ButtonDownFcn',@(hObject,eventdata)linked_cursor);
        set(h2,'ButtonDownFcn',@(hObject,eventdata)linked_cursor);
    case 'Contact numbers'
        contacts=model.GNM(model.current_structure).contacts;
        plot(contacts,'k.');
        hold on;
        plot(contacts(contacts==0),'ro');
end;

% --- Executes during object creation, after setting all properties.
function popupmenu_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes_fluctuations_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_fluctuations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function corr_click(hObject, eventdata)

global model
global hMain

attract=0; % range around clicked point in which residue pair with maximum correlation is searched

handles=guidata(hObject);

old_units=get(handles.axes_fluctuations,'Units');
set(handles.axes_fluctuations,'Units','pixels');
point3=get(handles.axes_fluctuations,'Position');
set(handles.axes_fluctuations,'Units',old_units);

cp = get(handles.axes_fluctuations,'CurrentPoint');
x=round(cp(1,1));
y=round(cp(1,2));

if handles.crsr_mode==1, % zoom
    handles.xc=x;
    handles.yc=y;
    handles.zwindow=line('XData',x,'YData',y,'color','b');
    set(handles.axes_fluctuations,'DrawMode','fast');
    guidata(hObject,handles);
    return
end;

x1=x;
y1=y;
corr=-1;
for k1=x-attract:x+attract,
    for k2=y-attract:y+attract,
        if k1>=1 && k1<=model.GNM(model.current_structure).residues && k2>=1 && k2<=model.GNM(model.current_structure).residues,
            if model.GNM(model.current_structure).fluctuations(k1,k2)>corr,
                corr=model.GNM(model.current_structure).fluctuations(k1,k2);
                x1=k1;
                y1=k2;
            end;
        end;
    end;
end;

ind1=model.coarse(model.current_structure).indices(x1,:);
ind2=model.coarse(model.current_structure).indices(y1,:);
adr1=mk_address(ind1);
p=strfind(adr1,']');
adr1=adr1(p+1:end);
adr2=mk_address(ind2);
p=strfind(adr2,']');
adr2=adr2(p+1:end);
col='r';
if handles.crsr_mode==0,
    hMain=cmd(hMain,sprintf('select %s',adr1));
    hMain=cmd(hMain,sprintf('select %s',adr2));
    col=[1,0.75,0];
end;
set(handles.text_selected,'String',sprintf('%s; %s = %6.4f',adr1,adr2,model.GNM(model.current_structure).fluctuations(x1,y1)));


axes(handles.axes_fluctuations);
if ishandle(handles.crs1),
    delete(handles.crs1);
end;
if ishandle(handles.crs2),
    delete(handles.crs2);
end;

xcl=4*(handles.xe-handles.xa)/point3(3);
ycl=4*(handles.ye-handles.ya)/point3(4);
handles.crs1=plot([x1-xcl,x1+xcl],[y1,y1],'color',col);
handles.crs2=plot([x1,x1],[y1-ycl,y1+ycl],'color',col);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_darker.
function pushbutton_darker_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_darker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_fluctuations);
handles.saturation=2*handles.saturation;
handles=update_correlation_plot(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_lighter.
function pushbutton_lighter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lighter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_fluctuations);
handles.saturation=handles.saturation/2;
handles=update_correlation_plot(handles);
guidata(hObject,handles);

function [sc,rmsd]=optimum_scale(y,x)
% Optimum scaling of data y_i, so that the root mean square deviation
% rmsd = sum((sc*y_i-x_i)^2) is minimized
% vectors x and y must have the same length
%
% y     vector of data to be scaled
% x     vector of data that should be fitted
% sc    scaling factor
% rmsd  minimized root mean square deviation
%
% G. Jeschke, 2010

sc=sum(x.*y)/sum(y.^2);
diff=sc*y-x;
rmsd=sqrt(sum(diff.^2)/(length(diff)-1));


% --- Executes on button press in pushbutton_detach_fluctuations.
function pushbutton_detach_fluctuations_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_fluctuations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

snum=model.current_structure;
resnums=model.coarse(snum).resnums;
figure(2);
clf;

fluc2=model.GNM(model.current_structure).fluctuations;
[m,n]=size(fluc2);
for k=1:m,
    fluc2(k,k)=0;
end;
ma=max(max(fluc2));
for k=1:m,
    fluc2(k,k)=ma;
end;
fluc2=255*fluc2/max(max(fluc2));

image(fluc2);
set(gca,'YDir','normal');

axis equal
axis tight
colormap(handles.greyscale.^handles.saturation);
hold on;

state=get(handles.uitoggletool_chain,'State');
if strcmpi(state,'on'),
    chains=model.coarse(model.current_structure).chains;
    point3=get(handles.axes_fluctuations,'Position');
    dx=-2*(handles.xe-handles.xa)/point3(3);
    dy=0.5*(handles.ye-handles.ya)/point3(4);
    [mc,nc]=size(chains);
    for k=1:mc,
        xx=[chains(k,2),chains(k,2),chains(k,3),chains(k,3),chains(k,2)];
        yy=[chains(k,2),chains(k,3),chains(k,3),chains(k,2),chains(k,2)];
        plot(xx,yy,'m');
        [stag,ctag]=mk_address_parts([model.current_structure,chains(k,1)]);
        text(chains(k,3)+dx,chains(k,2)+dy,ctag,'color','m');
    end;
end;


% --- Executes on button press in pushbutton_detach_plot.
function pushbutton_detach_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

contents = get(handles.popupmenu_plot,'String');
type=contents{get(handles.popupmenu_plot,'Value')};

figure(2);
clf;

switch type
    case 'Mode distribution'
        plot(model.GNM(model.current_structure).lambda);
        xlabel('Mode number');
        ylabel('Normalized square frequency');
    case 'B factors'
        plot(model.coarse(model.current_structure).Bfactors,'k');
        hold on;
        msf=lin_fit(8*pi^2*model.GNM(model.current_structure).msf,model.coarse(model.current_structure).Bfactors);
        plot(msf,'r');
        xlabel('Residue number');
        ylabel('B factor');
    case 'Single mode'
        mode=model.GNM(model.current_structure).u(:,handles.mode+1);
        plot(mode,'k');
        hold on;
        plot([1,length(mode)],[0,0],'g');
        xlabel('Residue number');
        ylabel('Mode eigenvector');
end;


% --------------------------------------------------------------------
function uitoggletool_zoom_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state=get(hObject,'State');
if strcmpi(state,'on'),
    set(handles.uitoggletool_select,'State','off');
    set(handles.uitoggletool_inquire,'State','off');
    handles.crsr_mode=1;
end;
guidata(hObject,handles);

% --------------------------------------------------------------------
function uitoggletool_select_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state=get(hObject,'State');
if strcmpi(state,'on'),
    set(handles.uitoggletool_zoom,'State','off');
    set(handles.uitoggletool_inquire,'State','off');
    handles.crsr_mode=0;
end;
set(gcf,'Pointer','arrow');
guidata(hObject,handles);

% --------------------------------------------------------------------
function uitoggletool_inquire_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_inquire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state=get(hObject,'State');
if strcmpi(state,'on'),
    set(handles.uitoggletool_select,'State','off');
    set(handles.uitoggletool_zoom,'State','off');
    handles.crsr_mode=2;
end;

set(gcf,'Pointer','crosshair');
guidata(hObject,handles);

% --------------------------------------------------------------------
function uipushtool_full_scale_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_full_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

handles.xa=1;
handles.ya=1;
handles.xe=model.GNM(model.current_structure).residues;
handles.ye=model.GNM(model.current_structure).residues;
handles=update_correlation_plot(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function uipushtool_unselect_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain=cmd(hMain,'unselect *');
handles=update_correlation_plot(handles);
guidata(hObject,handles);


function handles=update_correlation_plot(handles)

global model

fluc2=model.GNM(model.current_structure).fluctuations;
[m,n]=size(fluc2);
for k=1:m,
    fluc2(k,k)=0;
end;
ma=max(max(fluc2));
for k=1:m,
    fluc2(k,k)=ma;
end;
fluc2=255*fluc2/max(max(fluc2));

if get(handles.radiobutton_neg_corr,'Value'),
    fluc2=-fluc2;
    colmap=handles.greyscale.^handles.saturation;
elseif get(handles.radiobutton_both_corr,'Value'),
    fluc2=(fluc2+100)/2;
    colmap=handles.tempscale;
    saturate=abs(linspace(0,1,50)).^handles.saturation;
    colmap(1:50,1)=colmap(1:50,1).*saturate';
    colmap(1:50,2)=colmap(1:50,2).*saturate';
    colmap(52:101,2)=colmap(52:101,2).*flipud(saturate');
    colmap(52:101,3)=colmap(52:101,3).*flipud(saturate');
else
    set(handles.radiobutton_pos_corr,'Value',1);
    colmap=handles.greyscale.^handles.saturation;
end;

snum=model.current_structure;
resnums=model.coarse(snum).resnums;

axes(handles.axes_fluctuations);
handles.corr_plot=image(fluc2);
set(gca,'YDir','normal');

axis equal
axis tight
colormap(colmap);
hold on;

state=get(handles.uitoggletool_chain,'State');
if strcmpi(state,'on'),
    chains=model.coarse(model.current_structure).chains;
    point3=get(handles.axes_fluctuations,'Position');
    dx=-2*(handles.xe-handles.xa)/point3(3);
    dy=0.5*(handles.ye-handles.ya)/point3(4);
    [mc,nc]=size(chains);
    for k=1:mc,
        xx=[chains(k,2),chains(k,2),chains(k,3),chains(k,3),chains(k,2)];
        yy=[chains(k,2),chains(k,3),chains(k,3),chains(k,2),chains(k,2)];
        plot(xx,yy,'m');
        [stag,ctag]=mk_address_parts([model.current_structure,chains(k,1)]);
        text(chains(k,3)+dx,chains(k,2)+dy,ctag,'color','m');
    end;
end;

axis([handles.xa,handles.xe,handles.ya,handles.ye]);
set(handles.corr_plot,'ButtonDownFcn',@corr_click);


function my_button_up(hObject, eventdata)

handles=guidata(hObject);
if handles.crsr_mode~=1,
    return
end;

old_units=get(0,'Units');
set(0,'Units','pixels');
point=get(0,'PointerLocation');
set(0,'Units',old_units);

old_units=get(gcf,'Units');
set(gcf,'Units','pixels');
point2=get(gcf,'Position');
set(gcf,'Units',old_units);

old_units=get(handles.axes_fluctuations,'Units');
set(handles.axes_fluctuations,'Units','pixels');
point3=get(handles.axes_fluctuations,'Position');
set(handles.axes_fluctuations,'Units',old_units);

whereami=point-point2(1:2)-point3(1:2);
xe=round(handles.xa+(handles.xe-handles.xa)*whereami(1)/point3(3));
ye=round(handles.ya+(handles.ye-handles.ya)*whereami(2)/point3(4));
xa=handles.xc;
ya=handles.yc;
handles.xa=min(xa,xe);
handles.xe=max(xa,xe);
handles.ya=min(ya,ye);
handles.ye=max(ya,ye);
handles=update_correlation_plot(handles);
guidata(hObject,handles);

function my_motion_fct(hObject, eventdata)

handles=guidata(hObject);
if handles.crsr_mode~=1 || ~isfield(handles,'xc') || ~isfield(handles,'yc'),
    return;
end;

cp = get(handles.axes_fluctuations,'CurrentPoint');
x=cp(1,1);
y=cp(1,2);
xdata=[handles.xc,handles.xc,x,x,handles.xc];
ydata=[handles.yc,y,y,handles.yc,handles.yc];

if isfield(handles,'zwindow')
    set(handles.zwindow,'XData',xdata,'YData',ydata);drawnow
end;


% --------------------------------------------------------------------
function uipushtool_3D_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

figure(1); clf;
fluc2=model.GNM(model.current_structure).fluctuations;
[m,n]=size(fluc2);
for k=1:m,
    fluc2(k,k)=0;
end;
ma=max(max(fluc2));
for k=1:m,
    fluc2(k,k)=ma;
end;

surf(fluc2);
colormap(handles.greyscale.^handles.saturation);
shading flat


% --------------------------------------------------------------------
function uitoggletool_chain_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_chain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=update_correlation_plot(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_minus.
function pushbutton_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.mode>1,
    handles.mode=handles.mode-1;
end;
set(handles.edit_mode_number,'String',sprintf('%5i',handles.mode));
handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_plus.
function pushbutton_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

if handles.mode<model.GNM(model.current_structure).residues-1,
    handles.mode=handles.mode+1;
end;
set(handles.edit_mode_number,'String',sprintf('%5i',handles.mode));
handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);


function edit_mode_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mode_number as text
%        str2double(get(hObject,'String')) returns contents of edit_mode_number as a double

global model

[v,handles]=edit_update_MMM(handles,hObject,1,model.GNM(model.current_structure).residues-1,handles.mode,'%5i',1);
handles.mode=v;
handles=update_auxiliary_plot(hObject,handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_mode_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mode_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [y1,p,rmsd,r]=lin_fit(x,y)

[y0,poi]=sort(y,'ascend');
x0=x(poi);
p=polyfit(x0,y0,1);
y1=polyval(p,x);
rmsd=y1-y;
rmsd=sqrt(sum(rmsd.*rmsd)/(length(rmsd)-1));
r0=corrcoef(x,y);
r=r0(1,2);

function linked_cursor

global model

stype=get(gcf,'SelectionType');
handles=guidata(gcf);

switch lower(stype)
    case 'normal'
        mark='b:';
        right=0;
        if ~isempty(handles.crsr_left),
            for k=1:length(handles.crsr_left),
                if ishandle(handles.crsr_left(k)),
                    delete(handles.crsr_left(k));
                end;
            end;
        end;
    case 'alt'
        mark='g:';
        right=1;
        if ~isempty(handles.crsr_right),
            for k=1:length(handles.crsr_right),
                if ishandle(handles.crsr_right(k)),
                    delete(handles.crsr_right(k));
                end;
            end;
        end;
end;
cp = get(handles.axes_plot,'CurrentPoint');
axes(handles.axes_plot);
hold on;
x=round(cp(1,1));
h1=plot([x,x],[handles.ya_aux,handles.ye_aux],mark);
axes(handles.axes_fluctuations);
h2=plot([x,x],[1,model.GNM(model.current_structure).residues],mark);
h3=plot([1,model.GNM(model.current_structure).residues],[x,x],mark);
if right,
    handles.crsr_right=[h1 h2 h3];
    handles.range(2)=x;
else
    handles.crsr_left=[h1 h2 h3];
    handles.range(1)=x;
end;

axes(handles.axes_fluctuations);
plot([x,x],[handles.ya_aux,handles.ye_aux],mark);

ind=model.coarse(model.current_structure).indices(x,:);
adr=mk_address(ind);
p=strfind(adr,']');
adr=adr(p+1:end);
set(handles.text_selected_2,'String',adr);

guidata(handles.axes_plot,handles);

% --- Executes on button press in pushbutton_clear_cursors.
function pushbutton_clear_cursors_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear_cursors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=update_auxiliary_plot(hObject,handles);
handles=update_correlation_plot(handles);
handles.range=[];
guidata(hObject,handles);


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

if ~isempty(handles.range) && length(handles.range)>=2,
    ind=model.coarse(model.current_structure).indices(handles.range(1),:);
    adr=mk_address(ind);
    cid=ind(2);
    for k=handles.range(1)+1:handles.range(2),
        ind=model.coarse(model.current_structure).indices(k,:);
        if ind(2)~=cid,
            [stag,ctag,modelnum,resnum]=mk_address_parts(model.coarse(model.current_structure).indices(k-1,:));
            adr=sprintf('%s-%i',adr,resnum);
            hMain=cmd(hMain,sprintf('select %s',adr));
            cid=ind(2);
            adr=mk_address(ind);
        end;
    end;
    [stag,ctag,modelnum,resnum]=mk_address_parts(ind);
    adr=sprintf('%s-%i',adr,resnum);
    hMain=cmd(hMain,sprintf('select %s',adr));
end;


% --------------------------------------------------------------------
function uipushtool_grey_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_grey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain=cmd(hMain,'color ! grey');
 
% --- Executes on button press in pushbutton_optimize_range.
function pushbutton_optimize_range_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_optimize_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model

contents = get(handles.popupmenu_plot,'String');
type=contents{get(handles.popupmenu_plot,'Value')};
if ~strcmpi(type,'Single mode'),
    add_msg_board('Range can be optimized only in single model display');
    return
end;
if isempty(handles.range) || length(handles.range)<2,
    add_msg_board('Range is not (fully) defined. Both cursors must be displayed.');
    return
end;
resa=handles.range(1);
rese=handles.range(2);
if resa>rese, % correct sequence, if required
    resa=handles.range(2);
    rese=handles.range(1);
end;
mode=model.GNM(model.current_structure).u(:,handles.mode+1);
mi=min(mode(resa:rese));
ma=max(mode(resa:rese));
if abs(mi)>ma,
    sgn=-1;
    amp=abs(mi);
else
    sgn=1;
    amp=ma;
end;
if sgn==1,
    if mode(resa)<amp/10,
        while mode(resa)<amp/10 && resa<length(mode),
            resa=resa+1;
        end;
    else
        while mode(resa)>=amp/10 && resa>1,
            resa=resa-1;
        end;
        resa=resa+1;
    end;
    if mode(rese)<amp/10,
        while mode(rese)<amp/10 && rese>1,
            rese=rese-1;
        end;
    else
        while mode(rese)>=amp/10 && rese<length(mode),
            rese=rese+1;
        end;
        rese=rese-1;
    end;
else
    if mode(resa)>-amp/10,
        while mode(resa)>-amp/10 && resa<length(mode),
            resa=resa+1;
        end;
    else
        while mode(resa)<=-amp/10 && resa>1,
            resa=resa-1;
        end;
        resa=resa+1;
    end;
    if mode(rese)>-amp/10,
        while mode(rese)>-amp/10 && rese>1,
            rese=rese-1;
        end;
    else
        while mode(rese)<=-amp/10 && rese<length(mode),
            rese=rese+1;
        end;
        rese=rese-1;
    end;    
end;
handles.range(1)=resa;
handles.range(2)=rese;

if ~isempty(handles.crsr_left),
    for k=1:length(handles.crsr_left),
        if ishandle(handles.crsr_left(k)),
            delete(handles.crsr_left(k));
        end;
    end;
end;
if ~isempty(handles.crsr_right),
    for k=1:length(handles.crsr_right),
        if ishandle(handles.crsr_right(k)),
            delete(handles.crsr_right(k));
        end;
    end;
end;

axes(handles.axes_plot);
hold on;
h1=plot([resa,resa],[handles.ya_aux,handles.ye_aux],'b:');
h1r=plot([rese,rese],[handles.ya_aux,handles.ye_aux],'g:');
axes(handles.axes_fluctuations);
h2=plot([resa,resa],[1,model.GNM(model.current_structure).residues],'b:');
h3=plot([1,model.GNM(model.current_structure).residues],[resa,resa],'b:');
h2r=plot([rese,rese],[1,model.GNM(model.current_structure).residues],'g:');
h3r=plot([1,model.GNM(model.current_structure).residues],[rese,rese],'g:');
handles.crsr_left=[h1 h2 h3];
handles.crsr_right=[h1r h2r h3r];
guidata(hObject,handles),

% --- Executes on button press in pushbutton_color_range.
function pushbutton_color_range_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_color_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

color_selection;
if isempty(hMain.color_selection),
    add_msg_board('Color selection cancelled by user.');
    return
else
    if isfloat(hMain.color_selection),
        com='color';
        arg=sprintf('%6.3f%6.3f%6.3f',hMain.color_selection);
    else
        com='colorscheme';
        arg=sprintf('%s',hMain.color_selection);
    end;
end;

if ~isempty(handles.range) && length(handles.range)>=2,
    ind=model.coarse(model.current_structure).indices(handles.range(1),:);
    adr=mk_address(ind);
    cid=ind(2);
    for k=handles.range(1)+1:handles.range(2),
        ind=model.coarse(model.current_structure).indices(k,:);
        if ind(2)~=cid,
            [stag,ctag,modelnum,resnum]=mk_address_parts(model.coarse(model.current_structure).indices(k-1,:));
            adr=sprintf('%s-%i',adr,resnum);
            hMain=cmd(hMain,sprintf('%s %s %s',com,adr,arg));
            cid=ind(2);
            adr=mk_address(ind);
        end;
    end;
    [stag,ctag,modelnum,resnum]=mk_address_parts(ind);
    adr=sprintf('%s-%i',adr,resnum);
    hMain=cmd(hMain,sprintf('%s %s %s',com,adr,arg));
end;


% --- Executes on button press in radiobutton_pos_corr.
function radiobutton_pos_corr_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_pos_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_pos_corr

if get(hObject,'Value'),
    set(handles.radiobutton_neg_corr,'Value',0);
    set(handles.radiobutton_both_corr,'Value',0);
end;
handles=update_correlation_plot(handles);
guidata(hObject,handles);

% --- Executes on button press in radiobutton_neg_corr.
function radiobutton_neg_corr_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_neg_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_neg_corr

if get(hObject,'Value'),
    set(handles.radiobutton_pos_corr,'Value',0);
    set(handles.radiobutton_both_corr,'Value',0);
end;
handles=update_correlation_plot(handles);
guidata(hObject,handles);


% --- Executes on button press in radiobutton_both_corr.
function radiobutton_both_corr_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_both_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_both_corr

if get(hObject,'Value'),
    set(handles.radiobutton_neg_corr,'Value',0);
    set(handles.radiobutton_pos_corr,'Value',0);
end;
handles=update_correlation_plot(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general

my_path=pwd;
cd(general.reports);

snum=model.current_structure;

Ca_coor=model.coarse(snum).Ca_coor;
rindices=model.coarse(snum).indices;
masses=model.coarse(snum).masses;
Bfactors=model.coarse(snum).Bfactors;
restypes=model.coarse(snum).restypes;
covar=model.GNM(snum).fluctuations;
[m,n]=size(rindices);
resnums=zeros(1,m);
for k=1:m,
    [stag,ctag,modelnum,resnum]=mk_address_parts(rindices(k,:));
    resnums(k)=resnum;
end;

[filename, pathname] = uiputfile('*.mat', 'Save coarse-grained model and covariance matrix');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Saving of model cancelled by user');
else
    reset_user_paths(pname);
    general.reports=pname;
    fname=fullfile(pathname, filename);
    msg=sprintf('Coarse-grained model saved as .mat file: %s',fname);
    add_msg_board(msg);
end
save(fname,'Ca_coor','rindices','masses','Bfactors','restypes','covar','resnums');

cd(my_path);
