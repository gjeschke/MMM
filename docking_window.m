function varargout = docking_window(varargin)
%DOCKING_WINDOW M-file for docking_window.fig
%      DOCKING_WINDOW, by itself, creates a new DOCKING_WINDOW or raises the existing
%      singleton*.
%
%      H = DOCKING_WINDOW returns the handle to a new DOCKING_WINDOW or the handle to
%      the existing singleton*.
%
%      DOCKING_WINDOW('Property','Value',...) creates a new DOCKING_WINDOW using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to docking_window_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DOCKING_WINDOW('CALLBACK') and DOCKING_WINDOW('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DOCKING_WINDOW.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help docking_window

% Last Modified by GUIDE v2.5 09-Sep-2013 16:07:18


% Y.Polyhach 09/2013

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @docking_window_OpeningFcn, ...
                   'gui_OutputFcn',  @docking_window_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before docking_window is made visible.
function docking_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)


global hMain
% global MMM_icon

% Old version with MMM figure icon, blocked because of warning
% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help_docking,'CData',cdata);

% global MMM_icon % check if needed!!!!!!!

% j = get(hObject,'javaframe');    
% j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  % create a java image and set the figure icon

% Choose default command line output for docking_window
% handles.output = hObject;
hMain.docking=[];

set(handles.save_stats,'Value',1.0);
set(handles.save_pdb,'Value',1.0);
set(handles.gridMethod,'Enable','off');
set(handles.fitMethod,'Enable','off');
set(handles.checkhomooligomer,'Value',0.0); % default behaviour: het dimer is more general
set(handles.checkhomooligomer,'Enable','off');
set(handles.valMultiplicity,'String',num2str(2));
set(handles.valMultiplicity,'String',' ');
set(handles.valMultiplicity,'Enable','off');
set(handles.LargeGridcheckbox,'Value',0.0); % default behaviour: a grid is assumed to be not too large at start
set(handles.finishexit,'Enable','off');
set(handles.rundocking,'Enable','off');
% set(handles.fitOnlyPanel,'Visible','off');
% set(handles.GridSearchPanel,'Visible','off');

% default values for the hetero grid search
def.al_grid_st=0;
def.al_grid_end=360;
def.al_grid_elem=9;
def.be_grid_st=0;
def.be_grid_end=180;
def.be_grid_end_homo=90;
def.be_grid_elem=8;
def.ga_grid_st=0;
def.ga_grid_end=360;
def.ga_grid_elem=9;

def.x_grid_st=-75;
def.x_grid_end=75;
def.x_grid_elem=9;
def.y_grid_st=-75;
def.y_grid_end=75;
def.y_grid_elem=9;
def.z_grid_st=-75;
def.z_grid_end=75;
def.z_grid_elem=9;

pstr=sprintf('%d',def.al_grid_st);
set(handles.get_Alpha_grid_st,'String',pstr);
pstr=sprintf('%d',def.al_grid_end);
set(handles.get_Alpha_grid_end,'String',pstr);
pstr=sprintf('%d',def.al_grid_elem);
set(handles.get_Alpha_grid_elem,'String',pstr);
pstr=sprintf('%d',def.be_grid_st);
set(handles.get_Beta_grid_st,'String',pstr);
pstr=sprintf('%d',def.be_grid_end);
set(handles.get_Beta_grid_end,'String',pstr);
pstr=sprintf('%d',def.be_grid_elem);
set(handles.get_Beta_grid_elem,'String',pstr);
pstr=sprintf('%d',def.ga_grid_st);
set(handles.get_Gama_grid_st,'String',pstr);
pstr=sprintf('%d',def.ga_grid_end);
set(handles.get_Gama_grid_end,'String',pstr);
pstr=sprintf('%d',def.ga_grid_elem);
set(handles.get_Gama_grid_elem,'String',pstr);

pstr=sprintf('%d',def.x_grid_st);
set(handles.get_x_grid_st,'String',pstr);
pstr=sprintf('%d',def.x_grid_end);
set(handles.get_x_grid_end,'String',pstr);
pstr=sprintf('%d',def.x_grid_elem);
set(handles.get_x_grid_elem,'String',pstr);
pstr=sprintf('%d',def.y_grid_st);
set(handles.get_y_grid_st,'String',pstr);
pstr=sprintf('%d',def.y_grid_end);
set(handles.get_y_grid_end,'String',pstr);
pstr=sprintf('%d',def.y_grid_elem);
set(handles.get_y_grid_elem,'String',pstr);
pstr=sprintf('%d',def.z_grid_st);
set(handles.get_z_grid_st,'String',pstr);
pstr=sprintf('%d',def.z_grid_end);
set(handles.get_z_grid_end,'String',pstr);
pstr=sprintf('%d',def.z_grid_elem);
set(handles.get_z_grid_elem,'String',pstr);

handles.grid.grSize=def.al_grid_elem*def.be_grid_elem*def.ga_grid_elem*def.x_grid_elem*def.y_grid_elem*def.z_grid_elem;
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);

if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
%     set(handles.LargeGridcheckbox,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
%     set(handles.LargeGridcheckbox,'Enable','on');
end

% stack value for large grids (same for homo and hetero cases)
def.stackSize=100000;
pstr=sprintf('%d',def.stackSize);
set(handles.get_stackSize,'String',pstr);
set(handles.get_stackSize,'Enable','off');

% default starting values for fitting when the window is initilized
% (after grid search is completed, the best grid point has to be taken as a starting point for fitting)
def.al_fit=0;
pstr=sprintf('%d',def.al_fit);
set(handles.get_Alpha_fitSt,'String',pstr);
def.be_fit=0;
pstr=sprintf('%d',def.be_fit);
set(handles.get_Beta_fitSt,'String',pstr);
def.ga_fit=0;
pstr=sprintf('%d',def.ga_fit);
set(handles.get_Gama_fitSt,'String',pstr);
def.x_fit=0;
pstr=sprintf('%d',def.x_fit);
set(handles.get_x_fitSt,'String',pstr);
def.y_fit=0;
pstr=sprintf('%d',def.y_fit);
set(handles.get_y_fitSt,'String',pstr);
def.z_fit=0;
pstr=sprintf('%d',def.z_fit);
set(handles.get_z_fitSt,'String',pstr);
% handles.dockFit=dockFit;
handles=appendFields(handles,def,1);
handles.def=def;

set(handles.bestfit_alpha,'Enable','off');
set(handles.bestfit_beta,'Enable','off');
set(handles.bestfit_gama,'Enable','off');
set(handles.bestfit_x,'Enable','off');
set(handles.bestfit_y,'Enable','off');
set(handles.bestfit_z,'Enable','off');

% Update handles structure
guidata(hObject, handles);
% guidata(hObject)
% uiwait(handles.docking_window);
% UIWAIT makes docking_window wait for user response (see UIRESUME)
% uiwait(handles.docking_window);


% --- Outputs from this function are returned to the command line.
function varargout = docking_window_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in load_pushbotton.
function load_pushbotton_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbotton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global model
global general % needed to enabling current path resetting (check if really needed!!!!!)

mypath=pwd;

cd(general.restraint_files);
% snum=model.current_structure;
% axes(handles.axes_model);
[fname,pname,findex]=uigetfile('*.dat','Load constraints from file');
% if isequal(fname,0) || isequal(pname,0) % to be used when no findex is returned
if isequal(findex,0)
    add_msg_board('Constraint loading cancelled by user');
else
%     reset_user_paths(pname); % !!!!!! check if convenient!!!!!
    general.restraint_files=pname; % check how this would work!!!!!
    restraints=rd_restraints_docking(fullfile(pname,fname));
    handles.ensemble=restraints.ensemble;
    handles.uncertainty=restraints.uncertainty;
    handles.exclude=restraints.exclude;
    handles.target=restraints.target;
    if isfield(restraints,'basis'),
        handles.basis=restraints.basis;
    end;
%     if isfield(restraints,'tif'), % auxilary functionality - think how to have it as well!!!!!!
%         handles.tif=restraints.tif;
%     else
%         handles.tif=ENM_param.tif;
%     end;
%     if isfield(restraints,'mmax'),
%         handles.mmax=restraints.mmax;
%     else
%         handles.mmax=ENM_param.mmax;
%     end;
    hfig=gcf;
    set(hfig,'Pointer','watch');
    drawnow;
    [DEER,cancelled]=process_DEER_restraints(restraints);
    set(hfig,'Pointer','arrow');
    if cancelled,
        add_msg_board('Processing of DEER constraints cancelled.');
        return
    end;
%     keyboard
    handles.dockBasis=DEER(1).dockBasis;
    handles.dockBasis.constrFile=fullfile(pname,fname);
    handles.dockBasis.PDB=restraints.PDB;
    handles.dockBasis.restraintsFile=fullfile(pname,fname);
    DEER=rmfield(DEER,'dockBasis');
    handles.DEER=DEER;

    set(handles.gridMethod,'Value',0.0);
    set(handles.fitMethod,'Value',1.0);
    set(handles.gridMethod,'Enable','on');
    set(handles.fitMethod,'Enable','on');
    set(handles.fitStartValues,'Visible','on');
    set(handles.fitBestValues,'Visible','on');
    set(handles.save_stats,'Visible','on');
    set(handles.save_pdb,'Visible','on');
    set(handles.rundocking,'Enable','on');
%     set(handles.append2model,'Visible','on');
%     set(handles.save_mat_file,'Visible','on');
    switch handles.dockBasis.oligomer
        case 1
            % % !!!!
            if isfield(handles.dockBasis,'multiplicity')
                multiplicity=handles.dockBasis.multiplicity;                
            else
                multiplicity=2;
            end
            Gama_oligomer=360/multiplicity;
      
            set(handles.checkhomooligomer,'Value',1.0);
            set(handles.checkhomooligomer,'Enable','on');
            set(handles.get_Beta_grid_end,'String',num2str(90));
            set(handles.valMultiplicity,'Enable','on');
            set(handles.valMultiplicity,'String',num2str(multiplicity));
            set(handles.get_z_grid_st,'Enable','off');
            set(handles.get_z_grid_st,'String',num2str(0));
            set(handles.get_z_grid_end,'Enable','off');
            set(handles.get_z_grid_end,'String',num2str(0));
            set(handles.get_z_grid_elem,'Enable','off');
            set(handles.get_z_grid_elem,'String',num2str(1));
            set(handles.get_Gama_grid_st,'Enable','off');
            set(handles.get_Gama_grid_st,'String',num2str(0));
            set(handles.get_Gama_grid_end,'Enable','off');
            set(handles.get_Gama_grid_end,'String',num2str(Gama_oligomer));
            set(handles.get_Gama_grid_elem,'Enable','off');
            set(handles.get_Gama_grid_elem,'String',num2str(1));
            set(handles.get_Gama_fitSt,'Enable','off');
            set(handles.get_Gama_fitSt,'String',num2str(Gama_oligomer));
            set(handles.get_z_fitSt,'Enable','off');
                                   
            % display here may be a warning that by default a homodimer
            % case is treated!!!
        case 2
            set(handles.checkhomooligomer,'Value',0.0);
            set(handles.checkhomooligomer,'Enable','off');
%             set(handles.valMultiplicity,'String',num2str(1));
            set(handles.valMultiplicity,'String',' ');
            set(handles.valMultiplicity,'Enable','off');
    end
    
    
%     [direct,cancelled]=process_direct_restraints(restraints);
%     if cancelled,
%         add_msg_board('Processing of direct constraints cancelled.');
%         return
%     end;
%     handles.direct=direct;
    
%     [displacements,cancelled]=process_displacement_restraints(restraints); % check if needed!!!!
%     if cancelled,
%         add_msg_board('Processing of displacement constraints cancelled.');
%         return
%     end;
%     handles.displacements=displacements;

%     set(handles.checkbox_restraints,'Enable','on'); % flag indicating that restraints have been loaded - to be implemented!!!!!
%     set(handles.checkbox_restraints,'Value',1);
%     

%     if ~isempty(DEER),
%     dvec=zeros(length(handles.DEER),3);
%     for k=1:length(handles.DEER),
% %         handles.DEER(k).l1=line(handles.DEER(k).xyz1(1),handles.DEER(k).xyz1(2),handles.DEER(k).xyz1(3),...
% %             'Parent',handles.restraint_graphics,'Marker','.','Color','b');
% %         handles.DEER(k).l2=line(handles.DEER(k).xyz2(1),handles.DEER(k).xyz2(2),handles.DEER(k).xyz2(3),...
% %             'Parent',handles.restraint_graphics,'Marker','.','Color','b');
%         x=[handles.DEER(k).xyz1(1) handles.DEER(k).xyz2(1)];
%         y=[handles.DEER(k).xyz1(2) handles.DEER(k).xyz2(2)];
%         z=[handles.DEER(k).xyz1(3) handles.DEER(k).xyz2(3)];
%         r0=norm(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/10;
%         dvec(k,:)=(handles.DEER(k).r-r0)*(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/(10*r0);
%         det=abs(r0-handles.DEER(k).r)/handles.DEER(k).sigr;
%         if det>2,
%             col='r';
%         elseif det>1,
%             col=[255,190,0]/255; % yellow-orange
%         else
%             col='g';
%         end;
%         handles.DEER(k).ll=line(x,y,z,'Parent',handles.restraint_graphics,'Color',col,'LineWidth',1.5,'LineStyle',':');
%         cindices=handles.DEER(k).indices;
%         f1=false;
%         f2=false;
%         for l=1:length(model.coarse(snum).indices),
%             diff=cindices(1,:)-model.coarse(snum).indices(l,:);
%             if sum(abs(diff))==0,
%                 handles.DEER(k).res1=l;
%                 f1=true;
%                 x=[handles.DEER(k).xyz1(1) model.coarse(snum).Ca_coor(l,1)];
%                 y=[handles.DEER(k).xyz1(2) model.coarse(snum).Ca_coor(l,2)];
%                 z=[handles.DEER(k).xyz1(3) model.coarse(snum).Ca_coor(l,3)];
%                 handles.DEER(k).rl1=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
%             end;
%             diff=cindices(2,:)-model.coarse(snum).indices(l,:);
%             if sum(abs(diff))==0,
%                 handles.DEER(k).res2=l;
%                 f2=true;
%                 x=[handles.DEER(k).xyz2(1) model.coarse(snum).Ca_coor(l,1)];
%                 y=[handles.DEER(k).xyz2(2) model.coarse(snum).Ca_coor(l,2)];
%                 z=[handles.DEER(k).xyz2(3) model.coarse(snum).Ca_coor(l,3)];
%                 handles.DEER(k).rl2=line(x,y,z,'Parent',handles.restraint_graphics,'Color','b','LineWidth',1);
%             end;
%         end;
%         if ~f1,
%             add_msg_board(sprintf('Warning: Residue for first label of restraint %i not in network model.',k));
%         end;
%         if ~f2,
%             add_msg_board(sprintf('Warning: Residue for second label of restraint %i not in network model.',k));
%         end;
%     end;
%     dvec=dvec/sqrt(sum(sum(dvec.^2)));
%     elseif ~isempty(handles.direct),
%         dvec=zeros(length(handles.direct),3);
%         network0=model.coarse(snum).Ca_coor;
%         for k=1:length(handles.direct),
%             xyz1=network0(handles.direct(k,1),:);
%             xyz2=network0(handles.direct(k,2),:);
%             line(xyz1(1),xyz1(2),xyz1(3),...
%                 'Parent',handles.restraint_graphics,'Marker','.','Color','b');
%             line(xyz2(1),xyz2(2),xyz2(3),...
%                 'Parent',handles.restraint_graphics,'Marker','.','Color','b');
%             x=[xyz1(1) xyz2(1)];
%             y=[xyz1(2) xyz2(2)];
%             z=[xyz1(3) xyz2(3)];
%             r0=norm(xyz1-xyz2)/10;
%             dvec(k,:)=(handles.direct(k,4)-r0)*(xyz1-xyz2)/(10*r0);
%             det=abs(r0-handles.direct(k,4))/handles.direct(k,5);
%             if det>2,
%                 col='r';
%             elseif det>1,
%                 col=[255,190,0]/255; % yellow-orange
%             else
%                 col='g';
%             end;
%             line(x,y,z,'Parent',handles.restraint_graphics,'Color',col,'LineWidth',1.5,'LineStyle',':');
%         end;
%     end;
end;

% if ~isempty(DEER),
%     overlaps=zeros(1,100);
%     for k=1:100,
%         evec=model.ANM(model.current_structure).u(:,k+6);
%         m=length(evec)/3;
%         mode=reshape(evec,3,m);
%         mdvec=zeros(length(handles.DEER),3);
%         for kk=1:length(handles.DEER),
%             dchange=mode(:,handles.DEER(kk).res1)-mode(:,handles.DEER(kk).res2);
%             mdvec(kk,:)=dchange';
%         end;
%         mdvec=mdvec/sqrt(sum(sum(mdvec.^2)));
%         overlaps(k)=abs(sum(sum(dvec.*mdvec)));
%     end;
%     % figure(1); clf;
%     % plot(overlaps,'k.');
%     % figure(2); clf;
%     % plot(cumsum(overlaps),'k');
%     axes(handles.axes_plot);
%     cla;
%     set(handles.text_auxiliary,'String','Constraint matching');
%     ma=0;
%     rmsd=0;
%     fom=0; % figure of merit
%     for k=1:length(handles.DEER),
%         r0=norm(handles.DEER(k).xyz1-handles.DEER(k).xyz2)/10;
%         rmsd=rmsd+(r0-handles.DEER(k).r)^2;
%         det=abs(r0-handles.DEER(k).r)/handles.DEER(k).sigr;
%         fom=fom+det^2;
%         if det>2,
%             col='r';
%         elseif det>1,
%             col=[255,190,0]/255; % yellow-orange
%         else
%             col='g';
%         end;
%         errorbar(k,handles.DEER(k).r,handles.DEER(k).sigr,'k');
%         if handles.DEER(k).r+handles.DEER(k).sigr>ma,
%             ma=handles.DEER(k).r+handles.DEER(k).sigr;
%         end;
%         if r0>ma,
%             ma=r0;
%         end;
%         hold on;
%         plot(k,r0,'.','Color',col);
%     end;
%     handles.rmsd=sqrt(rmsd/length(handles.DEER));
%     handles.fom=sqrt(fom/length(handles.DEER));
%     axis([0,length(handles.DEER)+1,0,1.05*ma]);
%     xlabel('Constraint number');
%     ylabel('Distance (nm)');
%     set(handles.text_info,'String',sprintf('Loaded %i DEER restraints with rmsd of %5.2f nm to template',length(handles.DEER),handles.rmsd));
%     set(handles.text_auxiliary_msg,'String',sprintf('Figure of merit: %6.4f',handles.fom));
% end;
% set(handles.pushbutton_fit,'Enable','on');
% set(handles.pushbutton_parametrize,'Enable','on');
% cd(my_path);
% keyboard
cd(mypath);
guidata(hObject,handles);

       
function [direct,cancelled]=process_direct_restraints(restraints)

global model

cancelled=false;
if ~isfield(restraints,'direct'),
    direct=[];
    return;
end;

snum=model.current_structure;

md=length(restraints.direct);
direct=zeros(md,5);

cindices=model.coarse(snum).indices;
[mn,nn]=size(cindices);
for k=1:md,
    adr1=restraints.direct(k).adr1;
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Specified residue %s does not exist in template structure.',adr1));
        cancelled=true;
        direct=[];
        return;
    end;
    net1=0;
    for kk=1:mn,
        match=sum(abs(ind1-cindices(kk,:)));
        if match==0,
           net1=kk;
           break;
        end;
    end;
    direct(k,1)=net1;
    adr2=restraints.direct(k).adr2;
    ind2=resolve_address(adr2);
    net2=0;
    for kk=1:mn,
        match=sum(abs(ind2-cindices(kk,:)));
        if match==0,
           net2=kk;
           break;
        end;
    end;
    direct(k,2)=net2;    
    direct(k,4)=restraints.direct(k).r;
    direct(k,5)=restraints.direct(k).sigr;
    if net1==0 || net2==0,
        add_msg_board('ERROR: Constrained C_alpha-C_alpha distance does not exist in network model.');
        cancelled=true;
        direct=[];
        return;
    end;
end;

function get_Alpha_fit_Callback(hObject, eventdata, handles)
% hObject    handle to get_Alpha_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Alpha_fit as text
%        str2double(get(hObject,'String')) returns contents of get_Alpha_fit as a double

[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.dockFit.al_fit=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_Alpha_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Alpha_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Beta_fit_Callback(hObject, eventdata, handles)
% hObject    handle to get_Beta_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Beta_fit as text
%        str2double(get(hObject,'String')) returns contents of get_Beta_fit as a double

[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.dockFit.be_fit=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Beta_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Beta_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function get_Gama_fit_Callback(hObject, eventdata, handles)
% hObject    handle to get_Gama_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Gama_fit as text
%        str2double(get(hObject,'String')) returns contents of get_Gama_fit as a double

[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.dockFit.ga_fit=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Gama_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Gama_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function get_x_fit_Callback(hObject, eventdata, handles)
% hObject    handle to get_x_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_x_fit as text
%        str2double(get(hObject,'String')) returns contents of get_x_fit as a double


[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.x_fit=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_x_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_x_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_y_fit_Callback(hObject, eventdata, handles)
% hObject    handle to get_y_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_y_fit as text
%        str2double(get(hObject,'String')) returns contents of get_y_fit as a double

[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.y_fit=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_y_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_y_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_z_fit_Callback(hObject, eventdata, handles)
% hObject    handle to get_z_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_z_fit as text
%        str2double(get(hObject,'String')) returns contents of get_z_fit as a double

[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.z_fit=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_z_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_z_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function show_gridSize_Callback(hObject, eventdata, handles)
% hObject    handle to show_gridSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of show_gridSize as text
%        str2double(get(hObject,'String')) returns contents of show_gridSize as a double



% --- Executes during object creation, after setting all properties.
function show_gridSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_gridSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function get_stackSize_Callback(hObject, eventdata, handles)
% hObject    handle to get_stackSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_stackSize as text
%        str2double(get(hObject,'String')) returns contents of get_stackSize as a double

[v,handles]=edit_update_MMM(handles,hObject,100000,1000000,100000,'%d',1);
handles.stackSize=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_stackSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_stackSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function get_z_grid_elem_Callback(hObject, eventdata, handles)
% hObject    handle to get_z_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_z_grid_elem as text
%        str2double(get(hObject,'String')) returns contents of get_z_grid_elem as a double

[v,handles]=edit_update_MMM(handles,hObject,1,1e7,9,'%d',1);
handles.z_grid_elem=v;
if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
    set(handles.get_stackSize,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
    set(handles.LargeGridcheckbox,'Enable','on');
    set(handles.get_stackSize,'Enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_z_grid_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_z_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_z_grid_end_Callback(hObject, eventdata, handles)
% hObject    handle to get_z_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_z_grid_end as text
%        str2double(get(hObject,'String')) returns contents of get_z_grid_end as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.z_grid_end=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_z_grid_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_z_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_z_grid_st_Callback(hObject, eventdata, handles)
% hObject    handle to get_z_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_z_grid_st as text
%        str2double(get(hObject,'String')) returns contents of get_z_grid_st as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);

handles.z_grid_st=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_z_grid_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_z_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_y_grid_end_Callback(hObject, eventdata, handles)
% hObject    handle to get_y_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_y_grid_end as text
%        str2double(get(hObject,'String')) returns contents of get_y_grid_end as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.y_grid_end=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_y_grid_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_y_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_y_grid_elem_Callback(hObject, eventdata, handles)
% hObject    handle to get_y_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_y_grid_elem as text
%        str2double(get(hObject,'String')) returns contents of get_y_grid_elem as a double
[v,handles]=edit_update_MMM(handles,hObject,1,1e7,9,'%d',1);

handles.y_grid_elem=v;
if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
    set(handles.get_stackSize,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
    set(handles.LargeGridcheckbox,'Enable','on');
    set(handles.get_stackSize,'Enable','on');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_y_grid_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_y_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_y_grid_st_Callback(hObject, eventdata, handles)
% hObject    handle to get_y_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_y_grid_st as text
%        str2double(get(hObject,'String')) returns contents of get_y_grid_st as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.y_grid_st=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_y_grid_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_y_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_x_grid_elem_Callback(hObject, eventdata, handles)
% hObject    handle to get_x_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_x_grid_elem as text
%        str2double(get(hObject,'String')) returns contents of get_x_grid_elem as a double
[v,handles]=edit_update_MMM(handles,hObject,1,1e7,9,'%d',1);

handles.x_grid_elem=v;
if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
    set(handles.get_stackSize,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
    set(handles.LargeGridcheckbox,'Enable','on');
    set(handles.get_stackSize,'Enable','on');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_x_grid_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_x_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_x_grid_st_Callback(hObject, eventdata, handles)
% hObject    handle to get_x_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_x_grid_st as text
%        str2double(get(hObject,'String')) returns contents of get_x_grid_st as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.x_grid_st=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_x_grid_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_x_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_x_grid_end_Callback(hObject, eventdata, handles)
% hObject    handle to get_x_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_x_grid_end as text
%        str2double(get(hObject,'String')) returns contents of get_x_grid_end as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.x_grid_end=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_x_grid_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_x_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Alpha_grid_st_Callback(hObject, eventdata, handles)
% hObject    handle to get_Alpha_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Alpha_grid_st as text
%        str2double(get(hObject,'String')) returns contents of get_Alpha_grid_st as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.al_grid_st=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Alpha_grid_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Alpha_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% pstr=sprintf('%d',handles.stackSize);
% set(hObject,'String',pstr);


function get_Alpha_grid_end_Callback(hObject, eventdata, handles)
% hObject    handle to get_Alpha_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Alpha_grid_end as text
%        str2double(get(hObject,'String')) returns contents of get_Alpha_grid_end as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.al_grid_end=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Alpha_grid_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Alpha_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Alpha_grid_elem_Callback(hObject, eventdata, handles)
% hObject    handle to get_Alpha_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Alpha_grid_elem as text
%        str2double(get(hObject,'String')) returns contents of get_Alpha_grid_elem as a double
[v,handles]=edit_update_MMM(handles,hObject,1,1e7,0,'%d',1);

handles.al_grid_elem=v;
if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
    set(handles.get_stackSize,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
    set(handles.LargeGridcheckbox,'Enable','on');
    set(handles.get_stackSize,'Enable','on');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Alpha_grid_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Alpha_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% handles.al_grid_elem=9;
% pstr=sprintf('%d',handles.al_grid_elem);
% set(hObject,'String',pstr);
% guidata(hObject,handles);



function get_Beta_grid_st_Callback(hObject, eventdata, handles)
% hObject    handle to get_Beta_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Beta_grid_st as text
%        str2double(get(hObject,'String')) returns contents of get_Beta_grid_st as a double
[v,handles]=edit_update_MMM(handles,hObject,0,180,0,'%0.2f',0);
handles.be_grid_st=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Beta_grid_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Beta_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Beta_grid_elem_Callback(hObject, eventdata, handles)
% hObject    handle to get_Beta_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Beta_grid_elem as text
%        str2double(get(hObject,'String')) returns contents of get_Beta_grid_elem as a double
[v,handles]=edit_update_MMM(handles,hObject,1,1e7,0,'%d',1);
handles.be_grid_elem=v;
if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
    set(handles.get_stackSize,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
    set(handles.LargeGridcheckbox,'Enable','on');
    set(handles.get_stackSize,'Enable','on');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Beta_grid_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Beta_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% handles.be_grid_elem=8;
% pstr=sprintf('%d',handles.be_grid_elem);
% set(hObject,'String',pstr);
% guidata(hObject,handles);



function get_Beta_grid_end_Callback(hObject, eventdata, handles)
% hObject    handle to get_Beta_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Beta_grid_end as text
%        str2double(get(hObject,'String')) returns contents of get_Beta_grid_end as a double

if get(handles.checkhomooligomer,'Value')
    [v,handles]=edit_update_MMM(handles,hObject,0,90,0,'%0.2f',0);
else
    [v,handles]=edit_update_MMM(handles,hObject,0,180,0,'%0.2f',0);
end
handles.be_grid_end=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Beta_grid_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Beta_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Gama_grid_st_Callback(hObject, eventdata, handles)
% hObject    handle to get_Gama_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Gama_grid_st as text
%        str2double(get(hObject,'String')) returns contents of get_Gama_grid_st as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.ga_grid_st=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Gama_grid_st_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Gama_grid_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function get_Gama_grid_end_Callback(hObject, eventdata, handles)
% hObject    handle to get_Gama_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Gama_grid_end as text
%        str2double(get(hObject,'String')) returns contents of get_Gama_grid_end as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.ga_grid_end=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Gama_grid_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Gama_grid_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Gama_grid_elem_Callback(hObject, eventdata, handles)
% hObject    handle to get_Gama_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Gama_grid_elem as text
%        str2double(get(hObject,'String')) returns contents of get_Gama_grid_elem as a double
[v,handles]=edit_update_MMM(handles,hObject,1,1e7,0,'%d',1);
handles.ga_grid_elem=v;
if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
    set(handles.get_stackSize,'Enable','off');
else
    set(handles.LargeGridcheckbox,'Value',1.0);
    set(handles.LargeGridcheckbox,'Enable','on');
    set(handles.get_stackSize,'Enable','on');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Gama_grid_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Gama_grid_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close docking_window.
function docking_window_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to docking_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

% global hMain
delete(hObject);


% --- Executes on button press in autoFitbox.
function autoFitbox_Callback(hObject, eventdata, handles)
% hObject    handle to autoFitbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoFitbox
guidata(hObject,handles);




% --- Executes on button press in LargeGridcheckbox.
function LargeGridcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to LargeGridcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LargeGridcheckbox
checktest=get(hObject,'Value');
if checktest==1
    set(handles.get_stackSize,'Enable','on');
    guidata(hObject,handles);
else
    set(handles.get_stackSize,'Enable','off');
    guidata(hObject,handles);
end



function bestfit_gama_Callback(hObject, eventdata, handles)
% hObject    handle to bestfit_gama (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bestfit_gama as text
%        str2double(get(hObject,'String')) returns contents of bestfit_gama as a double


% --- Executes during object creation, after setting all properties.
function bestfit_gama_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bestfit_gama (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bestfit_z_Callback(hObject, eventdata, handles)
% hObject    handle to bestfit_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bestfit_z as text
%        str2double(get(hObject,'String')) returns contents of bestfit_z as a double


% --- Executes during object creation, after setting all properties.
function bestfit_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bestfit_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bestfit_y_Callback(hObject, eventdata, handles)
% hObject    handle to bestfit_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bestfit_y as text
%        str2double(get(hObject,'String')) returns contents of bestfit_y as a double


% --- Executes during object creation, after setting all properties.
function bestfit_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bestfit_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bestfit_x_Callback(hObject, eventdata, handles)
% hObject    handle to bestfit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bestfit_x as text
%        str2double(get(hObject,'String')) returns contents of bestfit_x as a double


% --- Executes during object creation, after setting all properties.
function bestfit_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bestfit_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bestfit_beta_Callback(hObject, eventdata, handles)
% hObject    handle to bestfit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bestfit_beta as text
%        str2double(get(hObject,'String')) returns contents of bestfit_beta as a double


% --- Executes during object creation, after setting all properties.
function bestfit_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bestfit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bestfit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to bestfit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bestfit_alpha as text
%        str2double(get(hObject,'String')) returns contents of bestfit_alpha as a double


% --- Executes during object creation, after setting all properties.
function bestfit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bestfit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Alpha_fitSt_Callback(hObject, eventdata, handles)
% hObject    handle to get_Alpha_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Alpha_fitSt as text
%        str2double(get(hObject,'String')) returns contents of get_Alpha_fitSt as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.al_fit=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_Alpha_fitSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Alpha_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Beta_fitSt_Callback(hObject, eventdata, handles)
% hObject    handle to get_Beta_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Beta_fitSt as text
%        str2double(get(hObject,'String')) returns contents of get_Beta_fitSt as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.al_fit=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_Beta_fitSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Beta_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_Gama_fitSt_Callback(hObject, eventdata, handles)
% hObject    handle to get_Gama_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_Gama_fitSt as text
%        str2double(get(hObject,'String')) returns contents of get_Gama_fitSt as a double
[v,handles]=edit_update_MMM(handles,hObject,0,360,0,'%0.2f',0);
handles.ga_fit=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_Gama_fitSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_Gama_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_x_fitSt_Callback(hObject, eventdata, handles)
% hObject    handle to get_x_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_x_fitSt as text
%        str2double(get(hObject,'String')) returns contents of get_x_fitSt as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.x_fit=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_x_fitSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_x_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_y_fitSt_Callback(hObject, eventdata, handles)
% hObject    handle to get_y_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_y_fitSt as text
%        str2double(get(hObject,'String')) returns contents of get_y_fitSt as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.y_fit=v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function get_y_fitSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_y_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function get_z_fitSt_Callback(hObject, eventdata, handles)
% hObject    handle to get_z_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_z_fitSt as text
%        str2double(get(hObject,'String')) returns contents of get_z_fitSt as a double
[v,handles]=edit_update_MMM(handles,hObject,-75,75,0,'%0.2f',0);
handles.z_fit=v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function get_z_fitSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_z_fitSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkhomooligomer.
function checkhomooligomer_Callback(hObject, eventdata, handles)
% hObject    handle to checkhomooligomer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkhomooligomer

% keyboard
checkhom=get(hObject,'Value');
% hetBeta=str2double(get(handles.get_Beta_grid_end,'String'));
valMult=str2double(get(handles.valMultiplicity,'String'));
if checkhom==1
    set(handles.valMultiplicity,'Enable','on');
%             set(handles.valMultiplicity,'String',num2str(valMult));
    set(handles.get_z_grid_st,'Enable','off');
    set(handles.get_z_grid_st,'String',num2str(0));
    set(handles.get_z_grid_end,'Enable','off');
    set(handles.get_z_grid_end,'String',num2str(0));
    set(handles.get_z_grid_elem,'Enable','off');
    set(handles.get_z_grid_elem,'String',num2str(1));
    set(handles.get_Beta_grid_end,'String',num2str(90));
    set(handles.get_Gama_grid_st,'Enable','off');
    set(handles.get_Gama_grid_st,'String',num2str(0));
    set(handles.get_Gama_grid_end,'Enable','off');
    Gama_oligomer=360/valMult;
    set(handles.get_Gama_grid_end,'String',num2str(Gama_oligomer));
    set(handles.get_Gama_grid_elem,'Enable','off');
    set(handles.get_Gama_grid_elem,'String',num2str(1));
    set(handles.get_Gama_fitSt,'Enable','off');
    set(handles.get_Gama_fitSt,'String',num2str(Gama_oligomer));
    set(handles.get_z_fitSt,'Enable','off');
    set(handles.get_z_fitSt,'String',num2str(0));
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
    pstr=sprintf('%g',handles.grid.grSize);
    set(handles.show_gridSize,'String',pstr);

%             set(handles.bestfit_gama,'Enable','off');
%             set(handles.bestfit_gama,'String',num2str(Gama_oligomer));
%             set(handles.bestfit_z,'Enable','off');
%             set(handles.bestfit_z,'String',num2str(0));

else
    set(handles.valMultiplicity,'Enable','off');
%             set(handles.valMultiplicity,'String',' ');
    set(handles.get_z_grid_st,'Enable','on');
    set(handles.get_z_grid_st,'String',num2str(-75));
    set(handles.get_z_grid_end,'Enable','on');
    set(handles.get_z_grid_end,'String',num2str(75));
    set(handles.get_z_grid_elem,'Enable','on');
    set(handles.get_z_grid_elem,'String',num2str(9));
    set(handles.get_Beta_grid_end,'String',num2str(180));
    set(handles.get_Gama_grid_st,'Enable','on');
    set(handles.get_Gama_grid_st,'String',num2str(0));
    set(handles.get_Gama_grid_end,'Enable','on');
    set(handles.get_Gama_grid_end,'String',num2str(360));
    set(handles.get_Gama_grid_elem,'Enable','on');
    set(handles.get_Gama_grid_elem,'String',num2str(9))
    set(handles.get_Gama_fitSt,'Enable','on');
    set(handles.get_Gama_fitSt,'String',num2str(0));
    set(handles.get_z_fitSt,'Enable','on');
    set(handles.get_z_fitSt,'String',num2str(0));
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
    pstr=sprintf('%g',handles.grid.grSize);
    set(handles.show_gridSize,'String',pstr);
%             set(handles.bestfit_gama,'Enable','off');
%             set(handles.bestfit_gama,'String',num2str(0));
%             set(handles.bestfit_z,'Enable','off');
%             set(handles.bestfit_z,'String',num2str(0));

end
guidata(hObject,handles);


function valMultiplicity_Callback(hObject, eventdata, handles)
% hObject    handle to valMultiplicity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valMultiplicity as text
%        str2double(get(hObject,'String')) returns contents of valMultiplicity as a double
% keyboard
[v,handles]=edit_update_MMM(handles,hObject,2,100,2,'%d',1);
handles.homMulti=v;
set(handles.get_Gama_grid_end,'String',num2str(360/v));
set(handles.get_Gama_fitSt,'String',num2str(360/v));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function valMultiplicity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valMultiplicity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_stats.
function save_stats_Callback(hObject, eventdata, handles)
% hObject    handle to save_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_stats
guidata(hObject,handles);


% --- Executes on button press in gridMethod.
function gridMethod_Callback(hObject, eventdata, handles)
% hObject    handle to gridMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gridMethod

if get(handles.checkhomooligomer,'Value')
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.x_grid_elem*handles.y_grid_elem;
    set(handles.get_Beta_grid_end,'String',num2str(90));
else
    handles.grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
end
pstr=sprintf('%g',handles.grid.grSize);
set(handles.show_gridSize,'String',pstr);
if handles.grid.grSize<1000000
    set(handles.LargeGridcheckbox,'Value',0.0);
    set(handles.LargeGridcheckbox,'Enable','off');
end
set(handles.fitMethod,'Value',0.0);
set(handles.fitOnlyPanel,'Visible','off');
set(handles.EulerRotPanel,'Visible','on');
set(handles.TransPanel,'Visible','on');
set(handles.SearchModePanel,'Visible','on');
set(handles.show_gridSize,'Visible','on');
set(handles.gridSize_text,'Visible','on');
set(handles.autoFitbox,'Visible','on');
set(handles.save_stats,'Visible','on');

guidata(hObject,handles);

% --- Executes on button press in fitMethod.
function fitMethod_Callback(hObject, eventdata, handles)
% hObject    handle to fitMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fitMethod
% checkM=get(hObject,'Value');
% if checkM==1
set(handles.gridMethod,'Value',0.0);
set(handles.fitOnlyPanel,'Visible','on');
set(handles.EulerRotPanel,'Visible','off');
set(handles.TransPanel,'Visible','off');
set(handles.SearchModePanel,'Visible','off');
set(handles.show_gridSize,'Visible','off');
set(handles.gridSize_text,'Visible','off');
set(handles.autoFitbox,'Visible','off');
% set(handles.save_stats,'Visible','off');
%  else
%     
% end
guidata(hObject,handles);


% --- Executes on button press in rundocking.
function rundocking_Callback(hObject, eventdata, handles)
% hObject    handle to rundocking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hMain
global model
global general

hfig=gcf;
set(hfig,'Pointer','watch');
drawnow;

my_dir=pwd;
cd(general.restraint_files);

dockBasis=handles.dockBasis;
sites=dockBasis.sites;
constraints.pairs_exp=dockBasis.pairs_exp;

fitORdock=get(handles.fitMethod,'Value');
homAShom=get(handles.checkhomooligomer,'Value');
dockcase=num2str([dockBasis.oligomer homAShom]);
fitOptions.MaxIter=2000; % (by the way default of fminsearch is given by 200*number of variables)
fitOptions.MaxFunEvals=2000; % (by the way default of fminsearch is given by 200*number of variables)
% fitOptions.Display='iter'; % no showing every fit step by default
fitOptions.PlotFcns=@optimplotx;
dockStat0=[];
dockStat0.constraints=constraints;
min_grid=[];
out_v=[];
v0=[];
if fitORdock
    switch dockcase
        case '1  1' % constraints for homooligomer loaded, docking follows homooligomer route
            constraints.multiplicity=str2double(get(handles.valMultiplicity,'String'));
            v0=[handles.al_fit, handles.be_fit, handles.x_fit, handles.y_fit];
            dockStat0.fitst=[v0(1:2), 360/constraints.multiplicity, v0(3:4), 0];
            v0(1:2)=v0(1:2)*pi/180;
            tic
            [out_v,best_rms,exitflag,output]=fminsearch(@(v) rms_sites_homdimer(constraints,sites,v),v0,fitOptions);
            fittime=toc;
            fittimestr=sprintf('pure fitting time %g sec',fittime);
            add_msg_board(fittimestr);
            % display the result in the docking window
            [~,~]=disp_update_MMM(handles,handles.bestfit_gama,360/constraints.multiplicity,'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_x,out_v(3),'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_y,out_v(4),'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_z,0,'%0.2f');
                        
%             set(handles.bestfit_gama,'String',num2str(360/constraints.multiplicity*(180/pi)));
%             set(handles.bestfit_x,'String',num2str(out_v(3)));
%             set(handles.bestfit_y,'String',num2str(out_v(4)));
%             set(handles.bestfit_z,'String',num2str(0));
            bestfitonly=[out_v(1:2)*180/pi, 360/constraints.multiplicity, out_v(3:4), 0, best_rms];
            [~,dockStat0.distBest]=rmsDist_sites_homdimer(constraints,sites,out_v);
            
        case '1  0' % constraints for homooligomer loaded, docking follows heterodimer route
            v0=[handles.al_fit, handles.be_fit, handles.ga_fit, handles.x_fit, handles.y_fit, handles.z_fit];
            dockStat0.fitst=v0;
            v0(1:3)=v0(1:3)*pi/180;
            sites.NOsite2=sites.NOsite1;
            tic
            [out_v,best_rms,exitflag,output]=fminsearch(@(v) rms_sites_hetdimer(constraints,sites,v),v0,fitOptions);
            fittime=toc;
            fittimestr=sprintf('pure fitting time %g sec',fittime);
            add_msg_board(fittimestr);
            [~,~]=disp_update_MMM(handles,handles.bestfit_gama,out_v(3)*180/pi,'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_x,out_v(4),'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_y,out_v(5),'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_z,out_v(6),'%0.2f');
            
%             set(handles.bestfit_gama,'String',num2str(out_v(3)*180/pi));
%             set(handles.bestfit_x,'String',num2str(out_v(4)));
%             set(handles.bestfit_y,'String',num2str(out_v(5)));
%             set(handles.bestfit_z,'String',num2str(out_v(6)));
            bestfitonly=[out_v(1:3)*180/pi, out_v(4:6), best_rms];
            [~,dockStat0.distBest]=rmsDist_sites_hetdimer(constraints,sites,out_v);
%             keyboard
        case '2  0' % constraints for heterodimer loaded, docking follows heterodimer route
            v0=[handles.al_fit, handles.be_fit, handles.ga_fit, handles.x_fit, handles.y_fit, handles.z_fit];
            dockStat0.fitst=v0;
            v0(1:3)=v0(1:3)*pi/180;
            tic
            [out_v,best_rms,exitflag,output]=fminsearch(@(v) rms_sites_hetdimer(constraints,sites,v),v0,fitOptions);
            fittime=toc;
            fittimestr=sprintf('pure fitting time %g sec',fittime);
            add_msg_board(fittimestr);
            [~,~]=disp_update_MMM(handles,handles.bestfit_gama,out_v(3)*180/pi,'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_x,out_v(4),'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_y,out_v(5),'%0.2f');
            [~,~]=disp_update_MMM(handles,handles.bestfit_z,out_v(6),'%0.2f');
%             
%             set(handles.bestfit_gama,'String',num2str(out_v(3)*180/pi));
%             set(handles.bestfit_x,'String',num2str(out_v(4)));
%             set(handles.bestfit_y,'String',num2str(out_v(5)));
%             set(handles.bestfit_z,'String',num2str(out_v(6)));
            bestfitonly=[out_v(1:3)*180/pi, out_v(4:6), best_rms];
            [~,dockStat0.distBest]=rmsDist_sites_hetdimer(constraints,sites,out_v);
%             keyboard
    end
    [~,~]=disp_update_MMM(handles,handles.bestfit_alpha,out_v(1)*180/pi,'%0.2f');
    [~,~]=disp_update_MMM(handles,handles.bestfit_beta,out_v(2)*180/pi,'%0.2f');
    
%     set(handles.bestfit_alpha,'String',num2str(out_v(1)*180/pi));
%     set(handles.bestfit_beta,'String',num2str(out_v(2)*180/pi));
    dockStat0.bestfitonly=bestfitonly;
    bestfitmsg=sprintf('bestfit model: %4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f',bestfitonly(1),bestfitonly(2),bestfitonly(3),bestfitonly(4),bestfitonly(5),bestfitonly(6));
    add_msg_board(bestfitmsg);
    add_msg_board(sprintf('at %4.4f rms',bestfitonly(7)));
    
else
    grid.alpha=linspace(handles.al_grid_st,handles.al_grid_end,handles.al_grid_elem); % still in grad!
    grid.beta=linspace(handles.be_grid_st,handles.be_grid_end,handles.be_grid_elem);
    grid.gama=linspace(handles.ga_grid_st,handles.ga_grid_end,handles.ga_grid_elem);
    grid.x=linspace(handles.x_grid_st,handles.x_grid_end,handles.x_grid_elem);
    grid.y=linspace(handles.y_grid_st,handles.y_grid_end,handles.y_grid_elem);
    grid.z=linspace(handles.z_grid_st,handles.z_grid_end,handles.z_grid_elem); 
    grid.grSize=handles.al_grid_elem*handles.be_grid_elem*handles.ga_grid_elem*handles.x_grid_elem*handles.y_grid_elem*handles.z_grid_elem;
    dockStat0.gridin=grid;
    
    switch dockcase
        case '1  1' % constraints for homooligomer loaded, docking follows homooligomer route
%             keyboard
            constraints.multiplicity=str2double(get(handles.valMultiplicity,'String'));
            grid0.alpha=grid.alpha*pi/180;
            grid0.beta=grid.beta*pi/180;
            grid0.gamaFixed=grid.gama*pi/180;
            grid0.x=grid.x;
            grid0.y=grid.y;
            tic
            if get(handles.LargeGridcheckbox,'Value')
                stackSize=handles.stackSize;
                [gridoutInd,min_grid0]=homdimer_grid_large(constraints,sites,grid0,stackSize);
                gridout=grid_index2values(gridoutInd,grid,1);
                dockStat0.largeGrid=1;
                dockStat0.stackSize=stackSize;
            else
                [gridoutInd,min_grid0]=homdimer_grid(constraints,sites,grid0);
                gridout=grid_index2values(gridoutInd,grid,1);
            end
            docktime=toc;
            docktimestr=sprintf('pure docking time %g sec',docktime);
            add_msg_board(docktimestr);
            gridoutStat(:,1:2)=gridout(:,1:2);
            gridoutStat(:,3)=ones(length(gridout(:,1)),1)*360/constraints.multiplicity;
            gridoutStat(:,4:5)=gridout(:,3:4);
            gridoutStat(:,6)=zeros(length(gridout(:,1)),1);
            gridoutStat(:,7)=gridout(:,5);
            dockStat0.gridout=gridoutStat;
            out_v=[grid.alpha(min_grid0(1)), grid.beta(min_grid0(2)), grid.x(min_grid0(3)), grid.y(min_grid0(4))];
            out_v(1:2)=out_v(1:2)*pi/180;
            if get(handles.autoFitbox,'Value')
                v0=out_v;
                [out_v,best_rms,exitflag,output]=fminsearch(@(v) rms_sites_homdimer(constraints,sites,v),v0,fitOptions);
                fitpostdock=[out_v(1:2)*180/pi, 360/constraints.multiplicity, out_v(3:4), 0, best_rms];
                dockStat0.fitpostdock=fitpostdock;
            end
            [~,dockStat0.distBest]=rmsDist_sites_homdimer(constraints,sites,out_v);                        

        case '1  0' % constraints for homooligomer loaded, docking follows heterodimer route
            grid0.alpha=grid.alpha*pi/180;
            grid0.beta=grid.beta*pi/180;
            grid0.gama=grid.gama*pi/180;
            grid0.x=grid.x;
            grid0.y=grid.y;
            grid0.z=grid.z;
            sites.NOsite2=sites.NOsite1;
            tic
            if get(handles.LargeGridcheckbox,'Value')
                stackSize=handles.stackSize;
                [gridoutInd,min_grid0]=hetdimer_grid_large(constraints,sites,grid0,stackSize);
                gridout=grid_index2values(gridoutInd,grid);
                dockStat0.largeGrid=1;
                dockStat0.stackSize=stackSize;
            else
                [gridoutInd,min_grid0]=hetdimer_grid(constraints,sites,grid0);
                gridout=grid_index2values(gridoutInd,grid);
            end
            docktime=toc;
            docktimestr=sprintf('pure docking time %g sec',docktime);
            add_msg_board(docktimestr);
            dockStat0.gridout=gridout;
            out_v=[grid.alpha(min_grid0(1)), grid.beta(min_grid0(2)), grid.gama(min_grid0(3)), grid.x(min_grid0(4)), grid.y(min_grid0(5)), grid.z(min_grid0(6))];
            out_v(1:3)=out_v(1:3)*pi/180;
            if get(handles.autoFitbox,'Value')
                v0=out_v;
                [out_v,best_rms,exitflag,output]=fminsearch(@(v) rms_sites_hetdimer(constraints,sites,v),v0,fitOptions);
                fitpostdock=[out_v(1:3)*180/pi, out_v(4:6), best_rms];
                dockStat0.fitpostdock=fitpostdock;
            end
            [~,dockStat0.distBest]=rmsDist_sites_hetdimer(constraints,sites,out_v);  
            
        case '2  0' % constraints for heterodimer loaded, docking follows heterodimer route
%             [grid_select,min_grid0]=hetdimer_grid_large(constraints,sites,grid,N);
            grid0=grid;
            grid0.alpha=grid.alpha*pi/180;
            grid0.beta=grid.beta*pi/180;
            grid0.gama=grid.gama*pi/180;
            tic
            if get(handles.LargeGridcheckbox,'Value')
                stackSize=handles.stackSize;
                [gridoutInd,min_grid0]=hetdimer_grid_large(constraints,sites,grid0,stackSize);
                gridout=grid_index2values(gridoutInd,grid);
                dockStat0.stackSize=stackSize;
            else
                [gridoutInd,min_grid0]=hetdimer_grid(constraints,sites,grid0);
                gridout=grid_index2values(gridoutInd,grid);
            end
            docktime=toc;
            docktimestr=sprintf('pure docking time %g sec',docktime);
            add_msg_board(docktimestr);
            dockStat0.gridout=gridout;
            out_v=[grid.alpha(min_grid0(1)), grid.beta(min_grid0(2)), grid.gama(min_grid0(3)), grid.x(min_grid0(4)), grid.y(min_grid0(5)), grid.z(min_grid0(6))];
            out_v(1:3)=out_v(1:3)*pi/180;
            if get(handles.autoFitbox,'Value')
                v0=out_v;
                [out_v,best_rms,exitflag,output]=fminsearch(@(v) rms_sites_hetdimer(constraints,sites,v),v0,fitOptions);
                fitpostdock=[out_v(1:3)*180/pi, out_v(4:6), best_rms];
                dockStat0.fitpostdock=fitpostdock;
            end
            [~,dockStat0.distBest]=rmsDist_sites_hetdimer(constraints,sites,out_v);
    end
%     keyboard
    bestdockmsg=sprintf('best grid model: %4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f',dockStat0.gridout(1,1),dockStat0.gridout(1,2),dockStat0.gridout(1,3),dockStat0.gridout(1,4),dockStat0.gridout(1,5),dockStat0.gridout(1,6));
%     [~,dockStat0.distBest]=rmsDist_sites_hetdimer(constraints,sites,out_v);
            
    add_msg_board(bestdockmsg);
    add_msg_board(sprintf('at %4.4f rms',dockStat0.gridout(1,7)));
    if get(handles.autoFitbox,'Value')
        bestfitmsg=sprintf('fitted model: %4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f',dockStat0.fitpostdock(1:6));
    
        add_msg_board(bestfitmsg);
        add_msg_board(sprintf('at %4.4f rms',dockStat0.fitpostdock(7)));
    end;
end

% for ii=1:numel(dockBasis.strind)
%     set_object(dockBasis.strind(ii),'hide');
% end
for ii=1:numel(model.structure_ids) % hide everything displayed before
    set_object(ii,'hide');
end
axes(hMain.axes_model);
cla reset;
axis equal
axis off
set(gca,'Clipping','off');
set(gcf,'Renderer','opengl');
hold on
hMain.camlight=camlight;
guidata(hMain.axes_model,hMain);
hMain.virgin=0;
if gca~=hMain.axes_model,
    axes(hMain.axes_model);
end;
% handles.strcount=0;


acttime=datestr(clock,'dd-mm-yyyy-HH-MM-SS');
switch dockcase
    case '1  1'
%         keyboard
        multiplicity=constraints.multiplicity;
        str0=model.structures{dockBasis.strind(1)};
        allxyz1new=str0(1).xyz{1};
        [allxyz1new0,~]=get_XYZcentre(allxyz1new);
        eulang1=[out_v(1:2),0];
%         eulang1=[min_grid.alm, min_grid.bem,0]*pi/180;
        Rp1=get_rotmat(eulang1);
        t1=[0 0 0]; % no translation for the first molecule in the dimer
        t2=[out_v(3:4),0];
%         t2=[min_grid.xm, min_grid.ym, 0]; % translation of second moelcule (already in the "dimer frame")
        
        eulang2=[out_v(1:2),360/multiplicity*pi/180];
%         eulang2=[min_grid.alm, min_grid.bem, grid0.gamaFixed]*pi/180;
        Rp2=get_rotmat(eulang2);
%         allxyz1new=rot_trans(allxyz1new0,Rp1,t1);
%         allxyz2new=rot_trans(allxyz1new0,Rp2,t2);
        transformed_str{1}=rot_trans(allxyz1new0,Rp1,t1);
        transformed_str{2}=rot_trans(allxyz1new0,Rp2,t2);
        
%         rotmat=affine('Euler',-eulang1);
%         transmat = affine('translation',t1);
%         transformed_str{1} = affine_trafo_coor(allxyz1new0,{rotmat,transmat});
%         rotmat=affine('Euler',-eulang2);
%         transmat = affine('translation',t2);
%         transformed_str{2} = affine_trafo_coor(allxyz1new0,{rotmat,transmat});
        
%             eulang2=[grid0.gamaFixed,0,0]; % rotation of the second molecule (basically, a reflection about the c2 axis)
%             Rp2=get_rotmat(eulang2);
%             allxyz1new=rot_trans(allxyz1new0,Rp1,t1);
%             allxyz2new=rot_trans(allxyz1new,Rp2,t2);
       for ii=1:2
           trsnum0=numel(model.structure_ids)+1;
           trs_name0=sprintf('trs%i',trsnum0);
           [snum0,message]=copy_structure(dockBasis.strind(1),trs_name0);
        %     model.structures{snum0}.xyz{1}=transformed_str{ii};
           struct0=model.structures{snum0};
           struct0(1).xyz{1}=transformed_str{ii};
           model.structures{snum0}=struct0;
           set_object(snum0,'show',{'ribbon'});
           if get(handles.save_pdb,'Value')
%                save_pdb0=sprintf('transformed_%i',ii);
               save_pdb0=sprintf('%s_docked_mol%i_%s',dockBasis.PDB(1).pdbid,ii,acttime);
               wr_pdb_docking(save_pdb0,trs_name0);
           end
           dockStat0.str2save{ii}=trs_name0;
        end

    case '1  0'
        str0=model.structures{dockBasis.strind(1)};
        allxyz1new=str0(1).xyz{1};
%             [allxyz1new,~]=get_XYZcentre(allxyz1new);
%             allxyz2new=allxyz1new;
        [transformed_str{1},~]=get_XYZcentre(allxyz1new);
        transformed_str{2}=transformed_str{1};
        t=out_v(4:6);
        eulang=out_v(1:3);
%         t=[min_grid.xm,min_grid.ym,min_grid.zm]; % translation of second molecule, Angstroems
%         eulang=[min_grid.alm,min_grid.bem,min_grid.gam]*pi/180;   % rotation of second molecule, grad
        Rp2=get_rotmat(eulang); % rotation matrix for second molecule
        transformed_str{2}=rot_trans(transformed_str{2},Rp2,t);   % transform molecule2
        for ii=1:2
            trsnum0=numel(model.structure_ids)+1;
            trs_name0=sprintf('trs%i',trsnum0);
            [snum0,message]=copy_structure(dockBasis.strind(1),trs_name0);
        %     model.structures{snum0}.xyz{1}=transformed_str{ii};
            struct0=model.structures{snum0};
            struct0(1).xyz{1}=transformed_str{ii};
            model.structures{snum0}=struct0;
            set_object(snum0,'show',{'ribbon'});
            if get(handles.save_pdb,'Value')
%                 save_pdb0=sprintf('transformed_%i',ii);
                save_pdb0=sprintf('%s_docked_mol%i_%s',dockBasis.PDB(1).pdbid,ii,acttime);
                wr_pdb_docking(save_pdb0,trs_name0);
            end
            dockStat0.str2save{ii}=trs_name0;
         end
            
    case '2  0'
        allxyz1new=model.structures{dockBasis.strind(1)}.xyz{1};
        [transformed_str{1},~]=get_XYZcentre(allxyz1new);
        allxyz2new=model.structures{dockBasis.strind(2)}.xyz{1};
        [transformed_str{2},~]=get_XYZcentre(allxyz2new);
        t=out_v(4:6);
        eulang=out_v(1:3);
    %         t=[min_grid.xm,min_grid.ym,min_grid.zm]; % translation of second molecule, Angstroems
    %         eulang=[min_grid.alm,min_grid.bem,min_grid.gam]*pi/180;   % rotation of second molecule, grad
        Rp2=get_rotmat(eulang); % rotation matrix for second molecule
%         allxyz2new=rot_trans(allxyz2new,Rp2,t);   % transform molecule2
        transformed_str{2}=rot_trans(transformed_str{2},Rp2,t);   % transform molecule2
        for ii=1:2
            trsnum0=numel(model.structure_ids)+1;
            trs_name0=sprintf('trs%i',trsnum0);
            [snum0,message]=copy_structure(dockBasis.strind(ii),trs_name0);
            %     model.structures{snum0}.xyz{1}=transformed_str{ii};
            struct0=model.structures{snum0};
            struct0(1).xyz{1}=transformed_str{ii};
            model.structures{snum0}=struct0;
            set_object(snum0,'show',{'ribbon'});
            if get(handles.save_pdb,'Value')
%                 save_pdb0=sprintf('transformed_%i',ii);
                
                save_pdb0=sprintf('%s_docked_mol%i_%s',dockBasis.PDB(ii).pdbid,ii,acttime);
                wr_pdb_docking(save_pdb0,trs_name0);
            end
            dockStat0.str2save{ii}=trs_name0;
        end
end

camlight(hMain.camlight);
camlookat(hMain.axes_model);

dockStat0.constrFile=dockBasis.constrFile;
% acttime=datestr(clock,'dd-mm-yyyy-HH-MM-SS');
defname=['Docking_session_',acttime,'.txt'];
dockStat0.acttime=acttime;
templstr='';
for ii=1:length(dockBasis.PDB)
    templstr0=dockBasis.PDB(ii).pdbid;
    templstr=[templstr,' ',templstr0];
end
dockStat0.pdbtemplates=templstr;
dockStat0.dockcase=dockcase;
dockStat0.fitORdock=fitORdock;
if get(handles.save_stats,'Value')
    [statID,statError]=fopen(defname,'w');
    wr_docking_session(statID,dockStat0);
    fclose(statID);
end
handles.dockStat0=dockStat0;

cd(my_dir);
set(handles.finishexit,'Enable','on');
set(hfig,'Pointer','arrow');

guidata(hObject, handles);

% --- Executes on button press in Cancel_pushbutton.
function Cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global hMain % check if needed!!!!!
% % think also if destroying any accidentally hMain.docking and hMain.DEER -
% % to avoid any possible future inconsistencies!
% if
% keyboard
delete(handles.docking_window);

% --- Executes on button press in finishexit.
function finishexit_Callback(hObject, eventdata, handles)
% hObject    handle to finishexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global general

if ~get(handles.save_stats,'Value') | ~get(handles.save_pdb,'Value')
    saveexit=questdlg('Save the last docking run before exit?','Save before exit');
    switch saveexit
        case 'Yes'
            my_dir=pwd;
            cd(general.restraint_files);
            dockStat0=handles.dockStat0;
            acttime=datestr(clock,'dd-mm-yyyy-HH-MM-SS');
            defname=['Docking_session_',acttime,'.txt'];
            dockStat0.acttime=acttime;
            [statID,statError]=fopen(defname,'w');
            hfig=gcf;
            set(hfig,'Pointer','watch');
            drawnow;
            wr_docking_session(statID,dockStat0);
            fclose(statID);
            for ii=1:numel(dockStat0.str2save)
                save_pdb0=sprintf('transformed_%i',ii);
                trs_name0=dockStat0.str2save{ii};
                wr_pdb_docking(save_pdb0,trs_name0);
            end
            set(hfig,'Pointer','arrow');
            cd(my_dir);
            delete(handles.docking_window);
        case 'No'
            delete(handles.docking_window);
            
        case 'Cancel'
            return
    end
else
    delete(handles.docking_window);
end

% --- Executes on button press in save_pdb.
function save_pdb_Callback(hObject, eventdata, handles)
% hObject    handle to save_pdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_pdb

% --- Executes on button press in append2model.
function append2model_Callback(hObject, eventdata, handles)
% hObject    handle to append2model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of append2model


% --- Executes on button press in save_mat_file.
function save_mat_file_Callback(hObject, eventdata, handles)
% hObject    handle to save_mat_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_mat_file



function [newxyz,shift]=get_XYZcentre(xyz)
%
% Transforms coordinates, so that the origin is geometrical centre (average
% coordinate)
% (basically, the same as the original center_protein.m function)
[m,n]=size(xyz);
cent=sum(xyz)/m;
shift=kron(ones(m,1),cent);
newxyz=xyz-shift;


% --- Executes on button press in pushbutton_help_docking.
function pushbutton_help_docking_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_docking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'docking.html');
webcall(entry,'-helpbrowser');

function [DEER,cancelled]=process_DEER_restraints(restraints)
% [DEER,cancelled]=process_DEER_restraints(restraints)
%
% Processes the DEER constraints in a list obtained with rd_restraints.m
% checks, whether rotamers for the requested label type and temperature
% were already computed, computes them, if not
% rotamers are not actually attached
% residue indices into the template structure, template label coordinates
% (mean N-O midpoint coordinates), and label position r.m.s.d. are computed
% and stored
%
% restraints    constraint structure (obtained with rd_restraints.m)
%
% DEER          array of result structure, length is number of constraints,
%               fields are:
%               .adr1   address of the first site
%               .ind1   residue indices of the first site
%               .adr2   address of the second site
%               .ind2   residue indices of the second site
%               .r      constraint on mean distance
%               .sigr   uncertainty of constraint
%               .xyz1   Cartesian coordinates of mean N-O midpoint of first
%                       site
%               .xyz2   Cartesian coordinates of mean N-O midpoint of
%                       second site
%               .rmsd1  r.m.s.d of mean N-O midpoint coordinate of first
%                       site
%               .rmsd2  r.m.s.d. of mean N-O midpoint coordinate of
%                       second site
% cancelled     flag that reporst on interactive cancellation of
%               processing, true cancelled by user, false not cancelled
%
% G. Jeschke 2010-2012
% (Y. Polyhach adjusted 08/2013)

global model
global hMain

cancelled=false;

if ~isfield(restraints,'DEER'),
    DEER=[];
    return;
end;

% % original version guje (before 190713)-----------------------
% snum=model.current_structure;

% if isfield(restraints,'PDB'),
%     if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
%         button = questdlg(sprintf('Constraint file specifies template %s, while current template is %s. Do you want to continue?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
%         if strcmp(button,'No'),
%             cancelled=true;
%             DEER=[];
%             return
%         end;
%     end;
% end; 

% yepo adjustment 220713---------
all_structures_str=textscan(model.structure_tags,'%s','delimiter',':','MultipleDelimsAsOne',1);
all_structures=all_structures_str{1};

str1name=[];
str2name=[];
switch length(restraints.PDB)
    case 1
         if strcmpi(restraints.PDB(1).pdbid,all_structures(1))
             str1name=all_structures{1};
%              allxyz1=model.structures{1}.xyz{1};
             str0=model.structures{1};
             allxyz1=str0(1).xyz{1};
             strind=1;
             [~,shift1]=get_XYZcentre(allxyz1);
             oligomer=1; % flag needed for identification of the actual docking case later on
         else
%         if isempty(str1name)
%             button=questdlg(sprintf('The %s structure in the constraint file, does not match with any structure in the MMM hierarchy. Do you want to continue?',restraints.PDB(1),'Yes','No','No');
%              if strcmp(button,'No'),
%                  cancelled=true;
%                  DEER=[];
%                  return
%              end
%         end 
             
             warningstr=warndlg(sprintf('The %s structure in the constraint file does not match with the first structure in the MMM hierarchy. Check constraints file.',restraints.PDB(1).pdbid));
             cancelled=true;
             DEER=[];
             return
         end
    case 2
        strind=zeros(1,2);
        oligomer=2;
        if length(all_structures)==1 & length(restraints.PDB)==2
            warningstr=warndlg(sprintf('There are two structures defined in the constraints file \n while only one structure is available in the MMM hierarchy.\n Docking is interrapted. Check constraints file.'));
            cancelled=true;
            DEER=[];
            return
        end
        for kk=1:2
            if strcmpi(restraints.PDB(1).pdbid,all_structures(kk))
                str1name=all_structures{kk};
                allxyz1=model.structures{kk}.xyz{1};
                strind(1)=kk;
                [~,shift1]=get_XYZcentre(allxyz1);
            end
        end
        if isempty(str1name)
            warningstr=warndlg(sprintf('The %s structure in the constraint file does not match with either the first or the second structure in the MMM hierarchy. Check constraints file.',restraints.PDB(1).pdbid));
            cancelled=true;
            DEER=[];
            return
        end
        for kk=1:2
            if strcmpi(restraints.PDB(2).pdbid,all_structures(kk))
                str2name=all_structures{kk};
                allxyz2=model.structures{kk}.xyz{1};
                strind(2)=kk;
                [~,shift2]=get_XYZcentre(allxyz2);
            end
        end
        if isempty(str2name)
            if isempty(str1name)
                warningstr=warndlg(sprintf('The %s structure in the constraint file does not match with either the first or the second structure in the MMM hierarchy. Check constraints file.',restraints.PDB(2).pdbid));
                cancelled=true;
                DEER=[];
                return
            end
        end
    case 3 
        disp('decide which of the 3 are to be docked')
        return
end
% -----------------------------------------------------------------

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

% check whether sites are already labelled and whether all restraint sites
% do exist
% identity of the label is checked
% labeling temperature is NOT checked

% keyboard

T_list=zeros(1,200);
if ~isempty(labels),
    lindices=zeros(length(labels),4);
    for k=1:length(labels),
        cindices=labels(k).indices;
        if ~isempty(cindices),
            lindices(k,:)=cindices;
        end;
    end;
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
%         adr1=restraints.DEER(k).adr1;
        adr1_0=restraints.DEER(k).adr1;
        adr1=['[',str1name,']',adr1_0];
        ind1=resolve_address(adr1);
        if isempty(ind1),
            add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
            
            % yepo modidified 220713 --------------------------------------
%             add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum))); 
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(1)));
            % -------------------------------------------------------------
            
            add_msg_board('Processing of DEER constraints cancelled');
            cancelled=true;
            DEER=[];
            return;
        end;
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,adr1));
            end;
        end;
        if oligomer==2
    %         adr2=restraints.DEER(k).adr2;
            adr2_0=restraints.DEER(k).adr2;
            adr2=['[',str2name,']',adr2_0];
            ind2=resolve_address(adr2);
            if isempty(ind2),
                add_msg_board(sprintf('ERROR: Constraint %i has second label at site %s',k,adr2));

                % yepo modidified 220713 --------------------------------------
    %             add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
                add_msg_board(sprintf('This site does not exist in current template %s',mk_address(2)));
                % -------------------------------------------------------------

                add_msg_board('Processing of DEER constraints cancelled');
                cancelled=true;
                DEER=[];
                return;
            end;
            found=false;
            for l=1:length(labels),
                diff=ind2-lindices(l,:);
                if sum(abs(diff))==0  && strcmpi(labels(l).name,restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                for l=1:length(to_do_list),
                    if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                        found=true;
                    end;
                end;
                if ~found,
                    poi=poi+1;
                    to_do_list{poi}=adr2;
                    label_list{poi}=restraints.DEER(k).label;
                    T_list(poi)=restraints.DEER(k).T;
                    add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,adr2));
                end;
            end;
        end;
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1_0=restraints.DEER(k).adr1;
        adr1=['[',str1name,']',adr1_0];
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            label_list{poi}=restraints.DEER(k).label;
            T_list(poi)=restraints.DEER(k).T;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
        if oligomer==2
            adr2_0=restraints.DEER(k).adr2;
            adr2=['[',str2name,']',adr2_0];
            found=false;
            for l=1:length(to_do_list),
                if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr2;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr2));
            end;
        end
    end;
end;

% keyboard

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;


labels=label_information(model.sites);

% yepo adjustment 220713------------
np=length(restraints.DEER);
label_pairs=zeros(np,5);
label_pairs(:,5)=(1:np)'; % the first two columns are not used!!
NOsite1=zeros(np,5); % the 5th column is not used for sure, the 4th is not used yet (and may be will never be used!)
NOsite2=zeros(np,5);
%--------------------------------

DEER(length(restraints.DEER)).ind1=[]; % initilize DEER structure
DEER(length(restraints.DEER)).adr1=[];
DEER(length(restraints.DEER)).r=[];
DEER(length(restraints.DEER)).sigr=[];
DEER(length(restraints.DEER)).xyz1=[];
DEER(length(restraints.DEER)).rmsd1=[];
for k=1:length(restraints.DEER),
    adr1_0=restraints.DEER(k).adr1;
    adr1=['[',str1name,']',adr1_0];
    ind1=resolve_address(adr1);
%     adr2_0=restraints.DEER(k).adr2;
%     adr2=['[',str2name,']',adr2_0];
%     ind2=resolve_address(adr2);
    DEER(k).r=restraints.DEER(k).r;
    DEER(k).sigr=restraints.DEER(k).sigr;
    DEER(k).ind1=ind1;
%     DEER(k).ind2=ind2;
    DEER(k).adr1=adr1;
%     DEER(k).adr2=adr2;
    f1=false;
    for l=1:length(labels),
        diff1=ind1-labels(l).indices;
        if sum(abs(diff1))==0,
            f1=true;
            DEER(k).xyz1=labels(l).xyz;
            DEER(k).rmsd1=labels(l).rmsd;
        end
    end
    label_pairs(k,3)=10*restraints.DEER(k).r;       % to make sure dist's are already in Angstroems!!!
    label_pairs(k,4)=10*restraints.DEER(k).sigr;    % to make sure dist's are already in Angstroems!!!
    NOsite1(k,1:3)=DEER(k).xyz1;
    NOsite1(k,4)=DEER(k).rmsd1;
%     NOsite2(k,1:3)=DEER(k).xyz2;
%     NOsite2(k,4)=DEER(k).rmsd2;
    if oligomer==2
    %     adr1_0=restraints.DEER(k).adr1;
    %     adr1=['[',str1name,']',adr1_0];
    %     ind1=resolve_address(adr1);
        adr2_0=restraints.DEER(k).adr2;
        adr2=['[',str2name,']',adr2_0];
        ind2=resolve_address(adr2);
    %     DEER(k).r=restraints.DEER(k).r;
    %     DEER(k).sigr=restraints.DEER(k).sigr;
    %     DEER(k).ind1=ind1;
        DEER(k).ind2=ind2;
    %     DEER(k).adr1=adr1;
        DEER(k).adr2=adr2;
        f2=false;
        for l=1:length(labels),
            diff2=ind2-labels(l).indices;
            if sum(abs(diff2))==0,
                f2=true;
                DEER(k).xyz2=labels(l).xyz;
                DEER(k).rmsd2=labels(l).rmsd;
            end
        end
        NOsite2(k,1:3)=DEER(k).xyz2;
        NOsite2(k,4)=DEER(k).rmsd2;
        if ~f1 || ~f2,
            add_msg_board('ERROR: Automatic rotamer computation error.');
            add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
            cancelled=true;
            DEER=[];
            return;
        end
    end
    % --------------------------------------------------------------------
    if ~f1
        add_msg_board('ERROR: Automatic rotamer computation error.');
        add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
        cancelled=true;
        DEER=[];
        return;
    end
end

% yepo adjustment 220713----------
dockBasis.pairs_exp=label_pairs;
dockBasis.sites.NOsite1=NOsite1(:,1:3)-shift1(1:np,:);
if oligomer==2
    dockBasis.sites.NOsite2=NOsite2(:,1:3)-shift2(1:np,:);
end
dockBasis.strind=strind;
dockBasis.oligomer=oligomer;
DEER(1).dockBasis=dockBasis;
%---------------------------------
function labels=label_information(sites)

global model
global label_defs

poi=0;
for k0=1:length(sites),
    for k1=1:length(sites{k0}),
        for k=1:length(sites{k0}(k1).residue),
            poi=poi+1;
            labels(poi).indices=sites{k0}(k1).residue(k).indices;
            id=tag2id(sites{k0}(k1).residue(k).label,label_defs.restags);
            labels(poi).name=label_defs.residues(id).short_name;
            labels(poi).T=sites{k0}(k1).residue(k).T;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            x=sum(NOpos(:,1).*NOpos(:,4));
            y=sum(NOpos(:,2).*NOpos(:,4));
            z=sum(NOpos(:,3).*NOpos(:,4));
            labels(poi).xyz=[x y z];
            labels(poi).rmsd=NOpos_rmsd(NOpos);
        end;
    end;
end;

function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for ? -> nm
