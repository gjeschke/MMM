function varargout = FRET_window(varargin)
% FRETWINDOW M-file for FRET_window.fig
%      FRET_WINDOW, by itself, creates a new FRET_WINDOW or raises the existing
%      singleton*.
%
%      H = FRET_WINDOW returns the handle to a new FRET_WINDOW or the handle to
%      the existing singleton*.
%
%      FRET_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRET_WINDOW.M with the given input arguments.
%
%      FRET_WINDOW('Property','Value',...) creates a new FRET_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deer_window_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deer_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deer_window

% Last Modified by GUIDE v2.5 14-Jul-2020 15:10:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deer_window_OpeningFcn, ...
                   'gui_OutputFcn',  @deer_window_OutputFcn, ...
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


% --- Executes just before deer_window is made visible.
function deer_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to deer_window (see VARARGIN)

% Choose default command line output for deer_window
handles.output = hObject;
handles.eventdata = eventdata;


handles.updated = 0;
handles.exp_FRET_efficiency = 0.5;
set(handles.edit_exp_FRET_efficiency,'String',sprintf('%5.3f',handles.exp_FRET_efficiency));
handles.sel_distr = 1;
handles.range = [2,8];
handles.chromophore_pair = [];
handles.labels = [];
handles.tweak_distr = [];
handles.tweak_rax = [];
handles.expanded = false;
handles.site_list = [];
handles.rsim = [];
handles.flex_color = [0,0.5,0.5];
handles.fit_color = [0,0.5,0];

set(handles.text_mean,'String','no distribution');
set(handles.text_stddev,'String',' ');
set(handles.text_coor,'String',' ');

handles.rmsd = 0;
handles.dr = 0;

% copy control flags
handles.copy_FRET = 0;
handles.copy_distr = 0;

% find all chromophore sites that have been computed, fill listbox
[handles,success] = mk_chromophore_list(handles);

s = load('helpicon.mat');
set(handles.pushbutton_help,'CData',s.cdata);

set(handles.edit_exp_FRET_efficiency,'Enable','on');

if ~success
    msgbox('Model must feature at least two chromophores. Use FRET rotamer computation first.','FRET simulation impossible','error');
    FRET_CloseRequestFcn(handles.FRET, eventdata, handles);
else
    guidata(hObject, handles);
end
% UIWAIT makes deer_window wait for user response (see UIRESUME)
% uiwait(handles.FRET);


% --- Outputs from this function are returned to the command line.
function varargout = deer_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = hObject;
varargout{1} = eventdata;
varargout{3} = handles;


% --- Executes on selection change in listbox_label.
function listbox_label_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_label contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_label

global model

handles.eventdata = eventdata;

sel = get(hObject,'Value');
kscan = handles.site_list(sel).item(1);
kres = handles.site_list(sel).item(2);
NOpos = model.sites{kscan}.residue(kres).NOpos;
pop = NOpos(:,4);
pop = pop/sum(pop);
x = sum(NOpos(:,1).*pop);
y = sum(NOpos(:,2).*pop);
z = sum(NOpos(:,3).*pop);
name = handles.site_list(sel).adr;
msg{1} = name;
codedstring = '\u212b';
AA = sprintf(strrep(codedstring, '\u', '\x'));
msg{2} = sprintf('at [%6.2f,%6.2f,%6.2f] %s',x,y,z,AA);
set(handles.text_coor,'String',msg);
addresses = get(hObject,'String');
adr = addresses{sel};
poi = strfind(adr,';');
if ~isempty(poi)
    adr = adr(1:poi-1);
end
[msg,absresnum] = get_object(adr,'absresnum');
if msg.error
    absresnum = 0;
end
set(handles.text_abs_number,'String',sprintf('%i',absresnum));

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


% --- Executes on button press in pushbutton_select_label.
function pushbutton_select_label_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
which_one = get(handles.listbox_label,'Value');
handles.labels = which_one;
handles.ff_multi = [];
handles.tweak_distr = [];
handles.tweak_rax = [];
set(handles.text_selected_distribution,'String','*** none ***');
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_add_label.
function pushbutton_add_label_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
which_one=get(handles.listbox_label,'Value');
if isempty(handles.labels)
    handles.labels= which_one;
    handles.ff_multi = [];
    handles.tweak_distr = [];
    handles.tweak_rax = [];
else
    curr_labels = handles.labels;
    double_sel = false;
    if min(abs(curr_labels-which_one)) == 0
        double_sel = true;
    end
    if ~double_sel
        handles.labels = [curr_labels which_one];
        handles.ff_multi = [];
        handles.tweak_distr = [];
        handles.tweak_rax = [];
    else
        msgbox('The same label cannot be selected twice.','Double selection ignored','warn');
    end
    poi=length(handles.labels);
    if poi == 2
        if handles.site_list(handles.labels(1)).chain == handles.site_list(handles.labels(2)).chain % same chain
            set(handles.checkbox_coupled_ensemble,'Value',1);
        else
            set(handles.checkbox_coupled_ensemble,'Value',0);
        end
    end
end
set(handles.text_selected_distribution,'String','*** none ***');
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_select_all.
function pushbutton_select_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
addresses = get(handles.listbox_label,'String');
handles.labels = 1:length(addresses);
handles.ff_multi = [];
handles.tweak_distr = [];
handles.tweak_rax = [];
poi=length(handles.labels);
if poi == 2
    ind1 = resolve_address(addresses{handles.labels(1)});
    ind2 = resolve_address(addresses{handles.labels(2)});
    if abs(sum(ind1(1:2)-ind2(1:2))) == 0 % same chain
        set(handles.checkbox_coupled_ensemble,'Value',1);
    else
        set(handles.checkbox_coupled_ensemble,'Value',0);
    end
end
set(handles.text_selected_distribution,'String','*** none ***');
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_remove_label.
function pushbutton_remove_label_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.evendata = eventdata;
new_labels= zeros(1,10);
which_one = get(handles.listbox_label,'Value');
addresses = get(handles.listbox_label,'String');
adr = addresses{which_one};
npoi=0;
removal=0;
if ~isempty(handles.labels)
    poi=length(handles.labels);
    for k=1:poi
        if strcmp(handles.site_list(handles.labels(k)).adr,adr)
            removal=1;
        else
            npoi=npoi+1;
            new_labels(npoi)=handles.labels(k);
        end % check for double selection of an existing label
    end
end
if removal
    handles.labels = new_labels(1:npoi);
    handles.ff_multi = [];
    handles.tweak_distr = [];
    handles.tweak_rax = [];
    if npoi == 2
        ind1 = resolve_address(handles.site_list(handles.labels(1)).adr);
        ind2 = resolve_address(handles.site_list(handles.labels(2)).adr);
        if abs(sum(ind1(1:2)-ind2(1:2))) == 0 % same chain
            set(handles.checkbox_coupled_ensemble,'Value',1);
        else
            set(handles.checkbox_coupled_ensemble,'Value',0);
        end
    end
else
    msgbox('Only a selected label can be deselected.','Deselection ignored','warn');
end
set(handles.text_selected_distribution,'String','*** none ***');
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general

handles.eventdata = eventdata;

my_path=pwd;
cd(general.DEER_files);

label_list = get(handles.listbox_label,'String');
if isfield(handles,'rsim')
    labstring='';
    for k=1:length(handles.labels)
        labstring=sprintf('%s%s_',labstring,label_list{handles.labels(k)});
    end
    suggestion=[labstring 'res.txt']; 
    [fname,pname]=uiputfile('*.txt','Save results',suggestion);
    if isequal(fname,0) || isequal(pname,0)
        add_msg_board('Saving of distance distribution and FRET simulation canceled by user.');
        return;
    end
    reset_user_paths(pname);
    general.DEER_files=pname;
    % Remove (last) extension, if any
    s=strfind(fname,'.');
    if ~isempty(s)
        fname=fname(1:s(length(s))-1);
    end
    % Remove suffix '_res', if present
    s=strfind(fname,'_res');
    if ~isempty(s)
        fname=fname(1:s(length(s))-1);
    end
    % ### the following function call needs to be replaced with a save
    % function for a FRET computation ###
    handles = save_FRET_result_MMM(handles,fname,pname);

else
    msgbox('Select two chromophores before saving.','Nothing to save','warn')
end

cd(my_path);
guidata(hObject,handles);

% --- Executes on button press in checkbox_diffusion_model.
function checkbox_diffusion_model_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_diffusion_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_diffusion_model

if hObject.Value
    handles.text_static_dynamic.String = 'Dynamic computation';
    handles.text_static_dynamic.Tooltip = 'Relative diffusion of chromophores is modelled';
else
    handles.text_static_dynamic.String = 'Static computation';
    handles.text_static_dynamic.Tooltip = 'Relative diffusion of chromophores is not considered';    
end
handles = update_plots(handles);
guidata(hObject,handles);




% --- Executes on button press in pushbutton_detach_distr.
function pushbutton_detach_distr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_distr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
handles.copy_distr = 1;
handles = update_plots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FRET_CloseRequestFcn(handles.FRET, eventdata, handles);


% --- Executes when user attempts to close FRET.
function FRET_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to FRET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

% ### if you want  user confirmation before closing or automatic saving of
% results, this needs to be done here before delete(hObject)
delete(hObject);


function edit_exp_FRET_efficiency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_exp_FRET_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_exp_FRET_efficiency as text
%        str2double(get(hObject,'String')) returns contents of edit_exp_FRET_efficiency as a double

handles.eventdata = eventdata;
[exp_E,handles] = edit_update_MMM(handles,hObject,0,1,0.5,'%5.3f',0);
handles.exp_FRET_efficiency = exp_E;
handles = update_plots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_exp_FRET_efficiency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_exp_FRET_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = update_plots(handles)
% Display update

global model

col1=[0,1,0];
col2=[0,0,1];

if handles.copy_distr
    figure(1); clf;
    set(gca,'FontSize',16);
else
    axes(handles.axes_distribution);
    cla;
    if isfield(handles,'left_crsr')
        handles=rmfield(handles,'left_crsr');
    end
    if isfield(handles,'right_crsr')
        handles=rmfield(handles,'right_crsr');
    end
end
hold on;

xlabel('distance (nm');
ylabel('Probability density (1/nm)');
if ~isempty(handles.labels) && ~isempty(handles.site_list)
    lab_string = handles.site_list(handles.labels(1)).adr;
    if length(handles.labels)>1
        for k=2:length(handles.labels)
            lab_string=sprintf('%s; %s',lab_string,handles.site_list(handles.labels(k)).adr);
        end
    end
    set(handles.text_chromophore_pair,'String',lab_string);
else
    set(handles.text_chromophore_pair,'String','');
end    

numlab=length(handles.labels);

do_ensemble = get(handles.checkbox_ensemble,'Value');
if numlab > 2 && do_ensemble
    add_msg_board('Warning: Ensemble computation not implemented for multi-spin effects.');
    do_ensemble = false;
    set(handles.checkbox_ensemble,'Value',0);
    set(handles.checkbox_coupled_ensemble,'Enable','off');
end

rmin0=1e6;
rmax0=0;
% Display distance distributions
if numlab>1
    poi=0;
    pairs=numlab*(numlab-1)/2; % number of spin pairs
    for k=1:numlab-1
        currlabel = handles.labels(k);
        % find the Calpha coordinate
        switch handles.site_list(currlabel).type
            case 'label'
                adr1 = cut_address(handles.site_list(currlabel).adr);     
                dadr1 = adr1;
                if do_ensemble
                    ind1 = resolve_address(adr1);
                    ne1 = length(model.structures{ind1(1)}(ind1(2)).xyz);
                    ind1e = ind1;
                    CA1 = [0,0,0];
                    found = false;
                    for ke = 1:ne1
                        ind1e(3) = ke;
                        adr1e = mk_address(ind1e);
                        [~,CA01]=get_object(sprintf('%s.CA',adr1e),'coor');
                        if isempty(CA01) % for nucleotides
                            [~,CA01]=get_object(sprintf('%s.C1''',adr1e),'coor');
                        end
                        if ~isempty(CA01)
                            CA1 = CA1 + CA01;
                            found = true;
                        end
                    end
                    if found
                        CA1 = CA1/ne1;
                    else
                        CA1 = [];
                    end
                else
                    [~,CA1]=get_object(sprintf('%s.CA',adr1),'coor');
                    if isempty(CA1) % for nucleotides
                        [~,CA1]=get_object(sprintf('%s.C1''',adr1),'coor');
                    end
                end
        end
        for kk=k+1:numlab
            poi=poi+1;
            if pairs>1
                cf=(poi-1)/(pairs-1);
            else
                cf=0;
            end
            col=(1-cf)*col1+cf*col2;
            col=col/max(col);
            currlabel2 = handles.labels(kk);
            % find the Calpha coordinate
            switch handles.site_list(currlabel2).type
                case 'label'
                    adr2= cut_address(handles.site_list(currlabel2).adr);
                    dadr2 = adr2;
                    if do_ensemble
                        ind2 = resolve_address(adr1);
                        ne2 = length(model.structures{ind2(1)}(ind2(2)).xyz);
                        ind2e = ind2;
                        CA2 = [0,0,0];
                        found = false;
                        for ke = 1:ne2
                            ind2e(3) = ke;
                            adr2e = mk_address(ind2e);
                            [~,CA02]=get_object(sprintf('%s.CA',adr2e),'coor');
                            if isempty(CA02) % for nucleotides
                                [~,CA02]=get_object(sprintf('%s.C1''',adr2e),'coor');
                            end
                            if ~isempty(CA02)
                                CA2 = CA2 + CA02;
                                found = true;
                            end
                        end
                        if found
                            CA2 = CA2/ne2;
                        else
                            CA2 = [];
                        end
                    else
                        [~,CA2]=get_object(sprintf('%s.CA',adr2),'coor');
                        if isempty(CA2) % for nucleotides
                            [~,CA2]=get_object(sprintf('%s.C1''',adr2),'coor');
                        end
                    end
            end
            if ~isempty(CA1) && ~isempty(CA2)
                r0=norm(CA1-CA2)/10; % C-alpha/C-alpha distance in nm
                if r0-0.1<rmin0, rmin0=r0-0.1; end
                if r0+0.1>rmax0, rmax0=r0+0.1; end
            else
                r0=[];
            end
            adr1 = dadr1;
            adr2 = dadr2;
            handles.pairs{poi}=[adr1 '-' adr2];
            set(handles.FRET,'Pointer','watch');
            if do_ensemble
               [rax,act_distr] = pair_distribution(handles,adr1,adr2,currlabel,currlabel2);
               act_distr = zeros(size(act_distr));
               coupled = get(handles.checkbox_coupled_ensemble,'Value');
               if coupled % use same coordinate set numbers for both labels
                   for ke = 1:ne1
                       [~,act_distr0] = pair_distribution(handles,adr1,adr2,currlabel,currlabel2,ke,ke);
                       if ~isempty(act_distr0)
                        act_distr = act_distr + act_distr0;
                       end
                   end
                   act_distr = act_distr/ne1;
               else % combine labels from all pairs of coordinate set numbers
                   for ke = 1:ne1
                       for ke2 = 1:ne2
                           [~,act_distr0] = pair_distribution(handles,adr1,adr2,currlabel,currlabel2,ke,ke2);
                           if ~isempty(act_distr0)
                            act_distr = act_distr + act_distr0;
                           end
                       end
                   end
                   act_distr = act_distr/(ne1*ne2);
               end
            else
                [rax,act_distr] = pair_distribution(handles,adr1,adr2,currlabel,currlabel2);
            end
            set(handles.FRET,'Pointer','arrow');
            h=plot(rax,act_distr,'Color',col);
            set(h,'ButtonDownFcn',{@pair_ButtonDownFcn,poi});
            if isempty(r0)
                ndis=act_distr/sum(act_distr);
                r0=sum(rax.*ndis);
            end
            handles.pair_plots{poi}=h;
            handles.pair_distributions{poi}=act_distr;
            handles.pair_r0{poi}=r0;
            h=plot([r0,r0],[0,max(act_distr)],':','Color',col);
            handles.pair_Ca_plots{poi}=h;
            if poi==1
                full_distr = act_distr;
            else
                full_distr = full_distr+act_distr;
            end
        end
    end
    handles.rsim=rax;
    handles.dsim=full_distr;
    ndistr=full_distr/sum(full_distr);
    rmean=sum(ndistr.*rax);
    set(handles.text_mean,'String',sprintf('<r> = %4.2f nm',rmean));
    dr=rax-rmean;
    stddev=sqrt(sum(dr.^2.*ndistr));
    set(handles.text_stddev,'String',sprintf('sr = %4.2f nm',stddev));
    
    % ### Here you need to compute minimimum and optimum R0 from the
    % distance distribution ###
    tmax_min = 2*((rmean+stddev)/5)^3; % ### DEER example
    handles.R0_min = tmax_min; % ### DEER example
    tmax_opt = 2*((rmean+stddev)/3.2)^3; % ### DEER example
    handles.R0_opt = tmax_opt; % ### DEER example
    set(handles.text_minimum_R0,'String',sprintf('%3.1f us',tmax_min));
    set(handles.text_optimum_R0,'String',sprintf('%3.1f us',tmax_opt));
    
    rfilled=rax(full_distr>0.01*max(full_distr));
    rmin=min([min(rfilled) rmin0]);
    rmax=max([max(rfilled) rmax0]);
    bsl=0.1*(rmax-rmin);
    rmin=rmin-bsl;
    rmax=rmax+bsl;
    dmax=1.1*max(full_distr);
    plot(rax,full_distr,'r','LineWidth',1);
    if ~isempty(handles.tweak_distr) && ~isempty(handles.tweak_rax)
        plot(handles.tweak_rax,handles.tweak_distr,'r:','LineWidth',2);
        dmax=1.1*max([max(full_distr) max(handles.tweak_distr)]);
    end
    dmin=-0.09*dmax;
    if ~handles.expanded
        if rmin> handles.range(1)
            handles.range(1)=rmin;
        end
        if rmax< handles.range(2)
            handles.range(2)=rmax;
        end
    end
    handles = range_update(handles);
    if handles.expanded
        axis([handles.range(1),handles.range(2),dmin,dmax]);
    else
        axis([rmin,rmax,dmin,dmax]);
    end
    if isfield(handles,'AV_rax') && isfield(handles,'AVdistr')...
            && handles.checkbox_display_AV.Value
        sc = max(full_distr)/max(handles.AV_distr);
        plot(handles.AV_rax,sc*handles.AV_distr,':','Color',[0 0.75 0.25],'LineWidth',1.5);
    end
end

if handles.copy_FRET
    figure(2); clf;
    set(gca,'FontSize',16);
    hold on;
else
    axes(handles.axes_FRET);
    cla;
    hold on;
end

% ### only a dummy image is plotted here, you need to plot the R0/D plot
% ###

if handles.checkbox_plot_rmsd.Value
    C = imread('ngc6543a.jpg');
    image(C);
    xlabel('R0 (nm)');
    ylabel('Diff. const. (nm^2/ns)');
else
    % ### please check whether you want to plot something if the checkbox is
    % inactive, this would go here ###
end


handles.copy_distr = false;
handles.copy_FRET = false;

% Update handles structure
guidata(handles.listbox_label, handles);

function [handles,success] = mk_chromophore_list(handles)
% Fills listbox_label with the list of available labels
% returns a success flag success=1 if there are at least two labels,
% success=0 otherwise

global model

% ### for testing only, we use nitroxide labels ###
% this needs to be replace by 'chromophore'
label_class = 'nitroxide';

success=0;
sitelistpoi = 0;
label_list = cell(1,1000); % more than generous initialization
site_list(1000).type = '';
site_list(1000).item = [];
site_list(1000).chain = [];
site_list(1000).adr = '';

% return immediately, if no site scan had been performed
if ~isfield(model,'sites') || isempty(model.sites)
    return
end

for kscan = 1:length(model.sites) % loop over all site scans
    if strcmpi(model.sites{kscan}.class,label_class) % only for selected label class
        for kres = 1:length(model.sites{kscan}.residue) % loop over sites in scan
            sitelistpoi = sitelistpoi + 1;
            adr = mk_address(model.sites{kscan}.residue(kres).indices,true);
            label_list{sitelistpoi} = adr;
            site_list(sitelistpoi).type = 'label';
            site_list(sitelistpoi).item = [kscan,kres];
            site_list(sitelistpoi).chain = model.sites{kscan}.residue(kres).indices(2);
            site_list(sitelistpoi).adr = adr;
            site_list(sitelistpoi).label = model.sites{kscan}.residue(kres).label;
        end
    end
end
label_list = label_list(1:sitelistpoi);
site_list = site_list(1:sitelistpoi);

if sitelistpoi >= 2
    success=1;
else
    set(handles.listbox_label,'String',' ');
    set(handles.listbox_label,'Value',1);
    return;
end

handles.site_list = site_list;

set(handles.listbox_label,'String',label_list);
set(handles.listbox_label,'Value',1);

if sitelistpoi > 0
    kscan = site_list(1).item(1);
    kres = site_list(1).item(2);
    NOpos = model.sites{kscan}.residue(kres).NOpos;
    pop = NOpos(:,4);
    pop = pop/sum(pop);
    x = sum(NOpos(:,1).*pop);
    y = sum(NOpos(:,2).*pop);
    z = sum(NOpos(:,3).*pop);
    name = site_list(1).adr;
    msg{1}=name;
    codedstring = '\u212b';
    AA = sprintf(strrep(codedstring, '\u', '\x'));
    msg{2}=sprintf('at [%6.2f,%6.2f,%6.2f] %s',x,y,z,AA);
end
set(handles.text_coor,'String',msg);

function [rax,distr] = pair_distribution(handles,adr1,adr2,lab1,lab2,km1,km2)
% Distance distribution function for a pair of spin labels
% km1, km2 replace the model index, if present

global model

label1 = handles.site_list(lab1).label;
label2 = handles.site_list(lab2).label;

ind1 = resolve_address(adr1);

if exist('km1','var')
    ind1(3) = km1;
end

ind2 = resolve_address(adr2);

if exist('km2','var')
    ind2(3) = km2;
end

rax = [];
distr = [];

sig=0.1; % Gaussian broadening of the distance distribution

options.trajectory=[0,0];

NOpos1 = [];
NOpos2 = [];

for kscan = 1:length(model.sites)
    for kres = 1:length(model.sites{kscan}.residue)
        if strcmpi(label1,model.sites{kscan}.residue(kres).label)
            di = model.sites{kscan}.residue(kres).indices - ind1;
            if sum(abs(di)) == 0
                NOpos1 =  model.sites{kscan}.residue(kres).NOpos;
            end
        end
        if strcmpi(label2,model.sites{kscan}.residue(kres).label)
            di = model.sites{kscan}.residue(kres).indices - ind2;
            if sum(abs(di)) == 0
                NOpos2 =  model.sites{kscan}.residue(kres).NOpos;
            end
        end
    end
end

if exist('NOpos1','var') && exist('NOpos2','var')
    [rax,distr]=get_distribution(NOpos1,NOpos2,sig,[],[],[],options);
end

% --- Executes on button press in pushbutton_detach_FRET.
function pushbutton_detach_FRET_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_detach_FRET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
handles.copy_FRET = 1;
handles=update_plots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

handles.eventdata = eventdata;
entry=strcat(help_files,'deer_window.html');
webcall(entry,'-helpbrowser');
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_range_analysis.
function pushbutton_range_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_range_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

numlab=length(handles.labels);

if numlab~=2
    add_msg_board('Error: Range analysis requires exactly two labels.');
    return
else
    label_list = get(handles.listbox_label,'String');
    adr1 = label_list{handles.labels(1)};
    adr2 = label_list{handles.labels(2)};
    set(handles.FRET,'Pointer','watch');
    [~,fname] = analyze_range(handles,adr1,adr2);
    set(handles.FRET,'Pointer','arrow');
    if ~isempty(fname)
        hMain.report_file=fname;
        report_editor;
    end
end
% Update handles structure
guidata(hObject, handles);

function edit_lower_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower as text
%        str2double(get(hObject,'String')) returns contents of edit_lower as a double

handles.eventdata = eventdata;
if isempty(handles.rsim)
    add_msg_board('Error: Range can be updated only if distribution is available.');
else
    [v,handles] = edit_update_MMM(handles,hObject,min(handles.rsim),handles.range(2)-0.05,handles.range(1),'%4.2f',0);
    handles.range = [v,handles.range(2)];
    handles = range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end

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

handles.eventdata = eventdata;
if isempty(handles.rsim)
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(1)-0.05 >= min(handles.rsim)
    handles.range(1) = handles.range(1)-0.05;
    handles = range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_lower_right.
function pushbutton_lower_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lower_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
if isempty(handles.rsim)
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(1)+0.05 <= handles.range(2)-0.05
    handles.range(1) = handles.range(1)+0.05;
    handles = range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end


function edit_upper_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_upper as text
%        str2double(get(hObject,'String')) returns contents of edit_upper as a double

handles.eventdata = eventdata;

if isempty(handles.rsim)
    add_msg_board('Error: Range can be updated only if distribution is available.');
else
    [v,handles] = edit_update_MMM(handles,hObject,handles.range(1)+0.05,max(handles.rsim),handles.range(2),'%4.2f',0);
    handles.range = [handles.range(1),v];
    handles = range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end

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

handles.eventdata = eventdata;

if isempty(handles.rsim)
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(2)-0.05 >= handles.range(1)+0.05
    handles.range(2) = handles.range(2)-0.05;
    handles = range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on button press in pushbutton_upper_right.
function pushbutton_upper_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_upper_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;

if isempty(handles.rsim)
    add_msg_board('Error: Range can be updated only if distribution is available.');
elseif handles.range(2)+0.0<max(handles.rsim)
    handles.range(2) = handles.range(2)+0.05;
    handles=range_update(handles);
    % Update handles structure
    guidata(hObject, handles);
end

function handles = range_update(handles)

if isfield(handles,'left_crsr') && ishandle(handles.left_crsr)
    set(handles.left_crsr,'XData',[handles.range(1) handles.range(1)]);
else
    axes(handles.axes_distribution);
    handles.left_crsr = plot([handles.range(1) handles.range(1)],[min(handles.dsim),max(handles.dsim)],'b','LineWidth',0.5);
end
set(handles.edit_lower,'String',sprintf('%4.2f',handles.range(1)));

if isfield(handles,'right_crsr') && ishandle(handles.right_crsr)
    set(handles.right_crsr,'XData',[handles.range(2) handles.range(2)]);
else
    axes(handles.axes_distribution);
    handles.right_crsr = plot([handles.range(2) handles.range(2)],[min(handles.dsim),max(handles.dsim)],'m','LineWidth',0.5);
end
set(handles.edit_upper,'String',sprintf('%4.2f',handles.range(2)));


function [pairs,fname]=analyze_range(handles,adr1,adr2)
% Analyze which rotamer pairs contribute to distance distribution in the
% range speciefied by handles.range

global model
global general

pairs=[];
fname='';

rotamers1=0;
rotamers2=0;

library1='';
library2='';
if isfield(model,'labels')
    for k=1:length(model.labels)
        if strcmp(model.labels(k).adr,adr1)
            NOpos1=model.labels(k).NOpos;
            [significant1,n1]=size(NOpos1);
            rotamers1=true;
            [library1,NOpos1] = find_library_NOpos(adr1);
        end
        if strcmp(model.labels(k).adr,adr2)
            NOpos2=model.labels(k).NOpos;
            [significant2,n2]=size(NOpos2);
            rotamers2=true;
            [library2,NOpos2] = find_library_NOpos(adr2);
        end
    end
end
if ~isempty(library1)
    char_table1 = rotamer_char_table(library1);
end
if ~isempty(library2)
    char_table2 = rotamer_char_table(library2);
end

if ~rotamers1 && ~rotamers2
    add_msg_board('Error: For range analysis, at least one site must be a spin label with rotamers');
    return
end

[pairs,popsum] = pairs_in_range(handles.range,NOpos1,NOpos2); % range must be supplied in nm

if isempty(pairs)
    add_msg_board('No rotamer pairs contribute to distance distribution in selected range.');
else
    frac=sum(pairs(:,3))/popsum;
    add_msg_board(sprintf('Selected range corresponds to %6.1f%% of total pair population',100*frac));
    [sortpop,poppoi]=sort(pairs(:,3),'descend');
    sortpop=sortpop/sum(sortpop);
    [~,significant]=min(abs(cumsum(sortpop)-0.99));
    add_msg_board(sprintf('%i pairs contribute in this range, %i of those contribute 99%% population of the range',length(poppoi),significant));
    spairs=pairs(poppoi(1:significant),:); % leading 99% of pairs sorted by descending population
    poi1=0;
    [m1,~]=size(NOpos1);
    rot1=zeros(m1,2);
    for k=1:m1 % analyze site 1
        cpairs=find(spairs(:,1)==NOpos1(k,5));
        if ~isempty(cpairs)
            pops=sum(spairs(cpairs,3));
            poi1=poi1+1;
            rot1(poi1,1)=NOpos1(k,5);
            rot1(poi1,2)=pops;
        end
    end
    rot1=rot1(1:poi1,:);
    add_msg_board(sprintf('At site 1: %i out of %i rotamers contribute in selected range.',poi1,m1));
    poi2=0;
    [m2,~]=size(NOpos2);
    rot2=zeros(m2,2);
    for k=1:m2 % analyze site 2
        cpairs=find(spairs(:,2)==NOpos2(k,5));
        if ~isempty(cpairs)
            pops=sum(spairs(cpairs,3));
            poi2=poi2+1;
            rot2(poi2,1)=NOpos2(k,5);
            rot2(poi2,2)=pops;
        end
    end
    rot2=rot2(1:poi2,:);
    add_msg_board(sprintf('At site 2: %i out of %i rotamers contribute in selected range.',poi2,m2));
    fname=[general.tmp_files 'chi1_chi2_analysis.dat'];
    fid=fopen(fname,'w');
    fprintf(fid,'--- MMM chi1/chi2 analysis of rotamers ---\nSite 1: %s, Site 2: %s\ncontributions to the distance range [%4.2f, %4.2f] nm\n\n',adr1,adr2,handles.range(1),handles.range(2));
    fprintf(fid,'Selected range corresponds to %6.1f%% of total pair population\n',100*frac);
    fprintf(fid,'%i pairs contribute in this range,\n%i of those contribute 99%% population of the range\n\n',length(poppoi),significant);
    fprintf(fid,'At %s: %i out of %i rotamers contribute in selected range.\n',adr1,poi1,m1);
    if ~isempty(library1)
        fprintf(fid,'The chi1/chi2 analysis of rotamers at %s follows:\n\n',adr1);
        fprintf(fid,'rotamer   total population (%%)   pop. in range (%%)\n');
        [pops,chi1_chi2]=chi1_chi2_analysis(NOpos1,rot1,char_table1);
        [m,~]=size(pops);
        for k=1:m
            fprintf(fid,'(%s) %15.1f%20.1f\n',chi1_chi2{k},100*pops(k,1)/sum(pops(:,1)),100*pops(k,2)/sum(pops(:,2)));
        end
    end
    fprintf(fid,'\nAt %s: %i out of %i rotamers contribute in selected range.\n',adr2,poi2,m2);
    if ~isempty(library2)
        [pops,chi1_chi2]=chi1_chi2_analysis(NOpos2,rot2,char_table2);
        fprintf(fid,'The chi1/chi2 analysis of rotamers at %s follows:\n\n',adr2);
        fprintf(fid,'rotamer   total population (%%)   pop. in range (%%)\n');
        [m,~]=size(pops);
        for k=1:m
            fprintf(fid,'(%s) %15.1f%20.1f\n',chi1_chi2{k},100*pops(k,1)/sum(pops(:,1)),100*pops(k,2)/sum(pops(:,2)));
        end
    end
    fclose(fid);
    % Now color rotamers in range
    rotlist=NOpos1(:,5);
    poplist=NOpos1(:,4);
    [~,poppoi]=sort(poplist,'descend');
    attached1=rotlist(poppoi(1:significant1)); % list of numbers of significant rotamers
    for k=1:length(attached1)
        rotamer=attached1(k);
        in_range=sum(find(spairs(:,1)==rotamer));
        if in_range
            loctag=get_loctag(k);
            locadr=sprintf('%s.:%s',adr1,loctag);
            set_object(locadr,'color',{[0,0,0]});
        end
    end
    rotlist=NOpos2(:,5);
    poplist=NOpos2(:,4);
    [~,poppoi]=sort(poplist,'descend');
    attached2=rotlist(poppoi(1:significant2)); % list of numbers of significant rotamers
    for k=1:length(attached2)
        rotamer=attached2(k);
        in_range=sum(find(spairs(:,2)==rotamer));
        if in_range
            loctag=get_loctag(k);
            locadr=sprintf('%s.:%s',adr2,loctag);
            set_object(locadr,'color',{[0,0,0]});
        end
    end
end

function [library,NOpos] = find_library_NOpos(adr)
% returns the rotamer library for a labeled residue at a given address
% empty string is returned, if residue is not found or library was not
% stored

global model

library='';
NOpos=[];
indices=resolve_address(adr);
if isempty(indices) || length(indices)~=4,
    return;
end;


nscan=length(model.sites);
for k=1:nscan,
    if isfield(model.sites{k},'library'),
        currlib=model.sites{k}.library;
    else
        currlib='';
    end;    
    for kc=1:length(model.sites{k}),
        for kk=1:length(model.sites{k}(kc).residue),
            diff=sum(abs(indices-model.sites{k}(kc).residue(kk).indices));
            if diff==0,
                library=currlib;
                NOpos=model.sites{k}(kc).residue(kk).NOpos;
            end
        end;
    end;
end;

function [pops,chi1_chi2]=chi1_chi2_analysis(NOpos,rot,char_table)
% Populations of chi1-chi2 rotamer groups in full distribution NOpos and
% range-selected distribution rot, the character table char_table for the
% corresponding rotamer library must be provided

chi1_chi2{1}='m,m';
chi1_chi2{2}='m,t';
chi1_chi2{3}='m,p';
chi1_chi2{4}='t,m';
chi1_chi2{5}='t,t';
chi1_chi2{6}='t,p';
chi1_chi2{7}='p,m';
chi1_chi2{8}='p,t';
chi1_chi2{9}='p,p';

pops=zeros(9,2);

[m,n]=size(NOpos);

for k=1:m,
    indy=NOpos(k,5);
    character=char_table{indy};
    short_char=character(1:3);
    for kk=1:9,
        if strcmp(short_char,chi1_chi2{kk}),
            pops(kk,1)=pops(kk,1)+NOpos(k,4);
        end;
    end;
end;

[mr,nr]=size(rot);

for k=1:mr,
    indy=rot(k,1);
    character=char_table{indy};
    short_char=character(1:3);
    for kk=1:9,
        if strcmp(short_char,chi1_chi2{kk}),
            pops(kk,2)=pops(kk,2)+rot(k,2);
        end;
    end;
end;

function loctag=get_loctag(rotamer)
% return the location tag for a numbered rotamer

loctags0=':A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:';

locid=mod(rotamer,26);
if locid==0, locid=26; end
locmid=floor((rotamer-1)/26);
if locmid>0
    loctag=strcat(id2tag(locmid,loctags0),id2tag(locid,loctags0));
else
    loctag=id2tag(locid,loctags0);
end


% --- Executes on button press in togglebutton_expand.
function togglebutton_expand_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_expand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
state=get(hObject,'Value');
handles.expanded=state;
if state
    set(hObject,'TooltipString','Return to full view');
    set(hObject,'String','Full view');
else
    set(hObject,'TooltipString','Expand or shrink view to selected range');
    set(hObject,'String','Expand');
end

handles = update_plots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_uniform.
function pushbutton_uniform_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_uniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;

numlab=length(handles.labels);

if numlab~=2
    add_msg_board('Error: "Uniform!" requires exactly two labels.');
    return
else
    label_list = get(handles.listbox_label,'String');
    adr1=label_list{handles.labels(1)};
    adr2=label_list{handles.labels(2)};
    set(handles.FRET,'Pointer','watch');
    [rax,distr]=uniform_rotamers(adr1,adr2,handles);
    handles.tweak_rax = rax;
    handles.tweak_distr = distr;
    set(handles.FRET,'Pointer','arrow');
    handles = update_plots(handles);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_AV.
function pushbutton_AV_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
% ### extension should be adapted, multiple ###
[filename,pathname] = uigetfile('.add','Load AV distance distribution');
fname = fullfile(pathname, filename);
set(handles.FRET,'Pointer','watch');
drawnow
% ### the following must be replaced by a reader that provides distanance
% axis rax and distance distribution disrt form an AV file ###
[rax,distr] = rd_PRONOX(fname);
% ### end of section to be edited
set(handles.FRET,'Pointer','arrow');
handles.AV_rax = rax;
handles.AV_distr = distr;
handles = update_plots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_ensemble.
function checkbox_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ensemble

handles.eventdata = eventdata;
if get(hObject,'Value')
    set(handles.checkbox_coupled_ensemble,'Enable','on');
else
    set(handles.checkbox_coupled_ensemble,'Enable','off');
end
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes on button press in checkbox_coupled_ensemble.
function checkbox_coupled_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_coupled_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_coupled_ensemble

handles.eventdata = eventdata;
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes on button press in checkbox_display_AV.
function checkbox_display_AV_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_display_AV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_display_AV

handles.eventdata = eventdata;
handles = update_plots(handles);
guidata(hObject,handles);


function edit_diff_const_lower_bound_Callback(hObject, eventdata, handles)
% hObject    handle to edit_diff_const_lower_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_diff_const_lower_bound as text
%        str2double(get(hObject,'String')) returns contents of edit_diff_const_lower_bound as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_diff_const_lower_bound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_diff_const_lower_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_diff_const_upper_bound_Callback(hObject, eventdata, handles)
% hObject    handle to edit_diff_const_upper_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_diff_const_upper_bound as text
%        str2double(get(hObject,'String')) returns contents of edit_diff_const_upper_bound as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_diff_const_upper_bound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_diff_const_upper_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_exci_donor_lifetime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_exci_donor_lifetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_exci_donor_lifetime as text
%        str2double(get(hObject,'String')) returns contents of edit_exci_donor_lifetime as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_exci_donor_lifetime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_exci_donor_lifetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_diffusion_grid_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_diffusion_grid_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_diffusion_grid_size as text
%        str2double(get(hObject,'String')) returns contents of edit_diffusion_grid_size as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_diffusion_grid_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_diffusion_grid_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_diff_matrix_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_diff_matrix_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_diff_matrix_size as text
%        str2double(get(hObject,'String')) returns contents of edit_diff_matrix_size as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_diff_matrix_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_diff_matrix_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_update_model.
function pushbutton_update_model_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eventdata = eventdata;
% ### at this point you need to call the model update computation ###
% you best privide handles as an input and report back the results as field
% of structure handles
handles = update_plots(handles);
guidata(hObject,handles);


function edit_sim_FRET_efficiency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sim_FRET_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sim_FRET_efficiency as text
%        str2double(get(hObject,'String')) returns contents of edit_sim_FRET_efficiency as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_sim_FRET_efficiency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sim_FRET_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_F_radius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_F_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_F_radius as text
%        str2double(get(hObject,'String')) returns contents of edit_F_radius as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_F_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_F_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_error_F_radius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_error_F_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_error_F_radius as text
%        str2double(get(hObject,'String')) returns contents of edit_error_F_radius as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_error_F_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_error_F_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_R0_grid_points_Callback(hObject, eventdata, handles)
% hObject    handle to edit_R0_grid_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_R0_grid_points as text
%        str2double(get(hObject,'String')) returns contents of edit_R0_grid_points as a double

handles.eventdata = eventdata;
% ### &&& specify edit field &&& ###
handles = update_plots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_R0_grid_points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_R0_grid_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_fitted_diffusion_constant_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_fitted_diffusion_constant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_plot_rmsd.
function checkbox_plot_rmsd_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_rmsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_rmsd

handles.eventdata = eventdata;
handles = update_plots(handles);
guidata(hObject,handles);

function adr = cut_address(adr)
% removes residue name from address, if any

poi = strfind(adr,';');
if ~isempty(poi)
    adr = adr(1:poi-1);
end