function varargout = hierarchy_window_large(varargin)
% HIERARCHY_WINDOW_LARGE M-file for hierarchy_window_large.fig
%      HIERARCHY_WINDOW_LARGE, by itself, creates a new HIERARCHY_WINDOW_LARGE or raises the existing
%      singleton*.
%
%      H = HIERARCHY_WINDOW_LARGE returns the handle to a new HIERARCHY_WINDOW_LARGE or the handle to
%      the existing singleton*.
%
%      HIERARCHY_WINDOW_LARGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HIERARCHY_WINDOW_LARGE.M with the given input arguments.
%
%      HIERARCHY_WINDOW_LARGE('Property','Value',...) creates a new HIERARCHY_WINDOW_LARGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hierarchy_window_large_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hierarchy_window_large_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hierarchy_window_large

% Last Modified by GUIDE v2.5 25-Jul-2013 10:02:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hierarchy_window_large_OpeningFcn, ...
                   'gui_OutputFcn',  @hierarchy_window_large_OutputFcn, ...
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


% --- Executes just before hierarchy_window_large is made visible.
function hierarchy_window_large_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hierarchy_window_large (see VARARGIN)

% Choose default command line output for hierarchy_window_large
handles.output = hObject;

global hMain
global MMM_icon
% global proteopedia_icon
global PDBwiki_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

mhandles=guidata(hMain.figure);
locked=get(mhandles.uitoggletool_lock,'State');

hMain.hierarchy_display=1;
hMain.hierarchy_window_large=hObject;
hMain.hierarchy_last_selected=[];
hMain.hierarchy_atom_selected=[];

% load('down_button.mat');
% set(handles.down_button,'CData',icon_data);
% load('up_button.mat');
% set(handles.up_button,'CData',icon_data);

load helpicon
set(handles.pushbutton_help_hierarchy,'CData',cdata);

load pdbicon
set(handles.pushbutton_PDB_explore,'CData',cdata);

set(handles.pushbutton_pdbwiki,'CData',PDBwiki_icon);

% set(handles.pushbutton_EDS,'CData',proteopedia_icon);

load uniproticon
set(handles.pushbutton_UniProt,'CData',cdata);

handles=popup_setup(handles);
handles=sequence_display(handles);

if strcmpi(locked,'on'),
    set(handles.pushbutton_loop,'Enable','off');
    set(handles.pushbutton_helix,'Enable','off');
    set(handles.pushbutton_strand,'Enable','off');
else
    set(handles.pushbutton_loop,'Enable','on');
    set(handles.pushbutton_helix,'Enable','on');
    set(handles.pushbutton_strand,'Enable','on');
end;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hierarchy_window_large wait for user response (see UIRESUME)
% uiwait(handles.hierarchy_window_large);


% --- Outputs from this function are returned to the command line.
function varargout = hierarchy_window_large_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in up_button.
function up_button_Callback(hObject, eventdata, handles)
% hObject    handle to up_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function hierarchy_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to hierarchy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

handles=guidata(hObject);
set(handles.text_info,'String','All deselected.');

mouseclick=get(gcf,'SelectionType');

if strcmpi(mouseclick,'open')
    hMain.hierarchy_selected=[];
    hMain=cmd(hMain,'unselect *');
end;


function my_ButtonDownFcn(hObject,eventdata,mynum,icode)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

handles=guidata(hObject);

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;

cset=get(handles.popupmenu_chain_model,'Value');
address=sprintf('[%s](%s){%i}%i%s',stag,ctag,cset,mynum);
indices=resolve_address(address);
address=sprintf('[%s](%s){:}%i%s',stag,ctag,mynum);
resname=model.structures{snum}(indices(2)).residues{1}.info(indices(4)).name;
resnum=model.structures{snum}(indices(2)).residues{1}.info(indices(4)).number;
if isfield(model.structures{snum}(indices(2)).residues{1}.info(indices(4)),'insertion_code'),
    resicode=model.structures{snum}(indices(2)).residues{1}.info(indices(4)).insertion_code;
end;

my_handle=hMain.hierarchy_items(mynum);

mouseclick=get(gcf,'SelectionType');

if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')

    % Deselect old selection
    if ~isempty(hMain.hierarchy_selected), % was anything selected
        if mynum==hMain.hierarchy_last_selected,
            hMain.hierarchy_last_selected=[];
        end;
        hMain.hierarchy_selected=[];
    end;
    hMain=cmd(hMain,'unselect *');
end;

if strcmpi(mouseclick,'open')
    add_msg_board('--- All deselected ---');
    set(handles.text_info,'String','All deselected.');
    hMain.hierarchy_last_selected=[];
    return
end;

if strcmpi(mouseclick,'alt') % unselect click command
    hMain.hierarchy_itemsel(mynum)=0;
    if mynum==hMain.hierarchy_last_selected,
        hMain.hierarchy_last_selected=[];
    end;
    sel=[];
    if ~isempty(hMain.hierarchy_selected),
        for k=1:length(hMain.hierarchy_selected),
            if hMain.hierarchy_selected(k)~=my_handle,
                sel=[sel hMain.hierarchy_selected(k)];
            end;
        end;
    end;
    hMain.hierarchy_selected=sel;
    command=sprintf('unselect %s',address);
else % select click command
    hMain.hierarchy_itemsel(mynum)=1;
    hMain.hierarchy_last_selected=mynum;
    hMain.hierarchy_selected=[hMain.hierarchy_selected my_handle];
    command=sprintf('select %s',address);
    if sum(find(handles.incomplete==mynum))>0,
        add_msg_board(sprintf('Residue %s has missing atoms.',address));
    end;
end;
hMain=cmd(hMain,command);
handles=update_atom_selection(handles);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in down_button.
function down_button_Callback(hObject, eventdata, handles)
% hObject    handle to down_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function handles=sequence_display(handles)

global model
global hMain
global graph_settings
global residue_defs

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
cnum=tag2id(ctag,model.chain_tags{snum});
sequence=model.structures{snum}(cnum).sequence;

type=-1;
if isfield(model.structures{snum}(cnum),'seqtype'), % determine sequence type, information is present
    type=model.structures{snum}(cnum).seqtype;
    seqtest=double(sequence)>=double('a');
    if type==2 && sum(seqtest),
        type=3;
    end;
end;

switch type
    case -1
        seqtype='Sequence';
        panelback=[240,240,240]/255;
    case 0
        seqtype='Heterosequence';
        panelback=[255,255,255]/255;
    case 1
        seqtype='Peptide sequence';
        panelback=[240,240,240]/255;
    case 2
        seqtype='DNA sequence';
        panelback=[255,242,221]/255;
    case 3
        seqtype='RNA sequence';
        panelback=[245,235,235]/255;
end;
set(handles.hierarchy,'Title',sprintf('%s display - [%s](%s)',seqtype,stag,ctag));
set(handles.hierarchy,'BackgroundColor',panelback);

set(handles.hierarchy_window_large,'Pointer','watch');
drawnow

hydropathy_color=get(handles.radiobutton_hydropathy,'Value');
hydropathy_range=graph_settings.max_hydropathy-graph_settings.min_hydropathy;
helix_propensity_color=get(handles.radiobutton_helix_propensity,'Value');
helix_propensity_range=graph_settings.max_helix_propensity-graph_settings.min_helix_propensity;
charge_color=get(handles.radiobutton_charge,'Value');
DSSP_color=get(handles.radiobutton_DSSP,'Value');
if ~isfield(model.info{snum},'dssp') || isempty(model.info{snum}.dssp),
    DSSP_color=false;
    set(handles.radiobutton_DSSP,'Value',0);
    set(handles.radiobutton_DSSP,'Enable','off');
else
    set(handles.radiobutton_DSSP,'Enable','on');
end;

if isfield(model.info{snum},'resolution') && ~isempty(model.info{snum}.resolution),
    resstring=sprintf('%4.2f Å',model.info{snum}.resolution);
else
    resstring='not specified';
end;
    
set(hMain.MMM,'Name',sprintf('MMM - [%s](%s) Resolution %s',stag,ctag,resstring));


if ~isfield(model.structures{snum}(cnum),'seqexist') || length(model.structures{snum}(cnum).seqexist)~=length(sequence),
    seqexist=zeros(1,length(sequence));
    mres=length(model.structures{snum}(cnum).residues{1}.info);
    numbers=zeros(1,mres);
    for k=1:mres,
        numbers(k)=model.structures{snum}(cnum).residues{1}.info(k).number;
    end;
    for k=1:length(sequence),
        if ~isempty(find(numbers==k, 1)),
            seqexist(k)=1;
        end;
    end;
    model.structures{snum}(cnum).seqexist=seqexist;
else
    seqexist=model.structures{snum}(cnum).seqexist;
end;

info=model.structures{snum}(cnum).residues{1}.info;
% find out whether there are residues with missing atoms
incomplete=[];
if isfield(model.info{snum},'incomplete'),
    all_incomplete=model.info{snum}.incomplete;
    if ~isempty(all_incomplete),
        [m,n]=size(all_incomplete);
        for k=1:m,
            if all_incomplete(k,2)==cnum,
                num=model.structures{snum}(cnum).residues{all_incomplete(k,3)}.info(all_incomplete(k,4)).number;
                incomplete=[incomplete num];
            end;
        end;
    end;
end;
handles.incomplete=incomplete;
minres=16384;
maxres=0;
for k=1:length(info),
    if info(k).number<minres, minres=info(k).number; end;
    if info(k).number>maxres && info(k).type==1, maxres=info(k).number; end;
    if isfield(info,'terminal'),
        if info(k).terminal, break; end;
    end;
end;
lseq=length(sequence);

axes(handles.axes2);
cla
hold on

hMain.hierarchy_items=zeros(1,lseq);
hMain.hierarchy_itemsel=zeros(1,lseq);
hMain.hierarchy_selected=[];

set(handles.axes2,'Units','pixels');
posi0=get(handles.axes2,'Position');
xmax0=posi0(3);
ymax0=posi0(4);
set(handles.axes2,'Units','characters');
posi=get(handles.axes2,'Position');
xmin=0;
xmax=posi(3); %480;
ymin=0;
ymax=posi(4); %300;
lines=36;
cpl=50;
lneeded=ceil(lseq/cpl);
dw0=xmax0/xmax;
dh0=ymax0/ymax;
x00=4*dw0;
dw=1.6; % 8;
dh=1.3; %15;
x0=4*dw;

for k=1:cpl/10+1,
    text((k-1)*10*dw+x0,0,'|','Color',[0.25,0.25,0.75],'FontName','FixedWidth');
end;
for k=1:lneeded,
    text(0,k*dh,sprintf('%i',cpl*(k-1)+1),'Color',[0.25,0.25,0.75],'FontName','FixedWidth');
end;
for k=1:lneeded,
    text((cpl+1.5)*dw+x0,k*dh,sprintf('%i',cpl*k),'Color',[0.25,0.25,0.75],'FontName','FixedWidth');
end;
for k=1:cpl/10+1,
    text((k-1)*10*dw+x0,(lneeded+1)*dh,'|','Color',[0.25,0.25,0.75],'FontName','FixedWidth');
end;

poi=0;
for ln=1:lines,
    for cl=1:cpl,
        poi=poi+1;
        if poi>lseq, break; end;
        h=text(xmin+cl*dw+x0,ymin+ln*dh,sequence(poi));
        if ~seqexist(poi),
            set(h,'Color',[0.5,0.5,0.5]);
        else
            set(h,'ButtonDownFcn',{@my_ButtonDownFcn,poi});
            if sum(find(incomplete==poi))>0,
                set(h,'FontWeight','bold');
                set(h,'Color',[0,0.5,0]);
            end;
        end;
        set(h,'FontName','FixedWidth');
        if abs(type)==1,
            if hydropathy_color,
                resnum=findstr(sequence(poi),residue_defs.single_letter_code);
                if ~isempty(resnum) && resnum<=length(residue_defs.residues),
                    info=residue_defs.residues(resnum);
                    col=(info.hydropathy-graph_settings.min_hydropathy)/hydropathy_range;
                    set(h,'BackgroundColor',[col,col,(1-col)]);
                end;
            end;
            if helix_propensity_color,
                resnum=findstr(sequence(poi),residue_defs.single_letter_code);
                if ~isempty(resnum) && resnum<=length(residue_defs.residues),
                    info=residue_defs.residues(resnum);
                    col=1-sqrt((info.helix_propensity-graph_settings.min_helix_propensity)/helix_propensity_range);
                    set(h,'BackgroundColor',[col,col,(1-col)]);
                end;
            end;
            if charge_color,
                resnum=findstr(sequence(poi),residue_defs.single_letter_code);
                if ~isempty(resnum) && resnum<=length(residue_defs.residues),
                    info=residue_defs.residues(resnum);
                    switch info.charge,
                        case 2, col=[0.5,0.5,1];
                        case 1, col=[0.75,0.75,1];
                        case -1, col=[1,0.75,0.75];
                        case -2, col=[1,0.5,0.5];
                        otherwise col=[];
                    end;
                    if ~isempty(col),
                        set(h,'BackgroundColor',col);
                    end;
                end;
            end;
            if DSSP_color,
                snum=model.current_structure;
                ids=model.structure_ids;
                sid=[];
                for k=1:length(ids),
                    if ids(k)==snum,
                        sid=k;
                    end;
                end;
                stag=id2tag(sid,model.structure_tags);
                ctag=model.current_chain;

                cset=get(handles.popupmenu_chain_model,'Value');
                address=sprintf('[%s](%s){%i}%i%s',stag,ctag,cset,poi);
                indices=resolve_address(address);
                if ~isempty(indices) && length(indices)>=4 && isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'dssp'),
                    dssp_sec=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).dssp;
                else
                    col=[];
                    dssp_sec='';
                end;
                switch dssp_sec
                    case 'H'
                        col=graph_settings.DSSP_color_H;
                    case 'B'
                        col=graph_settings.DSSP_color_B;
                    case 'E'
                        col=graph_settings.DSSP_color_E;
                    case 'G'
                        col=graph_settings.DSSP_color_G;
                    case 'I'
                        col=graph_settings.DSSP_color_I;
                    case 'T'
                        col=graph_settings.DSSP_color_T;
                    case 'S'
                        col=graph_settings.DSSP_color_S;
                    case ' '
                        col=graph_settings.DSSP_color_C;
                    otherwise
                        col=[];
                end;
                if ~isempty(col),
                    set(h,'Color',col);
                end;
            end;
        end;
        hMain.hierarchy_items(poi)=h;
    end;
end;

% Plot secondary structure information

% Loops
if isfield(model.structures{snum}(cnum),'loop_defs'),
    loopnum=length(model.structures{snum}(cnum).loop_defs);
    for k=1:loopnum, 
        hrange=model.structures{snum}(cnum).loop_defs{k}.range;
        % disp(sprintf('Loop %i from %i to %i',k,hrange(1),hrange(2)));
        la=ceil(hrange(1)/cpl); % lines and columns of residues
        ca=mod(hrange(1),cpl);
        if ca==0, ca=cpl; end;
        le=ceil(hrange(2)/cpl);
        ce=mod(hrange(2),cpl);
        if ce==0, ce=cpl; end;
        while le>la, % loop over linebreak
            len=cpl-ca+1;
            x=linspace(0,len,8*len);
            h=plot(dw*x+x0+ca*dw,1.5*mod(3*x,2)/dh0+la*dh-0.6*dh,'k','LineWidth',1);
            set(h,'Color',graph_settings.coil_color,'ButtonDownFcn',{@loop_ButtonDownFcn,k});
            la=la+1;
            ca=1;
        end;
        len=ce-ca+1;
        x=linspace(0,len,8*len);
        h=plot(dw*x+x0+ca*dw,1.5*mod(3*x,2)/dh0+la*dh-0.6*dh,'k','LineWidth',1);
        set(h,'Color',graph_settings.coil_color,'ButtonDownFcn',{@loop_ButtonDownFcn,k});
    end;
end;

% Helices
if isfield(model.structures{snum}(cnum),'helix_defs'),
    helnum=length(model.structures{snum}(cnum).helix_defs);
    for k=1:helnum, 
        hrange=model.structures{snum}(cnum).helix_defs{k}.range;
        % disp(sprintf('Helix %i from %i to %i',k,hrange(1),hrange(2)));
        la=ceil(hrange(1)/cpl); % lines and columns of residues
        ca=mod(hrange(1),cpl);
        if ca==0, ca=cpl; end;
        le=ceil(hrange(2)/cpl);
        ce=mod(hrange(2),cpl);
        if ce==0, ce=cpl; end;
        while le>la, % helix over linebreak
            len=cpl-ca+1;
            x=linspace(0,len,8*len);
            y=sin(2*pi*x/3.6);
            h=plot(dw*x+x0+ca*dw,-1.5*y/dh0+la*dh-0.45*dh,'k','LineWidth',2);
            set(h,'Color',graph_settings.helix_color,'ButtonDownFcn',{@helix_ButtonDownFcn,k});
            la=la+1;
            ca=1;
        end;
        len=ce-ca+1;
        x=linspace(0,len,8*len);
        y=sin(2*pi*x/3.6);
        h=plot(dw*x+x0+ca*dw,-1.5*y/dh0+la*dh-0.45*dh','k','LineWidth',2);
        set(h,'Color',graph_settings.helix_color,'ButtonDownFcn',{@helix_ButtonDownFcn,k});
    end;
end;

% Strands
if isfield(model.structures{snum}(cnum),'sheet_defs'),
    helnum=length(model.structures{snum}(cnum).sheet_defs);
    for k=1:helnum, 
        hrange=model.structures{snum}(cnum).sheet_defs{k}.range;
        % disp(sprintf('Sheet %i from %i to %i',k,hrange(1),hrange(2)));
        la=ceil(hrange(1)/cpl); % lines and columns of residues
        ca=mod(hrange(1),cpl);
        if ca==0, ca=cpl; end;
        le=ceil(hrange(2)/cpl);
        ce=mod(hrange(2),cpl);
        if ce==0, ce=cpl; end;
        while le>la, % helix over linebreak
            len=cpl-ca+1;
            x=linspace(0,len,dw*len);
            h=plot(dw*x+x0+ca*dw,0*x+la*dh-0.50*dh,'k','LineWidth',2);
            set(h,'Color',graph_settings.sheet_color,'ButtonDownFcn',{@strand_ButtonDownFcn,k});
            la=la+1;
            ca=1;
        end;
        len=ce-ca+1;
        x=linspace(0,len,dw*len);
        h=plot(dw*x+x0+ca*dw,0*x+la*dh-0.50*dh,'k','LineWidth',2);
        set(h,'Color',graph_settings.sheet_color,'ButtonDownFcn',{@strand_ButtonDownFcn,k});
    end;
end;


axis([xmin,xmax,ymin,ymax]);

% Now check for existing selections

hMain.hierarchy_selected=[];

if isfield(model,'selections') && model.selections>0,
    for k=1:model.selections,
        indices=model.selected{k};
        indices=indices(indices>0);
        if length(indices)==4 && indices(1)==snum && indices(2)==cnum, % residue in current chain selected
            rnum=model.structures{snum}(cnum).residues{1}.info(indices(4)).number;
            inserted=0;
            if isfield(model.structures{snum}(cnum).residues{1}.info(indices(4)),'insertion_code'),
                if ~isempty(strtrim(model.structures{snum}(cnum).residues{1}.info(indices(4)).insertion_code)),
                    inserted=1;
                end;
            end;
            if rnum>0 && rnum<=length(hMain.hierarchy_items),
                if ~inserted,
                    my_handle=hMain.hierarchy_items(rnum);
                    hMain.hierarchy_selected=[hMain.hierarchy_selected my_handle];
                    hMain.hierarchy_itemsel(rnum)=1;
                end;
            end;
        end;
    end;
    if ~isempty(hMain.hierarchy_selected),
        for k=1:length(hMain.hierarchy_selected),
            set(hMain.hierarchy_selected(k),'Color',[0.75,0,0]);
        end;
        set(hMain.hierarchy_selected(end),'Color',[1,0,0]);
    end;
end;
 
set(handles.hierarchy_window_large,'Pointer','arrow');


% --- Executes when user attempts to close hierarchy_window_large.
function hierarchy_window_large_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to hierarchy_window_large (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global hMain

hMain.hierarchy_display=0;

delete(hObject);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

hMain.hierarchy_display=0;

delete(handles.hierarchy_window_large);


% --- Executes on selection change in popupmenu_structure.
function popupmenu_structure_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_structure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_structure

global model

selected=get(hObject,'Value');
snum=model.structure_ids(selected);
model.current_structure=snum;
ctag=id2tag(1,model.chain_tags{model.current_structure});
set(handles.popupmenu_chain,'Value',1);
set(handles.popupmenu_chain_model,'Value',1);
model.current_chain=ctag;
handles=popup_setup(handles);
handles=sequence_display(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_structure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_chain.
function popupmenu_chain_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_chain contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chain

global model

selected=get(hObject,'Value');
ctag=id2tag(selected,model.chain_tags{model.current_structure});
model.current_chain=ctag;
set(handles.popupmenu_chain_model,'Value',1);
handles=popup_setup(handles);
handles=sequence_display(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_chain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_chain_model.
function popupmenu_chain_model_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chain_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_chain_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chain_model


% --- Executes during object creation, after setting all properties.
function popupmenu_chain_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chain_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_atoms.
function listbox_atoms_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_atoms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_atoms contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_atoms


% --- Executes during object creation, after setting all properties.
function listbox_atoms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_atoms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_locations.
function listbox_locations_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_locations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_locations contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_locations


% --- Executes during object creation, after setting all properties.
function listbox_locations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_locations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=popup_setup(handles)

global model
global hMain

stags=model.structure_tags;
stags=stags(2:end); % remove leading delimiter
nonsense=textscan(stags,'%s','Delimiter',':');
stag_list=nonsense{1};
set(handles.popupmenu_structure,'String',stag_list);
snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
set(handles.popupmenu_structure,'Value',sid);

ctags=model.chain_tags{snum};
cnum=tag2id(model.current_chain,ctags);
ctags=ctags(2:end); % remove leading delimiter
nonsense=textscan(ctags,'%s','Delimiter',':');
ctag_list=nonsense{1};
set(handles.popupmenu_chain,'String',ctag_list);
set(handles.popupmenu_chain,'Value',cnum);

models=length(model.structures{snum}(cnum).residues);
for k=1:models,
    model_list(k) = cellstr(sprintf('%i',k));
end;
set(handles.popupmenu_chain_model,'String',model_list);
set(handles.popupmenu_chain_model,'Value',1);

info=model.structures{snum}(cnum).residues{1}.info;
poi=0;
cofacs=zeros(1,1000);
for k=1:length(info),
    if info(k).hetflag && ~strcmpi(info(k).name,'HOH'),
        poi=poi+1;
        cofac_list(poi)=cellstr(sprintf('%s %i',info(k).name,info(k).number));
        cofacs(poi)=info(k).number;
    end;
end;
if poi>0,
    hMain.hierarchy_cofactors=cofacs(1:poi);
    set(handles.popupmenu_cofactors,'String',cofac_list);
    set(handles.popupmenu_cofactors,'Value',1);
else
    hMain.hierarchy_cofactors=[];
    set(handles.popupmenu_cofactors,'String','*** No cofactors ***');
    set(handles.popupmenu_cofactors,'Value',1);
end;

% --- Executes on selection change in popupmenu_cofactors.
function popupmenu_cofactors_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cofactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_cofactors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupmenu_cofactors

% --- Executes during object creation, after setting all properties.
function popupmenu_cofactors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cofactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function loop_ButtonDownFcn(hObject,eventdata,mynum)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

handles=guidata(hObject);
snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
cnum=tag2id(ctag,model.chain_tags{snum});

address=sprintf('[%s](%s){:}<L.%s>',stag,ctag,model.structures{snum}(cnum).loop_defs{mynum}.name);

mouseclick=get(gcf,'SelectionType');

if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')
    hMain=cmd(hMain,'unselect *');
end;

if strcmpi(mouseclick,'alt') % unselect click command
    command=sprintf('unselect %s',address);
    set(handles.text_info,'String',sprintf('Loop %s deselected.',model.structures{snum}(cnum).loop_defs{mynum}.name));
else % select click command
    command=sprintf('select %s',address);
    set(handles.text_info,'String',sprintf('Loop %s selected.',model.structures{snum}(cnum).loop_defs{mynum}.name));
end;
hMain=cmd(hMain,command);

indices=resolve_address(address);
range=[min(indices(:,4)) max(indices(:,4))]; % range of residues

hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles,1);

function helix_ButtonDownFcn(hObject,eventdata,mynum)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

handles=guidata(hObject);
snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
cnum=tag2id(ctag,model.chain_tags{snum});

address=sprintf('[%s](%s){:}<H.%s>',stag,ctag,model.structures{snum}(cnum).helix_defs{mynum}.name);

mouseclick=get(gcf,'SelectionType');

if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')
    hMain=cmd(hMain,'unselect *');
end;

if strcmpi(mouseclick,'alt') % unselect click command
    command=sprintf('unselect %s',address);
    set(handles.text_info,'String',sprintf('Helix %s deselected.',model.structures{snum}(cnum).helix_defs{mynum}.name));
else % select click command
    command=sprintf('select %s',address);
    set(handles.text_info,'String',sprintf('Helix %s selected.',model.structures{snum}(cnum).helix_defs{mynum}.name));
end;
hMain=cmd(hMain,command);

indices=resolve_address(address);
range=[min(indices(:,4)) max(indices(:,4))]; % range of residues

hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles,1);

function strand_ButtonDownFcn(hObject,eventdata,mynum)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

handles=guidata(hObject);
snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
cnum=tag2id(ctag,model.chain_tags{snum});

address=sprintf('[%s](%s){:}<E.%s>',stag,ctag,model.structures{snum}(cnum).sheet_defs{mynum}.name);

mouseclick=get(gcf,'SelectionType');

if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')
    hMain=cmd(hMain,'unselect *');
end;

if strcmpi(mouseclick,'alt') % unselect click command
    command=sprintf('unselect %s',address);
    set(handles.text_info,'String',sprintf('Strand %s deselected.',model.structures{snum}(cnum).sheet_defs{mynum}.name));
else % select click command
    command=sprintf('select %s',address);
    set(handles.text_info,'String',sprintf('Strand %s selected.',model.structures{snum}(cnum).sheet_defs{mynum}.name));
end;
hMain=cmd(hMain,command);
hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles,1);

% --- Executes on button press in pushbutton_cofactor_plus.
function pushbutton_cofactor_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cofactor_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

selected=get(handles.popupmenu_cofactors,'Value');
selected=hMain.hierarchy_cofactors(selected);

hMain.hierarchy_last_selected=selected;
handles=update_atom_selection(handles);
guidata(hObject,handles);

cset=get(handles.popupmenu_chain_model,'Value');

address=mk_residue_address(selected);
command=sprintf('select %s',address);
hMain=cmd(hMain,command);

% --- Executes on button press in pushbutton_cofactor_minus.
function pushbutton_cofactor_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cofactor_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

selected=get(handles.popupmenu_cofactors,'Value');
selected=hMain.hierarchy_cofactors(selected);

address=mk_residue_address(selected);
command=sprintf('unselect %s',address);
hMain=cmd(hMain,command);

if hMain.hierarchy_last_selected==selected,
    hMain.hierarchy_last_selected=[];
    handles=update_atom_selection(handles);
    guidata(hObject,handles);
end;

% --- Executes on button press in pushbutton_cofactor_select.
function pushbutton_cofactor_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cofactor_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

selected=get(handles.popupmenu_cofactors,'Value');
selected=hMain.hierarchy_cofactors(selected);

hMain.hierarchy_last_selected=selected;
handles=update_atom_selection(handles);
guidata(hObject,handles);

hMain=cmd(hMain,'unselect *');
cset=get(handles.popupmenu_chain_model,'Value');
address=mk_residue_address(selected);
command=sprintf('select %s',address);
hMain=cmd(hMain,command);


% --- Executes on button press in pushbutton_atom_plus.
function pushbutton_atom_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_atom_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=atom_selection(hObject, eventdata, handles, 0);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_atom_minus.
function pushbutton_atom_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_atom_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

hMain.hierarchy_atom_selected=[];

if isempty(hMain.hierarchy_last_selected),
    set(handles.text_info,'String','Atom deselection impossible. Select residue or cofactor first.');
end;

selected=get(handles.listbox_atoms,'Value');

cset=get(handles.popupmenu_chain_model,'Value');

rnum=hMain.hierarchy_last_selected;
address=mk_residue_address(rnum);
indices=resolve_address(address);
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
set(handles.text_atom_list,'String',sprintf('Atoms - %s %i',info.name,info.number));
atag=id2tag(selected,info.atom_tags);
command=sprintf('unselect %s.%s',address,atag);

hMain=cmd(hMain,command);

handles=update_location_selection(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_atom_select.
function pushbutton_atom_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_atom_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=atom_selection(hObject, eventdata, handles, 1);
guidata(hObject,handles);

function handles=atom_selection(hObject, eventdata, handles, resetflag)
% Exclusive or additive atom selection
%
% resetflag 1: exclusive selection, 0: add to existing selection,
%           defaults to 0

if nargin<4, resetflag=0; end;

global model
global hMain

if isempty(hMain.hierarchy_last_selected),
    set(handles.text_info,'String','Atom selection impossible. Select residue or cofactor first.');
end;

selected=get(handles.listbox_atoms,'Value');

cset=get(handles.popupmenu_chain_model,'Value');

rnum=hMain.hierarchy_last_selected;
address=mk_residue_address(rnum);
indices=resolve_address(address);
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
set(handles.text_atom_list,'String',sprintf('Atoms - %s %i',info.name,info.number));
atag=id2tag(selected,info.atom_tags);
command=sprintf('select %s.%s',address,atag);

if resetflag,
    hMain=cmd(hMain,'unselect *');
end;

hMain=cmd(hMain,command);

indices=resolve_address([address '.' atag]);
hMain.hierarchy_atom_selected=indices;
handles=update_location_selection(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_location_plus.
function pushbutton_location_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_location_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=location_selection(hObject, eventdata, handles, 0);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_location_minus.
function pushbutton_location_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_location_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

if nargin<4, resetflag=0; end;

indices=hMain.hierarchy_atom_selected;

if isempty(indices),
    set(handles.text_info,'String','No location was selected.');
else
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
    rnum=hMain.hierarchy_last_selected;
    address=mk_residue_address(rnum);
    selected=get(handles.listbox_atoms,'Value');
    atag=id2tag(selected,info.atom_tags);
    selected=get(handles.listbox_locations,'Value');
    ltag=id2tag(selected,info.location_tags);
    address=sprintf('%s.%s:%s',address,atag,ltag);
    hMain=cmd(hMain,sprintf('unselect %s',address));
end;



% --- Executes on button press in pushbutton_location_select.
function pushbutton_location_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_location_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=location_selection(hObject, eventdata, handles, 1);
guidata(hObject,handles);

function handles=location_selection(hObject, eventdata, handles, resetflag)
% Exclusive or additive location selection
%
% resetflag 1: exclusive selection, 0: add to existing selection,
%           defaults to 0

global hMain
global model

if nargin<4, resetflag=0; end;

indices=hMain.hierarchy_atom_selected;

if isempty(indices),
    set(handles.text_info,'String','Location selection impossible. Select atom first.');
else
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
    rnum=hMain.hierarchy_last_selected;
    address=mk_residue_address(rnum);
    selected=get(handles.listbox_atoms,'Value');
    atag=id2tag(selected,info.atom_tags);
    selected=get(handles.listbox_locations,'Value');
    ltag=id2tag(selected,info.location_tags);
    address=sprintf('%s.%s:%s',address,atag,ltag);
    if resetflag,
        hMain=cmd(hMain,'unselect *');
    end;
    hMain=cmd(hMain,sprintf('select %s',address));
end;


function handles=update_atom_selection(handles,secflag)
% Updates the atom selection listbox if a new residue (or cofactor) was
% selected or the previously selected residue was unselected

global model
global hMain

if nargin<2, secflag=0; end;

if isempty(hMain.hierarchy_last_selected),
    set(handles.text_atom_list,'String','Atoms');
    set(handles.listbox_atoms,'Value',1);
    set(handles.listbox_atoms,'String','*** Nothing selected ***');
    if ~secflag,
        set(handles.text_info,'String',sprintf('No residue selected.'));
    end;
else
    rnum=hMain.hierarchy_last_selected;
    address=mk_residue_address(rnum);
    indices=resolve_address(address);
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
    set(handles.text_atom_list,'String',sprintf('Atoms - %s %i',info.name,info.number));
    atags=info.atom_tags(2:end);
    nonsense=textscan(atags,'%s','Delimiter',':');
    atag_list=nonsense{1};
    set(handles.listbox_atoms,'String',atag_list);
    set(handles.listbox_atoms,'Value',1);
    if ~secflag,
        set(handles.text_info,'String',sprintf('Residue %s-%i selected.',info.name,rnum));
    end;
end;

function handles=update_location_selection(handles)
% Updates the location selection listbox if a new residue (or cofactor) was
% selected or the previously selected residue was unselected

global model
global hMain

if isempty(hMain.hierarchy_atom_selected),
    set(handles.text_location_list,'String','Locations');
    set(handles.listbox_locations,'Value',1);
    set(handles.listbox_locations,'String','*** Nothing selected ***');
else
    indices=hMain.hierarchy_atom_selected;
    info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
    atag=id2tag(indices(5),info.atom_tags);
    set(handles.text_location_list,'String',sprintf('Locations - %s',atag));
    if isfield(info,'location_tags'),
        anum=info.atom_numbers{indices(5)};
        [m,n]=size(anum);
        if m<=1,
            set(handles.listbox_locations,'String','*** No alternate locations ***');
        else
            ltags=info.location_tags(2:end);
            nonsense=textscan(ltags,'%s','Delimiter',':');
            ltag_list=nonsense{1};
            set(handles.listbox_locations,'String',ltag_list);
            set(handles.listbox_locations,'Value',1);
        end;
    else
        set(handles.listbox_locations,'String','*** No alternate locations ***');
    end;
end;


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_hydropathy.
function radiobutton_hydropathy_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_hydropathy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=sequence_display(handles);
guidata(hObject,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_helix_propensity.
function radiobutton_helix_propensity_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_helix_propensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=sequence_display(handles);
guidata(hObject,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_nocolor.
function radiobutton_nocolor_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_nocolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=sequence_display(handles);
guidata(hObject,handles);


% --- Executes when selected object is changed in uipanel_color.
function uipanel_color_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_color 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles=sequence_display(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_cset_plus.
function pushbutton_cset_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cset_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

cset=get(handles.popupmenu_chain_model,'Value');

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('select [%s](%s){%i}',stag,ctag,cset);
hMain=cmd(hMain,command);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_cset_minus.
function pushbutton_cset_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cset_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

cset=get(handles.popupmenu_chain_model,'Value');

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('unselect [%s](%s){%i}',stag,ctag,cset);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_cset_select.
function pushbutton_cset_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cset_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

cset=get(handles.popupmenu_chain_model,'Value');

hMain=cmd(hMain,'unselect *');
hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles);

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('select [%s](%s){%i}',stag,ctag,cset);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_chain_plus.
function pushbutton_chain_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_chain_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('select [%s](%s)',stag,ctag);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_chain_minus.
function pushbutton_chain_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_chain_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('unselect [%s](%s)',stag,ctag);
hMain=cmd(hMain,command);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_chain_select.
function pushbutton_chain_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_chain_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

hMain=cmd(hMain,'unselect *');

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('select [%s](%s)',stag,ctag);
hMain=cmd(hMain,command);
hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_structure_plus.
function pushbutton_structure_plus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_structure_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
command=sprintf('select [%s]',stag);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_structure_minus.
function pushbutton_structure_minus_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_structure_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
command=sprintf('unselect [%s]',stag);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_structure_select.
function pushbutton_structure_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_structure_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain
global model

hMain=cmd(hMain,'unselect *');
hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles);

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
command=sprintf('select [%s]',stag);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_water.
function pushbutton_water_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global hMain

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
command=sprintf('select [%s](%s)"HOH"',stag,ctag);
hMain=cmd(hMain,command);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_plus_all_cofactors.
function pushbutton_plus_all_cofactors_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plus_all_cofactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

cset=get(handles.popupmenu_chain_model,'Value');

for k=1:length(hMain.hierarchy_cofactors),
    selected=hMain.hierarchy_cofactors(k);
    address=mk_residue_address(selected);
    command=sprintf('select %s',address);
    hMain=cmd(hMain,command);
end;

hMain.hierarchy_last_selected=selected;
handles=update_atom_selection(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_minus_all_cofactors.
function pushbutton_minus_all_cofactors_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_minus_all_cofactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hMain

cset=get(handles.popupmenu_chain_model,'Value');

for k=1:length(hMain.hierarchy_cofactors),
    selected=hMain.hierarchy_cofactors(k);
    address=mk_residue_address(selected);
    command=sprintf('unselect %s',address);
    hMain=cmd(hMain,command);
end;

hMain.hierarchy_last_selected=[];
handles=update_atom_selection(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_loop.
function pushbutton_loop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_object('*','secondary',{'loop'});
update_secondary_plot(hObject,handles);

% --- Executes on button press in pushbutton_strand.
function pushbutton_strand_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_object('*','secondary',{'strand'});
update_secondary_plot(hObject,handles);

% --- Executes on button press in pushbutton_helix.
function pushbutton_helix_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_object('*','secondary',{'helix'});
update_secondary_plot(hObject,handles);

function update_secondary_plot(hObject,handles)

global hMain

figure(hMain.figure);

[indices,message]=resolve_address('*');
done_sec=zeros(1000,2);
done_graph=zeros(1000,3);
poi_sec=0;
poi_graph=0;
[m,n]=size(indices);
for k=1:m,
    secflag=1;
    if poi_sec>0,
        for kk=1:poi_sec,
            if done_sec(kk,:)==indices(1:2), secflag=0; end;
        end;
    end;
    if secflag,
        poi_sec=poi_sec+1;
        done_sec(poi_sec,:)=indices(k,1:2);
        pdb_secondary(indices(k,1),indices(k,2));
    end;
    graphflag=1;
    if poi_graph>0,
        for kk=1:poi_graph,
            if done_graph(kk,:)==indices(k,1:3), graphflag=0; end;
        end;
    end;
    if graphflag,
        poi_graph=poi_graph+1;
        done_graph(poi_graph,:)=indices(k,1:3);
        set_chain_model(indices(k,1:3),'update');
    end;
end;
handles=sequence_display(handles);
guidata(hObject,handles);
highlight_selection;


% --- Executes on button press in pushbutton_help_hierarchy.
function pushbutton_help_hierarchy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_hierarchy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'hierarchy_window.html');
webcall(entry,'-helpbrowser');


% --- Executes on mouse press over figure background.
function hierarchy_window_large_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to hierarchy_window_large (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles=sequence_display(handles);
% guidata(hObject,handles);


% --- Executes on button press in pushbutton_PDB_explore.
function pushbutton_PDB_explore_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PDB_explore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global model
global queries

snum=model.current_structure;
if isfield(model.info{snum},'idCode') && ~isempty(model.info{snum}.idCode),
    url=[web_adr.PDB_explore queries.PDB_explore_id model.info{snum}.idCode];
    webcall(url);
else
    add_msg_board('Warning: No PDB identifier available.');
end;

% --- Executes on button press in pushbutton_UniProt.
function pushbutton_UniProt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_UniProt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global model

snum=model.current_structure;
ctag=model.current_chain;
cnum=tag2id(ctag,model.chain_tags{snum});
if isfield(model.structures{snum}(cnum),'dbref'),
    access=model.structures{snum}(cnum).dbref;
    [dbase,code]=strtok(access,':');
    if strcmpi(dbase,'UNP'), % for UniProt, sequences are read
        url=[web_adr.UniProt code(2:end)];
        webcall(url);
    else
        add_msg_board('Warning: UniProt access code is unknown.');
    end;
else
    add_msg_board('Warning: No sequence database information available.');
end;


% --- Executes on button press in pushbutton_EDS.
function pushbutton_EDS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_EDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global queries
global model

snum=model.current_structure;
if isfield(model.info{snum},'idCode') && ~isempty(model.info{snum}.idCode),
    url=[web_adr.EDS queries.EDS lower(model.info{snum}.idCode)];
    webcall(url);
else
    add_msg_board('Warning: No PDB identifier available.');
end;

% --- Executes on button press in pushbutton_pdbwiki.
function pushbutton_pdbwiki_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pdbwiki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global web_adr
global model

snum=model.current_structure;
if isfield(model.info{snum},'idCode') && ~isempty(model.info{snum}.idCode),
    url=[web_adr.PDBwiki lower(model.info{snum}.idCode)];
    webcall(url);
else
    add_msg_board('Warning: No PDB identifier available.');
end;
