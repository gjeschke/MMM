function varargout = magic_fit(varargin)
% MAGIC_FIT M-file for magic_fit.fig
%      MAGIC_FIT, by itself, creates a new MAGIC_FIT or raises the existing
%      singleton*.
%
%      H = MAGIC_FIT returns the handle to a new MAGIC_FIT or the handle to
%      the existing singleton*.
%
%      MAGIC_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAGIC_FIT.M with the given input arguments.
%
%      MAGIC_FIT('Property','Value',...) creates a new MAGIC_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before magic_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to magic_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help magic_fit

% Last Modified by GUIDE v2.5 03-Jul-2015 11:16:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @magic_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @magic_fit_OutputFcn, ...
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


% --- Executes just before magic_fit is made visible.
function magic_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to magic_fit (see VARARGIN)

global model
global MMM_icon

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

nstruct=length(model.structures);
for k=1:nstruct,
    stag_list{k}=id2tag(k,model.structure_tags,model.structure_ids);
end;
set(handles.popupmenu_template,'String',stag_list);
set(handles.popupmenu_to_fit,'String',stag_list);
set(handles.popupmenu_template,'Value',model.current_structure);
for k=1:nstruct,
    if k~=model.current_structure, break; end;
end;
set(handles.popupmenu_to_fit,'Value',k);

% Choose default command line output for magic_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes magic_fit wait for user response (see UIRESUME)
uiwait(handles.figure_magic);


% --- Outputs from this function are returned to the command line.
function varargout = magic_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global help_files

entry=strcat(help_files,'magic_fit.html');
webcall(entry,'-helpbrowser');

% --- Executes on selection change in popupmenu_template.
function popupmenu_template_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_template contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_template

% --- Executes during object creation, after setting all properties.
function popupmenu_template_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_to_fit.
function popupmenu_to_fit_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_to_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_to_fit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_to_fit


% --- Executes during object creation, after setting all properties.
function popupmenu_to_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_to_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
global hMain

templates=get(handles.popupmenu_template,'String');
template=templates{get(handles.popupmenu_template,'Value')};
tobefitted=get(handles.popupmenu_to_fit,'String');
tofit=tobefitted{get(handles.popupmenu_to_fit,'Value')};
mode=0;
if get(handles.radiobutton_backbone,'Value'),
    mode=1;
end;
if get(handles.radiobutton_CA,'Value'),
    mode=2;
end;
selected=get(handles.checkbox_selected,'Value');
aligned=get(handles.checkbox_aligned_residues,'Value');
whole=get(handles.checkbox_whole_structure,'Value');
set(handles.figure_magic,'Pointer','watch');
drawnow;
tag=['motion:' tofit '_to_' template];
id_target=tag2id(tofit,model.structure_tags,model.structure_ids);

ensemble_mode = get(handles.checkbox_ensemble_mode,'Value');
if ensemble_mode
    num_models = length(model.structures{id_target}(1).xyz);
    template0 = template;
    for k = 1:num_models
        [transmat,template,target]=find_magic_fit(handles,template0,tofit,mode,selected,aligned,whole,k);
        matrices{k} = transmat;
    end
else
    [matrices,template,target]=find_magic_fit(handles,template,tofit,mode,selected,aligned,whole);
end

if get(handles.checkbox_motion,'Value') && get(handles.radiobutton_CA,'Value') 
    if isfield(model,'motion')
        k=length(model.motion)+1;
    else
        k=1;
    end
    [m,n]=size(template);
    model.motion(k).tag=tag;
    model.motion(k).color=[0.47,0.53,0.60];
    model.motion(k).transparency=ones(m,1);
    model.motion(k).radius=0.4;
    model.motion(k).template=template;
    model.motion(k).target=target;
    model.motion(k).active=1;
    model.motion(k).snum=id_target;
    model.motion(k).transform=0;
end;

if ~isempty(matrices),
    id_trafo=tag2id(tofit,model.structure_tags,model.structure_ids);
    if ensemble_mode,
        transform_structure_ensemble(id_trafo,matrices);
    else
        transform_structure(id_trafo,matrices);
    end;
end;
% chains=length(model.structures{id_trafo});
% for k=1:chains,
%     sets=length(model.structures{id_trafo}(k));
%     for ks=1:sets,
%         coor0=model.structures{id_trafo}(k).xyz{ks};
%         [m,n]=size(coor0);
%         for ka=1:m,
%             coor0(ka,:)=affine_trafo_point(coor0(ka,:),matrices);
%         end;
%         model.structures{id_trafo}(k).xyz{ks}=coor0;
%     end;
% end;
hMain.graph_update=1;
set(handles.figure_magic,'Pointer','arrow');
delete(handles.figure_magic);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure_magic);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global matrices

global general

my_path=pwd;
cd(general.reports);

templates=get(handles.popupmenu_template,'String');
template=templates{get(handles.popupmenu_template,'Value')};
tobefitted=get(handles.popupmenu_to_fit,'String');
tofit=tobefitted{get(handles.popupmenu_to_fit,'Value')};
aligned=get(handles.checkbox_aligned_residues,'Value');
whole=get(handles.checkbox_whole_structure,'Value');

mode=0;
mode_name='all';
if get(handles.radiobutton_backbone,'Value'),
    mode=1;
    mode_name='backbone';
end;
if get(handles.radiobutton_CA,'Value'),
    mode=2;
    mode_name='CA';
end;
selected=get(handles.checkbox_selected,'Value');

fname0=sprintf('superpose_%s_onto_%s_by_%s.dat',tofit,template,mode_name);
[fname,pname] = uiputfile('.dat','Save transformation matrix',fname0);
if isequal(fname,0) || isequal(pname,0),
    add_msg_board('Saving of the transformation matrix canceled by user.');
    return;
end;
reset_user_paths(pname);
general.exp_files=pname;
if fname~=0,
    oname2=[pname fname];        
    set(handles.figure_magic,'Pointer','watch');
    matrices=find_magic_fit(handles,template,tofit,mode,selected,aligned,whole);
    if ~isempty(matrices),
        fid=fopen(oname2,'wt');
        fprintf(fid,'%%Affine transformation matrix to %s.\n',fname0(1:end-4));
        fprintf(fid,'%%Use function affine_trafo_point.m to transform coordinates.\n');
        for k=1:4,
            fprintf(fid,'%8.4f%8.4f%8.4f%8.4f\n',matrices(k,:));
        end;
        fclose(fid);
    end;
    set(handles.figure_magic,'Pointer','arrow');
end;

cd(my_path);

function [transmat,template,target]=find_magic_fit(handles,template,tofit,mode,selected,aligned,whole,sel_model)

global model

if exist('sel_model','var'),
    ensemble_mode = true;
else
    ensemble_mode = false;
end;

transmat=eye(4,4);

id_template=tag2id(template,model.structure_tags,model.structure_ids);
id_trafo=tag2id(tofit,model.structure_tags,model.structure_ids);

if selected || aligned || whole, % in all these cases arrays sel1 and sel2 of selected residues are generated
    if ~aligned && ~whole, % user selection
        [selections,msg]=resolve_address('*');
        [msel,nsel]=size(selections);
        poi1=0;
        poi2=0;
        sel1=zeros(msel,4);
        sel2=zeros(msel,4);
        if msel>0,
            for k=1:msel,
                cindices=selections(k,:);
                cindices=cindices(cindices>0);
                if length(cindices)==4,
                    if selections(k,1)==id_template,
                            poi1=poi1+1;
                            sel1(poi1,:)=cindices;
                    end;
                    if selections(k,1)==id_trafo,
                        poi2=poi2+1;
                        sel2(poi2,:)=cindices;
                    end;
                end;
            end;
        end;
    end;
    cid1=model.chain_ids(id_template);
    cid1=cid1{1};
    cid2=model.chain_ids(id_trafo);
    cid2=cid2{1};
    if whole, % try to align whole structures
        if length(cid1)~=length(cid2),
            add_msg_board('ERROR: Different number of chains in template and moving structure.');
            add_msg_board('Deactivate "whole structure".');
            return
        end;
        sel1=zeros(2000,4);
        sel2=sel1;
        psel=0;
        for k=1:length(cid1),
            seqs{1}=model.structures{id_template}(cid1(k)).sequence;
            seqs{2}=model.structures{id_trafo}(cid2(k)).sequence;
            sindices=[id_template,cid1(k);id_trafo,cid2(k)];
            [message,inname]=align_sequences(seqs,sindices,true);
            if message.error,
                add_msg_board('ERROR: MUSCLE sequence alignment failed.');
                add_msg_board(message.text);
                add_msg_board('Deactivate "whole structure".');
                return
            end;
            alignment=get_multiple_clustal(inname);
            [asel1,asel2]=select_aligned(alignment,sindices);
            [msel,nsel]=size(asel1);
            sel1(psel+1:psel+msel,:)=asel1;
            sel2(psel+1:psel+msel,:)=asel2;
            psel=psel+msel;
        end;
        sel1=sel1(1:psel,:);
        sel2=sel2(1:psel,:);
        poi1=psel;
        poi2=psel;
    elseif aligned, % we have to establish if there is only one chain per structure
                    % or if the same number of chains per structure is selected
        if length(cid1)==1 && length(cid2)==1,
            seqs{1}=model.structures{id_template}(cid1).sequence;
            seqs{2}=model.structures{id_trafo}(cid2).sequence;
            sindices=[id_template,cid1;id_trafo,cid2];
            [message,inname]=align_sequences(seqs,sindices,true);
            if message.error,
                add_msg_board('ERROR: MUSCLE sequence alignment failed.');
                add_msg_board(message.text);
                add_msg_board('Deactivate "whole structure".');
                return
            end;
            alignment=get_multiple_clustal(inname);
            [sel1,sel2]=select_aligned(alignment,sindices);
            [poi1,nsel]=size(sel1);
            poi2=poi1;
        else
            [selections,msg]=resolve_address('*');
            [msel,nsel]=size(selections);
            csel1=zeros(msel,2);
            csel2=zeros(msel,2);
            poi1=0;
            poi2=0;
            for k=1:msel,
                cindices=selections(k,:);
                cindices=cindices(cindices>0);
                if length(cindices)==2,
                    if selections(k,1)==id_template,
                        poi1=poi1+1;
                        csel1(poi1,:)=cindices;
                    end;
                    if selections(k,1)==id_trafo,
                        poi2=poi2+1;
                        csel2(poi2,:)=cindices;
                    end;
                end;
            end;
            csel1=csel1(1:poi1,:);
            csel2=csel2(1:poi2,:);
            if poi1~=poi2,
                add_msg_board('ERROR: Different number of chains selected in template and moving structure.');
                add_msg_board('Deactivate "aligned" to fit whole structure or change selection.');
                return
            end;
            sel1=zeros(2000,4);
            sel2=sel1;
            psel=0;
            for k=1:poi1,
                seqs{1}=model.structures{csel1(k,1)}(csel1(k,2)).sequence;
                seqs{2}=model.structures{csel2(k,1)}(csel2(k,2)).sequence;
                sindices=[csel1(k,:);csel2(k,:)];
                [message,inname]=align_sequences(seqs,sindices,true);
                if message.error,
                    add_msg_board('ERROR: MUSCLE sequence alignment failed.');
                    add_msg_board(message.text);
                    add_msg_board('Deactivate "whole structure".');
                    return
                end;
                alignment=get_multiple_clustal(inname);
                [asel1,asel2]=select_aligned(alignment,sindices);
                [msel,nsel]=size(asel1);
                sel1(psel+1:psel+msel,:)=asel1;
                sel2(psel+1:psel+msel,:)=asel2;
                psel=psel+msel;
            end;
            sel1=sel1(1:psel,:);
            sel2=sel2(1:psel,:);
            poi1=psel;
            poi2=psel;
        end;
    end;
    if poi1==0,
        add_msg_board('ERROR: No residues selected or aligned in template structure.');
        add_msg_board('Deactivate "only selected" or "aligned" to fit whole structure.');
        return
    end;
    if poi2==0,
        add_msg_board('ERROR: No residues selected or aligned in structure to be fitted.');
        add_msg_board('Deactivate "only selected" or "aligned" to fit whole structure.');
        return
    end;
    sel1=sel1(1:poi1,:);
    sel2=sel2(1:poi2,:);
    coor0=zeros(20000,3);
    coor1=coor0;
    poi=0;
    test_corr=get(handles.checkbox_check_correspondence,'Value');
    if ~test_corr && poi1~=poi2,
        add_msg_board('ERROR: Different number of selected residues and no correspondence check.');
        return
    end;
    for k1=1:poi1,
        cindices1=sel1(k1,:);
        [stag1,ctag1,modelnum1,resnum1,icode1]=mk_address_parts(cindices1);
        cindices2=[];
        if test_corr,
            for k2=1:poi2,
                tindices=sel2(k2,:);
                [stag2,ctag2,modelnum2,resnum2,icode2]=mk_address_parts(tindices);
                if strcmpi(ctag1,ctag2) && resnum1==resnum2 && strcmpi(icode1,icode2),
                    cindices2=sel2(k2,:);
                end;
            end;
        else
            cindices2=sel2(k1,:);
        end;
        if isempty(cindices2),
        else
            adr1=mk_address(cindices1);
            switch mode
                case 0
                    adr1=sprintf('%s.:',adr1);
                case 1
                    if model.structures{cindices1(1)}(cindices1(2)).seqtype == 1
                        adr1=sprintf('%s.N,CA,C,O',adr1);
                    elseif model.structures{cindices1(1)}(cindices1(2)).seqtype == 2
                        adr1=sprintf('%s.P,O5'',C5'',C4'',C3'',O3''',adr1);
                    end
                case 2
                    if model.structures{cindices1(1)}(cindices1(2)).seqtype == 1
                        adr1=sprintf('%s.CA',adr1);
                    elseif model.structures{cindices1(1)}(cindices1(2)).seqtype == 2
                        adr1=sprintf('%s.C4''',adr1);
                        % cheat for nucleotide backbone and sugar
                        % adr1=sprintf('%s.P,O5'',C5'',C4'',C3'',O3'',C1'',C2'',O2'',O4''',adr1);
                  end
            end;
            [indices,message]=resolve_address(adr1);

            [m,n]=size(indices);
            for k=1:m,
                cindices0=indices(k,:);
                cindices0=cindices0(1:5);
                address=mk_address(cindices0);
                postpoi=findstr(address,'.');
                post=address(postpoi:end);
                [stag2,ctag2,modelnum2,resnum2,icode2]=mk_address_parts(cindices2);
                if ensemble_mode,
                    modelnum2 = sel_model;
                end;
                adr2=sprintf('[%s](%s){%i}%i%s%s',stag2,ctag2,modelnum2,resnum2,icode2,post);
                [cindices,msg]=resolve_address(adr2);
                if msg.error==0 && ~isempty(cindices),
                    poi=poi+1;
                    [m2,n2]=size(cindices);
                    [msg,cc0]=get_atom(cindices0,'coor');
                    [msg,cc1]=get_atom(cindices,'coor');
                    if ~isempty(cc0) && ~isempty(cc1),
                        coor0(poi,:)=cc0;
                        coor1(poi,:)=cc1;
                    else
                        poi=poi-1;
                    end;
                end;
            end;
        end;
    end;
    coor0=coor0(1:poi,:);
    coor1=coor1(1:poi,:);
    add_msg_board(sprintf('%i selected atoms are fitted',poi));
else
    adr1=sprintf('[%s](:).',template);
    switch mode
        case 0
            adr1=sprintf('%s:',adr1);
        case 1
            adr1=sprintf('%sN,CA,C,O,P,O5'',C5'',C4'',C3'',O3''',adr1);
        case 2
            adr1=sprintf('%sCA,C4''',adr1);
            % cheat for nucleotdie backbone and sugar
            % adr1=sprintf('%sCA,P,O5'',C5'',C4'',C3'',O3'',C1'',C2'',O2'',O4''',adr1);
    end;

    [indices,message]=resolve_address(adr1);

    [m,n]=size(indices);
    
    coor0=zeros(m,3);
    coor1=coor0;

    poi=0;
    poi2=0;
    for k=1:m,
        cindices0=indices(k,:);
        cindices0=cindices0(1:5);
        address=mk_address(cindices0);
        prepoi=findstr(address,'[');
        pre=address(1:prepoi);
        postpoi=findstr(address,']');
        post=address(postpoi:end);
        adr2=[pre tofit post];
        [cindices,msg]=resolve_address(adr2);
        if msg.error==0 && ~isempty(cindices),
            if ensemble_mode,
                cindices(3) = sel_model;
            end;
            poi=poi+1;
            [m2,n2]=size(cindices);
            if m2>1,
                poi2=poi2+1;
            end;
            [msg,cc0]=get_atom(cindices0,'coor');
            [msg,cc1]=get_atom(cindices,'coor');
            if ~isempty(cc0) && ~isempty(cc1),
                coor0(poi,:)=cc0;
                coor1(poi,:)=cc1;
            else
                poi=poi-1;
            end;
        end;
    end;
    coor0=coor0(1:poi,:);
    coor1=coor1(1:poi,:);
    add_msg_board(sprintf('%i atoms were selected',m));
    add_msg_board(sprintf('of which %i atoms were found in 2nd structure',poi));
    if poi2>0,
        add_msg_board(sprintf('with %i of them being ambiguous',poi2));
    end;
end;


[m,n]=size(coor0);
coor2b=[];
if m>0 && n>0,
    [rmsd,coor2b,transmat]=rmsd_superimpose(coor0,coor1);
    add_msg_board(sprintf('RMSD of both structures is: %4.2f Å',rmsd));
else
    add_msg_board('Warning! No transformation performed as selections do not match.');
    transmat=[];
end;

% store Calpha coordinates for move display
if mode==2,
    template=coor2b;
    target=coor0;
else
    template=[];
    target=[];
end;

% --- Executes on button press in checkbox_selected.
function checkbox_selected_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_selected


% --- Executes on button press in checkbox_check_correspondence.
function checkbox_check_correspondence_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_check_correspondence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_check_correspondence


% --- Executes on button press in checkbox_aligned_residues.
function checkbox_aligned_residues_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_aligned_residues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_aligned_residues


% --- Executes on button press in checkbox_whole_structure.
function checkbox_whole_structure_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_whole_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_whole_structure

function [sel1,sel2]=select_aligned(alignment,sindices)
% select only residues that are aligned

global model

seq1=alignment(1).sequence;
seq2=alignment(2).sequence;
sel1=zeros(length(seq1),4);
sel2=sel1;
poi1=0;
poi2=0;
pois=0;
for k=1:length(seq1),
    match1=false;
    match2=false;
    if char(seq1(k))~='-',
        poi1=poi1+1;
        match1=true;
    end;
    if char(seq2(k))~='-',
        poi2=poi2+1;
        match2=true;
    end;
    if match1 && match2
        tag1=sprintf('%i',poi1);
        id1=tag2id(tag1,model.structures{sindices(1,1)}(sindices(1,2)).residues{1}.residue_tags);
        tag2=sprintf('%i',poi2);
        id2=tag2id(tag2,model.structures{sindices(2,1)}(sindices(2,2)).residues{1}.residue_tags);
        if ~isempty(id1) && ~isempty(id2),
            pois=pois+1;
            sel1(pois,:)=[sindices(1,:) 1 id1];
            sel2(pois,:)=[sindices(2,:) 1 id2];
        end;
    end;
end;
sel1=sel1(1:pois,:);
sel2=sel2(1:pois,:);


% --- Executes on button press in checkbox_motion.
function checkbox_motion_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_motion

if ~get(handles.radiobutton_CA,'Value') && get(hObject,'Value'),
    set(hObject,'Value',0);
    msgbox('Select superposition of C-alpha atoms to allow for motion display','Motion display impossible.','warn');
end;
guidata(hObject,handles);


% --- Executes on button press in checkbox_ensemble_mode.
function checkbox_ensemble_mode_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ensemble_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ensemble_mode
