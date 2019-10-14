function varargout = rigi_flex(varargin)
% RIGI_FLEX MATLAB code for rigi_flex.fig
%      RIGI_FLEX, by itself, creates a new RIGI_FLEX or raises the existing
%      singleton*.
%
%      H = RIGI_FLEX returns the handle to a new RIGI_FLEX or the handle to
%      the existing singleton*.
%
%      RIGI_FLEX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIGI_FLEX.M with the given input arguments.
%
%      RIGI_FLEX('Property','Value',...) creates a new RIGI_FLEX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rigi_flex_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rigi_flex_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rigi_flex

% Last Modified by GUIDE v2.5 04-Oct-2017 08:19:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rigi_flex_OpeningFcn, ...
                   'gui_OutputFcn',  @rigi_flex_OutputFcn, ...
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


% --- Executes just before rigi_flex is made visible.
function rigi_flex_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rigi_flex (see VARARGIN)

global MMM_icon
global hMain
global model

% Choose default command line output for rigi_flex
handles.output = hObject;

j = get(hObject,'javaframe');    
j.setFigureIcon(javax.swing.ImageIcon(im2java(MMM_icon)));  %create a java image and set the figure icon

load helpicon
set(handles.pushbutton_help,'CData',cdata);

if ~exist('model','var') || ~isfield(model,'current_structure'),
    set(hMain.popupmenu_view,'Value',1);
    % initialize display
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
    adr = '[none]';
else
    adr=mk_address(model.current_structure);
end

set(handles.fig_rigi_flex,'Name',sprintf('RigiFlex structure model based on rigid bodies %s',adr));

handles.max_time = 12;
set(handles.edit_max_time,'String',sprintf('%5.2f',handles.max_time));

set(handles.text_time_left,'String','Idle.');
set(handles.text_time_left,'ForegroundColor',[0,127,0]/256);

handles.ensemble = 20;
handles.ensemble_user = handles.ensemble;
set(handles.edit_ensemble,'String',sprintf('%i',handles.ensemble_user));

handles.p_model = erf(1/sqrt(2));
handles.p_model_user = handles.p_model;
set(handles.edit_p_model,'String',sprintf('%5.3f',handles.p_model_user));
handles.edit_dr_std.String = sprintf('%5.3f',sqrt(2)*erfinv(handles.p_model));

handles.max_trials = 500000000;
set(handles.edit_max_trials,'String',sprintf('%i',handles.max_trials));

handles.SANS_threshold = 15;
set(handles.edit_SANS_chi2,'String',sprintf('%5.1f',handles.SANS_threshold));

handles.SAXS_threshold = 9;
set(handles.edit_SAXS_chi2,'String',sprintf('%5.1f',handles.SAXS_threshold));

handles.xlink_threshold = 30;
set(handles.edit_xlink_threshold,'String',sprintf('%5.1f',handles.xlink_threshold));

handles.xlink_percentage = 0;
set(handles.edit_xlink_percentage,'String',sprintf('%5.1f',handles.xlink_percentage));

handles.diagnostics.success = 1;

handles.pushbutton_save.Enable =  'off';
handles.radiobutton_DEER.Enable =  'off';

handles.copy = false;
handles.curr_model = 1;
handles.curr_restraint = 1;
handles.progress = 0;

handles.plot_position = handles.axes_multi_plot.OuterPosition;


hMain.auxiliary=[hMain.auxiliary hObject];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rigi_flex wait for user response (see UIRESUME)
% uiwait(handles.fig_rigi_flex);


% --- Outputs from this function are returned to the command line.
function varargout = rigi_flex_OutputFcn(hObject, eventdata, handles) 
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

entry=strcat(help_files,'domain_ensemble.html');
webcall(entry,'-helpbrowser');


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global hMain
% global model
% 
% hMain.domain_ensemble_plot=false;
% model.current_structure=handles.current_structure;
% update_current_structure;
delete(handles.fig_rigi_flex);

% --- Executes on button press in pushbutton_restraints.
function pushbutton_restraints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global model
global general


my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load restraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Restraint loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    [~,name] = fileparts(fname);
    handles.save_path = pname;
    handles.save_name = name;
    
    hfig=gcf;
    set(hfig,'Pointer','watch');
    [restraints,failed] = rd_restraints_rigiflex(fullfile(pname,fname));
    if failed
        add_msg_board('ERROR: Reading of restraint file failed');
        return
    end
    handles.restraint_file = fullfile(pname,fname);
    if ~isempty(restraints.models)
        handles.ensemble = restraints.models;
    else
        restraints.models = handles.ensemble;
    end
    if isempty(restraints.ensemble)
        restraints.ensemble = handles.ensemble;
    end
    
    if ~isempty(restraints.max_time)
        handles.max_time = restraints.max_time;
    else
        restraints.max_time = handles.max_time;
    end
    
    adr=mk_address(model.current_structure);
    set(handles.fig_rigi_flex,'Name',sprintf('RigiFlex structure model based on rigid bodies %s',adr));
    
    [handles,restraints] = analyze_exhaustive(handles,restraints);
     
end

if ~isempty(restraints.p_model)
    handles.p_model = restraints.p_model;
else
    restraints.p_model = handles.p_model;
end

set(handles.edit_ensemble,'String',sprintf('%i',handles.ensemble));
set(handles.edit_p_model,'String',sprintf('%5.3f',handles.p_model));
handles.edit_dr_std.String = sprintf('%5.3f',sqrt(2)*erfinv(handles.p_model));
handles.text_SANS_curves.String = sprintf('(%i curves)',length(restraints.SANS));
handles.text_SAXS_curves.String = sprintf('(%i curves)',length(restraints.SAXS));
handles.text_xlinks.String = sprintf('(%i restraints)',length(restraints.xlinks));
handles.text_xlinks_required.String = sprintf('(%i restraints)',ceil(handles.xlink_percentage*length(restraints.xlinks)/100));
handles.edit_max_time.String = sprintf('%5.2f',handles.max_time);

handles.restraints = restraints;
if isfield(restraints,'solutions') && ~isempty(restraints.solutions)
    handles.pushbutton_run.String = 'Run RigiFlex';
else
    handles.pushbutton_run.String = 'Run Rigi';
end
set(handles.pushbutton_run,'Enable','on');
set(handles.pushbutton_restraints,'Enable','off');
set(hfig,'Pointer','arrow');
cd(my_path);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global general
global model

if isfield(handles.restraints,'solutions') && ~isempty(handles.restraints.solutions)
    last_step = 2;
    handles.progress = 0;
else
    last_step = 0;
end

for runstep = 0:last_step   
    
    switch handles.progress

        case 0 % Rigi
            options.granularity = 10000;
            options.max_trials = handles.max_trials;
            options.max_time = handles.max_time;
            options.deterministic = handles.checkbox_std_seed.Value;
            options.SANS_threshold = handles.SANS_threshold;
            options.SAXS_threshold = handles.SAXS_threshold;
            options.xlink_threshold = handles.xlink_threshold;
            options.xlink_percentage = handles.xlink_percentage;
            options.display_SANS_fit = handles.radiobutton_SANS_fit.Value;
            options.display_SANS_chi2 = handles.radiobutton_SANS_chi2.Value;
            options.display_SAXS_fit = handles.radiobutton_SAXS_fit.Value;
            options.display_xlinks = handles.radiobutton_crosslinks.Value;
            options.exhaustive = handles.checkbox_exhaustive.Value;

            [filename, pathname] = uiputfile([handles.restraints.newID '.pdb'], 'Save final model as PDB');
            if isequal(filename,0) || isequal(pathname,0)
                add_msg_board('RigiFlex modelling cancelled by user');
                return
            else
                reset_user_paths(pathname);
                general.pdb_files=pathname;
                fname = fullfile(pathname, filename);
                options.fname = fname;
            end

            set(gcf,'Pointer','watch');
            handles.text_time_left.String = 'Started.';
            handles.text_time_left.ForegroundColor = [0,0,180]/256;
            drawnow
            snum0 = model.current_structure;
            handles.template_snum = snum0;
            diagnostics = rigi_flex_engine(handles.restraints,options,handles);
            handles.diagnostics = diagnostics;
            fprintf(1,'Maximum runtime         : %4.1f h\n',options.max_time);
            fprintf(1,'Used runtime            : %4.1f h\n',handles.diagnostics.runtime/3600);
            fprintf(1,'Maximum number of trials: %i\n',options.max_trials);
            fprintf(1,'Used number of trials   : %i\n',handles.diagnostics.trials);
            if diagnostics.success
                add_msg_board('Rigi step successfully completed.');
                handles.text_time_per_model.String = sprintf('%12.1f',diagnostics.runtime/diagnostics.success);
                [pathstr,basname] = fileparts(options.fname);
                handles.basname = basname;
                report_name = fullfile(pathstr,strcat(basname,'_rigi_report.txt'));
                result_name = fullfile(pathstr,strcat(basname,'_rigi_diagnostics.mat'));
                handles.options = options;
                if handles.restraints.rigi_distributions
                    distributions = mk_report_distributions(report_name,handles,options);
                    handles.distributions = distributions;
                    handles.radiobutton_DEER.Enable =  'on';
                else
                    handles.distributions = [];
                    distributions = [];
                end
                restraints = handles.restraints;
                save(result_name,'diagnostics','options','restraints','distributions','report_name');
                handles.pushbutton_save.Enable =  'on';
                handles.pushbutton_run.String = 'Run Flex';
                if isfield(handles.restraints,'RNA') && isfield(handles.restraints.RNA,'bind') && ~isempty(handles.restraints.RNA.bind)
                    handles.pushbutton_run.String = 'RNA links';
                    handles.progress = 1;
                else
                    handles.pushbutton_run.String = 'Run Flex';
                    handles.progress = 2;
                    if isfield(handles.restraints,'flex_time')
                        handles.max_time = handles.restraints.flex_time;
                        handles.edit_max_time.String = sprintf('%5.2f',handles.max_time);
                    end
                end
            else
                add_msg_board('RigiFlex modeling run failed.');
                handles.text_time_per_model.String = 'n.a.';
            end


            handles.text_success.String = sprintf('%i',diagnostics.success);
            handles.text_time_left.String = 'Completed.';
            handles.text_time_left.ForegroundColor = [0,127,0]/256;

            update_plot(handles);

            set(gcf,'Pointer','arrow');

            guidata(hObject,handles);

        case 1 % RNA Flex
            if isempty(handles.restraints.pflex)
                only_FlexRNA = true;
            else
                only_FlexRNA = false;
            end
            [pathstr,basname] = fileparts(handles.options.fname);
            solutions_given = true;
            soln_name = strcat(basname,'_solutions.dat');
            soln_name = fullfile(pathstr,soln_name);
            solutions = load(soln_name);
            rba_solutions = round(solutions);
            stag = mk_address_parts(handles.diagnostics.snum);
            [pathstr,basname] = fileparts(handles.options.fname);
            report_name = fullfile(pathstr,strcat(basname,'_RNA_link_report.txt'));
            linkable_name = fullfile(pathstr,strcat(basname,'_RNA_linkable.dat'));
            [solutions,trafo,stemlibs,dvecs] = fit_stemloop_combinations(handles.restraints,handles.template_snum,handles.diagnostics.snum,report_name);
            % assign the stemloops to rigid bodies, needed for RBA adjustment
            rb_assign = zeros(1,length(stemlibs));
            rna_chains = zeros(1,length(stemlibs));
            for klib = 1:length(stemlibs)
                slindi = resolve_address(sprintf('%s%s',stag,stemlibs{klib}.chaintag));
                rna_chains(klib) = slindi(2);
                for kr = 1:length(handles.restraints.rb)
                    for kcc = 1:length(handles.restraints.rb(kr).chains)
                        if slindi(2) == handles.restraints.rb(kr).chains(kcc)
                            rb_assign(klib) = kr;
                        end
                    end
                end
            end
            report_fid = fopen(report_name,'at');
            % handles.diagnostics.snum = 2; % ### only for debugging
            handles = set_progress_interface(handles);
            handles.pushbutton_run.String = 'RNA links';
            set(gcf,'Pointer','watch');
            handles.text_time_left.String = 'Started.';
            handles.text_time_left.ForegroundColor = [0,0,180]/256;
            rba = [handles.diagnostics.snum,1];
            secdefs = process_rna_domain_restraints(handles.restraints,rba);
            RNA_linkers = zeros(handles.diagnostics.success,length(secdefs));
            RNA_linker_time = zeros(handles.diagnostics.success,length(secdefs));
            if isfield(handles.restraints,'RNA') && isfield(handles.restraints.RNA,'bind')
                num_ch = length(model.structures{handles.diagnostics.snum});
                initflag = true;
                chain_indices = zeros(1,length(handles.restraints.stemlibs));
                for ksl = 1:length(handles.restraints.stemlibs)
                    ind = resolve_address(sprintf('[%s]%s',stag,handles.restraints.stemlibs{ksl}.chaintag));
                    chain_indices(ksl) = ind(2);
                end
                success_vec = zeros(1,handles.diagnostics.success);
                initialize_model = true;
                for km = 1:handles.diagnostics.success % loop over rigid-body models
                    csoln = solutions(solutions(:,1)==km,:);
                    stretch = csoln(:,2);
                    [stretch,strpoi] = sort(stretch);
                    csoln = csoln(strpoi,:);
                    success = false;
                    kcombi = 0;
                    snum = handles.diagnostics.snum;
                    while ~isempty(stretch) && ~success && kcombi < length(stretch) && kcombi < 1
                        kcombi = kcombi + 1;
                        ccombi = csoln(kcombi,3:end);
                        % make adjustment transformation matrices
                        adj_transmats = cell(1,length(handles.restraints.rb));
                        for kr = 1:length(handles.restraints.rb)
                            adj_transmats{kr} = eye(4);
                        end
                        for kr = 2:length(rb_assign)
                            if rb_assign(kr) > 1
                                baspoi6 = (rb_assign(kr)-2)*6;
                                dtrans = 10*dvecs(km,baspoi6+1:baspoi6+3);
                                % fprintf(1,'Translation(%i): (%3.1f, %3.1f, %3.1f) Å\n',kr,dtrans);
                                deuler = dvecs(km,baspoi6+4:baspoi6+6);
                                % fprintf(1,'Rotation(%i)   : (%3.1f, %3.1f, %3.1f)°\n',kr,180*deuler/pi);
                                adj_transmats{rb_assign(kr)} = transrot2affine(dtrans,deuler);
                            end
                        end
                        % apply adjustment of RBA
                        for kr = 1:length(handles.restraints.rb)                       
                            for kc = 1:length(handles.restraints.rb(kr).chains)
                                if min(abs(rna_chains-handles.restraints.rb(kr).chains(kc))) ~= 0
                                    % adr_rep = mk_address([snum handles.restraints.rb(kr).chains(kc) km]);
              %                      fprintf(1,'Replacing coordinates of %s\n',adr_rep);
                                    model.structures{snum}(handles.restraints.rb(kr).chains(kc)).xyz{km} = ...
                                        affine_coor_set(model.structures{snum}(handles.restraints.rb(kr).chains(kc)).xyz{km},adj_transmats{kr});
                                end
                            end
                        end
                        % replace binding motifs by stemloops from library if stemloop
                        % libraries are tested
                        for klib = 1:length(stemlibs)
                            chain_coor = stemlibs{klib}.chains{ccombi(klib)}.xyz{1};
                            chain_coor = affine_coor_set(chain_coor,trafo{km,klib});
                            if initialize_model
                                fields = fieldnames(model.structures{snum}(chain_indices(klib)));
                                for kfield = 1:length(fields)
                                    if isfield(stemlibs{klib}.chains{ccombi(klib)},fields{kfield})
                                        model.structures{snum}(chain_indices(klib)).(fields{kfield}) = stemlibs{klib}.chains{ccombi(klib)}.(fields{kfield});
                                    end
                                end
                            end
                            model.structures{snum}(chain_indices(klib)).xyz{km} = chain_coor;
                            model.structures{snum}(chain_indices(klib)).residues{km} = stemlibs{klib}.chains{ccombi(klib)}.residues{1};
                            model.structures{snum}(chain_indices(klib)).Bfactor{km} = stemlibs{klib}.chains{ccombi(klib)}.Bfactor{1};
                            model.structures{snum}(chain_indices(klib)).Btensor{km} = stemlibs{klib}.chains{ccombi(klib)}.Btensor{1};
                            model.structures{snum}(chain_indices(klib)).atoms{km} = stemlibs{klib}.chains{ccombi(klib)}.atoms{1};
                            model.structures{snum}(chain_indices(klib)).isotopes = stemlibs{klib}.chains{ccombi(klib)}.isotopes;
                        end
                        initialize_model = false;
                        env_sel = zeros(num_ch,2);
                        for kb = 1:length(handles.restraints.RNA.bind)
                            adra = sprintf('[%s]%s',stag,handles.restraints.RNA.bind(kb).anchora);
                            inda = resolve_address(adra);
                            adre = sprintf('[%s]%s',stag,handles.restraints.RNA.bind(kb).anchora);
                            inde = resolve_address(adre);
                            env_sel(inda(2),1) = inda(4);
                            env_sel(inda(2),2) = inde(4);
                        end

                        handles.text_success.String = sprintf('%i',km);
                        rba = [handles.diagnostics.snum,km];
                        model.current_structure = handles.diagnostics.snum;
                        [secdefs,RNA] = process_rna_domain_restraints(handles.restraints,rba);
                        parts = length(handles.restraints.RNA_tags) + length(secdefs);
                        pindices = ones(parts,3);
                        for sec = 0:length(secdefs)
                            indpart = resolve_address(sprintf('%s%s',stag,handles.restraints.RNA_tags{sec+1}));
                            pindices(2*sec+1,1:2) = indpart;
                            pindices(2*sec+1,3) = km;
                        end
                        % Make environment coordinates
                        maxatoms = 50000;
                        environ = zeros(maxatoms,3);
                        poi = 0;
                        for kc = 1:num_ch
                            cmind = [handles.diagnostics.snum,kc,km];
                            if env_sel(kc,1) == 0
                                [~,xyz] = get_chain_model(cmind,'xyz_heavy');
                                [m,~] = size(xyz);
                                environ(poi+1:poi+m,:) = xyz;
                                poi = poi+m;
                                stag = mk_address_parts(cmind);
                                for kr = env_sel(kc,1)+1:env_sel(kc,2)-1
                                    [~,xyz] = get_residue([cmind,kr],'xyz_heavy');
                                    [m,~] = size(xyz);
                                    environ(poi+1:poi+m,:) = xyz;
                                    poi = poi+m;
                                end
                            end
                        end
                        environ = environ(1:poi,:);
                        ensemble = 1; % for this application of rna_flex_engine, ensemble must be 1
                        maxtime = handles.restraints.RNA.maxtime;
                        full_connect = true;
                        for sec = 1:length(secdefs)
                            handles.text_dmg_fail.String = sprintf('%i',sec);
                            handles.text_auxiliary_fail.String = sprintf('%5.2f',0);
                            handles.text_core_fail.String = sprintf('%5.2f',0);
                            drawnow
                            diagnostics = rna_flex_engine(handles,secdefs(sec),RNA,environ,ensemble,maxtime,[km,sec]);
                            fprintf(report_fid,'Rigid body arrangement %i, section %i\n',km,sec);
                            if ~isempty(diagnostics.snum)
                                fprintf(1,'RNA model stored in structure %i\n',diagnostics.snum);
                                RNA_linkers(km,sec) = diagnostics.snum;
                                RNA_linker_time(km,sec) = diagnostics.time_per_model;
                                pindices(2*sec,1) = diagnostics.snum;
                                fprintf(report_fid,'Successful modelling (structure %i) with %i models and %5.1f s/model\n\n',diagnostics.snum,diagnostics.success,diagnostics.time_per_model);
                            else
                                fprintf(1,'RNA modelling failed for section %i in RBA %i\n',sec,km);
                                fprintf(report_fid,'Modelling of RNA for RBA %i failed. Remaining RNA sections (if any) are skipped.\n\n',km);
                                full_connect = false;
                                break
                            end
                        end
                        if full_connect
                            fprintf(report_fid,'Modelling of RNA for RBA %i succeeded.\n\n',km);
                            nindices = [handles.diagnostics.snum,num_ch+1,km];
                            chain_tag = handles.restraints.RNA.chain_id;
                            poia = strfind(chain_tag,'(');
                            if isempty(poia)
                                poia = 0;
                            end
                            poie = strfind(chain_tag,')');
                            if isempty(poie)
                                poie = length(chain_tag)+1;
                            end
                            chain_tag = chain_tag(poia+1:poie-1);
                            stitch_chain(pindices,nindices,chain_tag,initflag);
                            initflag = false;
                            success_vec(km) = 1;
                            success = true;
                            if solutions_given
                                link_fid = fopen(linkable_name,'at');
                                fprintf(link_fid,'%8i%6i\n',rba_solutions(km,1),rba_solutions(km,2));
                                fclose(link_fid);
                            end
                            if only_FlexRNA % conformer should be saved, if it is already complete
                                new_cindices = zeros(1,num_ch+1);
                                ncp = 0;
                                for kc = 1:num_ch
                                    [~,ctag] = mk_address_parts([handles.diagnostics.snum,kc]);
                                    isRNA = false;
                                    for krna = 1:length(handles.restraints.RNA_tags)
                                        if strcmpi(['(' ctag ')'],handles.restraints.RNA_tags{krna})
                                            isRNA = true;
                                        end
                                    end
                                    if ~isRNA
                                        ncp = ncp + 1;
                                        new_cindices(ncp) = kc;
                                    end
                                end
                                new_cindices(ncp+1) = num_ch + 1;
                                new_cindices = new_cindices(1:ncp+1);
                                handles.new_cindices = new_cindices;
                                if isfield(model,'selected')
                                    model = rmfield(model,'selected');
                                end
                                spoi = 0;
                                for kc = new_cindices
                                    spoi = spoi + 1;
                                    model.selected{spoi} = [snum kc km];
                                end
                                fmname = sprintf('%s_m%i_flex.pdb',handles.basname,km);
                                message = wr_pdb_selected(fmname,handles.restraints.newID);
                                if message.error
                                    if interactive
                                        add_msg_board(sprintf('Warning: Model %s_m%i could not be automatically saved. %s\n',handles.basname,km,message.text));
                                    end
                                end
                            end
                        end
                    end
                end
            end
            handles.RNA_link_success = success_vec;
            handles = mk_RNA_report_distributions(report_fid,handles);
            fclose(report_fid);
            handles.text_time_left.String = 'Completed.';
            handles.text_time_left.ForegroundColor = [0,127,0]/256;
            set(gcf,'Pointer','arrow');
            handles.progress = 2;
            if only_FlexRNA
                handles.progress = 3;
            end
            handles.pushbutton_run.String = 'Flex Peptide';
            guidata(hObject,handles);

        case 2 % Flex
            if isfield(handles.restraints,'protein')
                pctag = handles.restraints.protein(2);
            else
                pctag = 'X';
            end
            snum = handles.diagnostics.snum;
            num_ch = length(model.structures{handles.diagnostics.snum});
            RNA = num_ch;
            to_be_saved = zeros(handles.diagnostics.success,2);
            if isfield(handles.restraints,'solutions') && ~isempty(handles.restraints.solutions)
                list_name = sprintf('%s_model_list.dat',handles.restraints.solutions);
            else
            [pathstr,basname] = fileparts(handles.options.fname);
                list_name = sprintf('%s_model_list.dat',basname);
                list_name = fullfile(pathstr,list_name);
            end

            options.min_approach = 1.2;
            options.deterministic = handles.checkbox_std_seed.Value;
            options.max_trials = handles.max_trials;
            options.max_time = handles.max_time;
            options.fname_bas = handles.basname;
            %     report_name = fullfile(pathstr,strcat(handles.basname,'_report.txt'));
            handles.report_name = strcat(handles.basname,'_flex_report.txt');
            fidr = fopen(handles.report_name,'wt');
            fprintf(fidr,'--- Generation of flexible peptide termini and linkers for rigid-body model %s ---\n',handles.basname);
            fclose(fidr);
            all_flex_models = zeros(handles.diagnostics.success,length(handles.restraints.pflex));
            all_flex_model_times = zeros(handles.diagnostics.success,length(handles.restraints.pflex));
            model.current_structure = handles.diagnostics.snum;
            adr = mk_address(model.current_structure);
            if isfield(handles,'RNA_link_success')
                comp_mask = handles.RNA_link_success;
            else
                comp_mask = ones(1,handles.diagnostics.success);
            end
            flex_success = zeros(handles.diagnostics.success,length(handles.restraints.pflex));
            stag = mk_address_parts(handles.diagnostics.snum);
            connected_chains = zeros(length(handles.restraints.pflex),2);
            for kp = 1:length(handles.restraints.pflex)
                Nchain = handles.restraints.pflex(kp).Nanchor(1:3);
                Nadr = sprintf('[%s]%s',stag,Nchain);
                Nind = resolve_address(Nadr);
                Cchain = handles.restraints.pflex(kp).Canchor(1:3);
                Cadr = sprintf('[%s]%s',stag,Cchain);
                Cind = resolve_address(Cadr);
                connected_chains(kp,1) = Nind(2);
                connected_chains(kp,2) = Cind(2);
            end
            for km = 1:handles.diagnostics.success % loop over rigid-body models
                options.max_trials = handles.max_trials;
                options.rm = km;
                options.template = sprintf('%s(:){%i}',adr,km);
                if ~comp_mask(km)
                    add_msg_board(sprintf('Warning: No RNA connection in rigid-body arrangement %i',km));
                    add_msg_board('Skipping computation of flexible peptide sections in this rigid-body arrangement');
                    options.max_trials = -1;
                else
                    for kp = 1:length(handles.restraints.pflex) % loop over flexible peptide segments
                        options.fd = kp;
                        handles.uipanel_runtime.Title = sprintf('Flexible peptide section %i in rigid-body arrangment %i',kp,km);
                        handles = set_progress_interface(handles);
                        set(gcf,'Pointer','watch');
                        handles.text_time_left.String = 'Started.';
                        handles.text_time_left.ForegroundColor = [0,0,180]/256;
                        drawnow
                        model.current_structure = handles.diagnostics.snum; % make the result of rigid-body arrangement the current structure
                        [restrain,restraints,monitor,cancelled,number,number_monitor] = process_domain_restraints(handles.restraints.pflex(kp),km);
                        options.max_time = handles.restraints.pflex(kp).time;
                        options.max_trials = 1e6;
                        set(handles.edit_max_time,'String',sprintf('%5.2f',options.max_time));
                        drawnow
                        if cancelled
                            add_msg_board(sprintf('ERROR: Restraint processing failed for flexible domain %i in model %i',kp,kp));
                            options.max_trials = -1;
                        end
                        options.n_restraints = number;
                        options.n_monitor = number_monitor;
                        options.monitor = monitor;
                        fidr = fopen(handles.report_name,'at');
                        fprintf(fidr,'\n### Flexible peptide section %i in rigid-body arrangment %i ###\n',kp,km);
                        fclose(fidr);
                        flex_diagnostics = flex_engine(restraints,restrain,options,handles);
                        if flex_diagnostics.success == 0
                            add_msg_board(sprintf('Warning: No models found for flexible domain %i in rigid-body arrangement %i',kp,km));
                            fprintf(fidr,'Modelling of flexible domain %i in rigid-body arrangement %i failed after %i s\n',kp,km,flex_diagnostics.runtime);
                            add_msg_board('Skipping remaining flexible domains in this RBA');
                            break
                            options.max_trials = -1;
                        else
                            flex_success(km,kp) = 1;
                        end
                        all_flex_models(km,kp) = flex_diagnostics.snum;
                        all_flex_model_times(km,kp) = flex_diagnostics.time_per_model;
                        handles.flex_diagnostics{km,kp} = flex_diagnostics;
                        handles.text_time_left.String = 'Completed.';
                        handles.text_time_left.ForegroundColor = [0,127,0]/256;
                        drawnow
                        set(gcf,'Pointer','arrow');
                    end
                end
                linked_model = 0;
                chains = zeros(1,num_ch+1);
                ncp = 0;
                for kc = 1:num_ch
                    [~,ctag] = mk_address_parts([handles.diagnostics.snum,kc]);
                    isRNA = false;
                    for krna = 1:length(handles.restraints.RNA_tags)
                        if strcmpi(['(' ctag ')'],handles.restraints.RNA_tags{krna})
                            isRNA = true;
                        end
                    end
                    if ~isRNA
                        ncp = ncp + 1;
                        chains(ncp) = kc;
                    end
                end
                save_it = prod(flex_success(km,:));
                for kp = 1:length(handles.restraints.pflex)
                    if flex_success(km,kp)
                        for kc = 1:num_ch
                            if chains(kc) == connected_chains(kp,1) || chains(kc) == connected_chains(kp,2)
                                chains(kc) = 0;
                            end
                        end
                        linked_model = all_flex_models(km,kp); 
                    end
                end
                if save_it
                    cch = 0;
                    found = true;
                    indices = cell(1,100);
                    while found
                        cch = cch + 1;
                        ctag = id2tag(cch,handles.restraints.peptide_tags{1});
                        if isempty(ctag)
                            found = false;
                            cch = cch - 1;
                        else
                            [cind,msg] = resolve_address(sprintf('[%s](%s)',stag,ctag));
                            if ~isempty(cind) && ~msg.error % this is a rigid-body peptide chain
                                indices{cch} = [snum cind(2) km];
                            else
                                flexnum = str2double(ctag);
                                secstruct = all_flex_models(km,flexnum);
                                indices{cch} = [secstruct 1 1];
                            end
                        end
                    end
                    indices = indices(1:cch);
                    chain = combine_chains(indices,pctag);
                    info = set_info(chain,linked_model);
                    [~,linked_model] = add_pdb(chain,info);
                    to_be_saved(km,1) = linked_model;
                    to_be_saved(km,2) = RNA;
                    if isfield(model,'selected')
                        model = rmfield(model,'selected');
                    end
                    model.selected{1} = [linked_model 1 1];
                    if isfield(handles.restraints,'RNA') && isfield(handles.restraints.RNA,'bind') && ~isempty(handles.restraints.RNA.bind)
                        model.selected{2} = [handles.diagnostics.snum RNA km];
                    end
                    fmname = sprintf('%s_m%i_flex.pdb',handles.basname,km);
                    wr_pdb_selected(fmname,handles.restraints.newID,[],true);
                    fid_list = fopen(list_name,'at');
                    fprintf(fid_list,'%s\n',fmname);
                    fclose(fid_list);
                end
            end
            handles.all_flex_models = all_flex_models;
            handles.all_flex_model_times = all_flex_model_times;
            handles.flex_success = flex_success;
            handles.flex_saved = to_be_saved; 
            handles.progress = 3;
            add_msg_board('Flex step completed.');
            if isfield(handles.restraints,'build_time')
                handles.max_time = handles.restraints.build_time;
                handles.edit_max_time.String = sprintf('%5.2f',handles.max_time);
            end
            handles.pushbutton_run.String = 'Run assembler';
            guidata(hObject,handles);
            return
        otherwise
            add_msg_board('Full RigiFlex run completed');
    end
end

function edit_max_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_time as text
%        str2double(get(hObject,'String')) returns contents of edit_max_time as a double

[v,handles]=edit_update_MMM(handles,hObject,0.05,1000,2,'%5.2f',0);
handles.max_time = v;
if isfield(handles,'restraints')
    [handles,restraints] = analyze_exhaustive(handles,handles.restraints);
    handles.restraints = restraints;
else
    handles.text_exhaustive_resolution.String = sprintf('n.a.');
    handles.edit_max_trials.String = sprintf('%i',handles.max_trials);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_max_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ensemble as text
%        str2double(get(hObject,'String')) returns contents of edit_ensemble as a double

[v,handles]=edit_update_MMM(handles,hObject,1,10000,20,'%i',1);
handles.ensemble = v;
handles.restraints.models = v;
handles.restraints.ensemble = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_ensemble_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_trials_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_trials as text
%        str2double(get(hObject,'String')) returns contents of edit_max_trials as a double

[v,handles]=edit_update_MMM(handles,hObject,1000,100000000000,5000000000,'%i',1);
handles.max_trials = v;
if isfield(handles,'restraints')
    [handles,restraints] = analyze_exhaustive(handles,handles.restraints);
    handles.restraints = restraints;
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_max_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_p_model_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p_model as text
%        str2double(get(hObject,'String')) returns contents of edit_p_model as a double

[v,handles]=edit_update_MMM(handles,hObject,0.05,0.95,0.5,'%5.3f',0);
handles.p_model = v;
handles.edit_dr_std.String = sprintf('%5.3f',sqrt(2)*erfinv(v));
handles.restraints.p_model = v;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_p_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_adaptive.
function checkbox_adaptive_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_adaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_adaptive


% --- Executes on button press in checkbox_save_cdir.
function checkbox_save_cdir_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_save_cdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_save_cdir


function edit_min_approach_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_approach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_approach as text
%        str2double(get(hObject,'String')) returns contents of edit_min_approach as a double

[v,handles]=edit_update_MMM(handles,hObject,0.5,2.5,2.0,'%5.2f',0);
handles.min_approach = v;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_min_approach_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_approach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function distributions = mk_report_distributions(report_name,handles,options)

sadr = mk_address(handles.diagnostics.snum);

fid_report = fopen(report_name,'wt');

fprintf(fid_report,'--- DEER restraint report ---\n');

fprintf(fid_report,'\nRestraint file used: %s\n',handles.restraint_file);

fprintf(fid_report,'\nMaximum number of models: %i\n',handles.restraints.models);
fprintf(fid_report,'Probability threshold   : %5.2f\n',handles.restraints.p_model);
if ~isempty(handles.restraints.SANS)
    fprintf(fid_report,'SANS fit chi^2 threshold: %5.2f\n',options.SANS_threshold);
end
if ~isempty(handles.restraints.SAXS)
    fprintf(fid_report,'SAXS fit chi^2 threshold: %5.2f\n',options.SAXS_threshold);
end
if ~isempty(handles.restraints.xlinks)
    fprintf(fid_report,'Crosslink threshold     : %5.1f Å\n',options.xlink_threshold);
    fprintf(fid_report,'%4.1f%% crosslinks must match.\n',options.xlink_percentage);
end
fprintf(fid_report,'Probability threshold   : %i\n',handles.restraints.p_model);
fprintf(fid_report,'Maximum runtime         : %4.1f h\n',options.max_time);
fprintf(fid_report,'Used runtime            : %4.1f h\n',handles.diagnostics.runtime/3600);
fprintf(fid_report,'Maximum number of trials: %i\n',options.max_trials);
fprintf(fid_report,'Used number of trials   : %i\n',handles.diagnostics.trials);

if options.deterministic
    fprintf(fid_report,'Standard seed used for pseudo-random number generator\n');
else
    fprintf(fid_report,'Shuffled seed used for pseudo-random number generator\n');
end
fprintf(fid_report,'Parallelization granularity: %i\n',options.granularity);
fprintf(fid_report,'\n%i models were generated.\n',handles.diagnostics.success);

% if ~isempty(handles.restraints.SANS)
%     fprintf(fid_report,'--- SANS fit results ---\n');
%     for km = 1:handles.diagnostics.success
%         fprintf(fid_report,'Model %i, SANS chi^2 = %4.2f',km,handles.diagnostics.final_chi2_SANS(km));
%         if length(handles.restraints.SANS) > 1
%             fprintf(fid_report,'(');
%             for ks = 1:length(handles.restraints.SANS)
%                 fprintf(fid_report,'%4.2f',handles.diagnostics.chi_SANS(ks,km));
%                 if ks < length(handles.restraints.SANS)
%                     fprintf(fid_report,', ');
%                 end
%             end
%             fprintf(fid_report,')');
%         end
%         fprintf(fid_report,'\n');
%     end
% end
poi = 0;
restraints = handles.restraints;
is_core_restraint = zeros(1,length(restraints.DEER));
all_overlaps = zeros(1,length(restraints.DEER));
known = zeros(1,9*length(restraints.rb)*(length(restraints.rb)-1)/2);
for kr1 = 1:length(restraints.rb)-1
    for kr2 = kr1+1:length(restraints.rb)
        for kp1 = 1:3
            for kp2 = 1:3
                poi = poi + 1;
                distributions(poi).type = 'Core restraint';
                distributions(poi).label1 = restraints.rb(kr1).label{kp1};
                distributions(poi).label2 = restraints.rb(kr2).label{kp2};
                distributions(poi).ind1 = restraints.rb(kr1).indices(kp1,:);
                distributions(poi).ind2 = restraints.rb(kr2).indices(kp2,:);
                ind1 = restraints.rb(kr1).indices(kp1,:);
                ind1(1) = handles.diagnostics.snum;
                ind2 = restraints.rb(kr2).indices(kp2,:);
                ind2(1) = handles.diagnostics.snum;
                distr = zeros(handles.diagnostics.success,301); % assume that get_distribution returns a distribution with 300 points
                sum_distr = zeros(1,301);
                sum_var = 0;
                for km = 1:handles.diagnostics.success
                    ind1(3) = km;
                    ind2(3) = km;
%                     adr1 = mk_address(ind1);
%                     adr2 = mk_address(ind2);
                    NO_pos1 = get_NO_pos(ind1,distributions(poi).label1,298);
                    NO_pos2 = get_NO_pos(ind2,distributions(poi).label2,298);
                    [rax,sim_distr] = get_distribution(NO_pos1,NO_pos2,0.05);
                    distr(km,:) = sim_distr;
                    rmean = sum(rax.*sim_distr);
                    cvar = sum(sim_distr.*(rax-rmean).^2);
                    sum_var = sum_var + cvar*handles.diagnostics.probabilities(km);
                    sum_distr= sum_distr + sim_distr*handles.diagnostics.probabilities(km);
                end
                sum_distr = sum_distr/sum(handles.diagnostics.probabilities);
                mean_var = sum_var/sum(handles.diagnostics.probabilities);
                distributions(poi).all_distr = distr;
                distributions(poi).sum_distr = sum_distr;
                distributions(poi).r_req = -1;
                distributions(poi).sigr_req = -1;
                [~,ctag1,~,resnum1] = mk_address_parts(ind1);
                [~,ctag2,~,resnum2] = mk_address_parts(ind2);
                distributions(poi).adr1 = sprintf('(%s)%i',ctag1,resnum1);
                distributions(poi).adr2 = sprintf('(%s)%i',ctag2,resnum2);
                distributions(poi).rax = rax;
                fprintf(fid_report,'\nRigid body %i, ref. point %i (%s)%i to rigid body %i, ref. point %i (%s)%i\n',kr1,kp1,ctag1,resnum1,kr2,kp2,ctag2,resnum2);
                meanr = 10*sum(sum_distr.*rax);
                sigr = sqrt(sum(sum_distr.*(10*rax-meanr).^2));
                distributions(poi).meanr = meanr;
                distributions(poi).sigr = sigr;
                fprintf(fid_report,'<r> = %4.1f Å, sigr = %4.1f Å, Rigidity index for ensemble: %5.2f\n',meanr,sigr,sqrt(mean_var/sigr^2));
                match = false;
                for kr = 1:length(restraints.DEER)
                    if ~isempty(restraints.DEER(kr).indices1) && ~isempty(restraints.DEER(kr).indices2)
                        rind1 = restraints.DEER(kr).indices1([2,4]);
                        rind2 = restraints.DEER(kr).indices2([2,4]);
                        if sum(abs(rind1-ind1([2,4]))) == 0 && sum(abs(rind2-ind2([2,4]))) == 0
                            match = true;
                            is_core_restraint(kr) = 1;
                            r_restr = restraints.DEER(kr).r;
                            sigr_restr = restraints.DEER(kr).sigr;
                        end
                        if sum(abs(rind2-ind1([2,4]))) == 0 && sum(abs(rind1-ind2([2,4]))) == 0
                            match = true;
                            is_core_restraint(kr) = 1;
                            r_restr = restraints.DEER(kr).r;
                            sigr_restr = restraints.DEER(kr).sigr;
                        end
                    end
                end
                if match
                    distributions(poi).r_req = r_restr;
                    distributions(poi).sigr_req = sigr_restr;
                    sim_distr = (rax-r_restr/10).^2/(2*0.01*sigr_restr^2);
                    sim_distr = exp(-sim_distr);
                    sim_distr = sim_distr/sum(sim_distr);
                    known(poi) = 1;
                    overlap = sum(min([sum_distr;sim_distr]));
                    fprintf(fid_report,' for restraint <r> = %4.1f Å, sigr = %4.1f Å (overlap: %6.4f).\n',r_restr,sigr_restr,overlap);
                    fprintf(fid_report,'Rigidity index to restraint: %5.2f\n',sqrt(mean_var/sigr_restr^2));
                else
                    fprintf(fid_report,' (unrestrained).\n');
                end
                mean_overlap = 0;
                for km1 = 1:handles.diagnostics.success-1
                    for km2 = km1+1:handles.diagnostics.success
                        mean_overlap = mean_overlap + sum(min([distr(km1,:);distr(km2,:)])); 
                    end
                end
                mean_overlap = 2*mean_overlap/((handles.diagnostics.success-1)*handles.diagnostics.success);
                all_overlaps(poi) = mean_overlap;
                distributions(poi).overlap = mean_overlap;
                fprintf(fid_report,'Mean overlap of distributions for the ensemble: %6.4f.\n',mean_overlap);
            end
        end
    end
end

poi_core = poi;

for kr = 1:length(restraints.DEER)
    if ~is_core_restraint(kr)
        ind1 = restraints.DEER(kr).indices1;
        if isempty(ind1) && ~isempty(restraints.DEER(kr).adr1)
            ind1 = resolve_address(sprintf('%s%s',sadr,restraints.DEER(kr).adr1));
        end
        ind2 = restraints.DEER(kr).indices2;
        if isempty(ind2) && ~isempty(restraints.DEER(kr).adr2)
            ind2 = resolve_address(sprintf('%s%s',sadr,restraints.DEER(kr).adr2));
        end
        if ~isempty(ind1) && ~isempty(ind2)
            poi = poi + 1;
            distributions(poi).type = 'Auxiliary restraint';
            psep = strfind(restraints.DEER(kr).label,'|');
            if isempty(psep)
                distributions(poi).label1 = restraints.DEER(kr).label;
                distributions(poi).label2 = restraints.DEER(kr).label;
            else
                distributions(poi).label1 = restraints.DEER(kr).label(1:psep-1);
                distributions(poi).label2 = restraints.DEER(kr).label(psep+1:end);
            end
            distributions(poi).ind1 = restraints.DEER(kr).indices1;
            distributions(poi).ind2 = restraints.DEER(kr).indices2;
            ind1(1) = handles.diagnostics.snum;
            ind2(1) = handles.diagnostics.snum;
            distr = zeros(handles.diagnostics.success,301); % assume that get_distribution returns a distribution with 300 points
            sum_distr = zeros(1,301);
            sum_var = 0;
            for km = 1:handles.diagnostics.success
                ind1(3) = km;
                ind2(3) = km;
                NO_pos1 = get_NO_pos(ind1,distributions(poi).label1,298);
                NO_pos2 = get_NO_pos(ind2,distributions(poi).label2,298);
                [rax,sim_distr] = get_distribution(NO_pos1,NO_pos2,0.05);
                rmean = sum(rax.*sim_distr);
                cvar = sum(sim_distr.*(rax-rmean).^2);
                distr(km,:) = sim_distr;
                sum_var= sum_var + cvar*handles.diagnostics.probabilities(km);
                sum_distr= sum_distr + sim_distr*handles.diagnostics.probabilities(km);
            end
            sum_distr = sum_distr/sum(handles.diagnostics.probabilities);
            mean_var = sum_var/sum(handles.diagnostics.probabilities);
            distributions(poi).all_distr = distr;
            distributions(poi).sum_distr = sum_distr;
            distributions(poi).r_req = restraints.DEER(kr).r;
            distributions(poi).sigr_req = restraints.DEER(kr).sigr;
            [~,ctag1,~,resnum1] = mk_address_parts(ind1);
            [~,ctag2,~,resnum2] = mk_address_parts(ind2);
            distributions(poi).adr1 = sprintf('(%s)%i',ctag1,resnum1);
            distributions(poi).adr2 = sprintf('(%s)%i',ctag2,resnum2);
            distributions(poi).rax = rax;
            fprintf(fid_report,'\nAuxiliary restraint (%s)%i to (%s)%i\n',ctag1,resnum1,ctag2,resnum2);
            meanr = 10*sum(sum_distr.*rax);
            sigr = sqrt(sum(sum_distr.*(10*rax-meanr).^2));
            distributions(poi).meanr = meanr;
            distributions(poi).sigr = sigr;
            fprintf(fid_report,'<r> = %4.1f Å, sigr = %4.1f Å for restraint <r> = %4.1f Å, sigr = %4.1f Å\n',meanr,sigr,restraints.DEER(kr).r,restraints.DEER(kr).sigr);
            fprintf(fid_report,'Rigidity index for ensemble: %5.2f; to restraint: %5.2f\n',sqrt(mean_var/sigr^2),sqrt(mean_var/restraints.DEER(kr).sigr^2));
            sim_distr = (rax-restraints.DEER(kr).r/10).^2/(2*0.01*restraints.DEER(kr).sigr^2);
            sim_distr = exp(-sim_distr);
            sim_distr = sim_distr/sum(sim_distr);
            overlap = sum(min([sum_distr;sim_distr]));
            fprintf(fid_report,' (overlap: %6.4f).\n',overlap);
            mean_overlap = 0;
            for km1 = 1:handles.diagnostics.success-1
                for km2 = km1+1:handles.diagnostics.success
                    mean_overlap = mean_overlap + sum(min([distr(km1,:);distr(km2,:)]));
                end
            end
            mean_overlap = 2*mean_overlap/((handles.diagnostics.success-1)*handles.diagnostics.success);
            distributions(poi).overlap = mean_overlap;
            fprintf(fid_report,'Mean overlap of distributions for the ensemble: %6.4f.\n',mean_overlap);
        end
    end
end

if sum(known) < poi_core
    fprintf(fid_report,'\n--- Potential further core restraints ordered by expected restraint strength ---\n\n'); 
    [sort_overlaps,spoi] = sort(all_overlaps);
    for k = 1:length(spoi)
        if spoi(k) <= length(known) && ~known(spoi(k))
            fprintf(fid_report,'%s-%s with ensemble overlap %6.4f\n',distributions(spoi(k)).adr1,distributions(spoi(k)).adr2,sort_overlaps(k));
        end
    end
end

fprintf(fid_report,'\n--- Restraint file ---\n\n');

fid_restraints = fopen(handles.restraint_file,'r');

if fid_restraints == -1
    fprintf(fid_report,'### Warning: Restraint file no longer available. ###\n');
    fclose(fid_report);
    return;
end

while 1
    tline = fgetl(fid_restraints);
    if ~ischar(tline), break, end
    fprintf(fid_report,'%s\n',tline);
end

fclose(fid_restraints);
fclose(fid_report);

function NO_pos = get_NO_pos(indices,label,T)

global model
global label_defs
global hMain

if strcmpi(label,'CA'),
    adr = sprintf('%s.CA',mk_address(indices));
    [~,xyz] = get_object(adr,'coor');
    NO_pos = [xyz 1];
    return
end;

NO_pos = [];
if isfield(model,'sites'),
    for k0=1:length(model.sites),
        for k1=1:length(model.sites{k0}),
            for k=1:length(model.sites{k0}(k1).residue)
                if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0,
                    id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                    if strcmpi(label,label_defs.residues(id).short_name)  || strcmpi(label,label_defs.residues(id).tc),
                        if T == model.sites{k0}(k1).residue(k).T,
                            NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

if isempty(NO_pos),
    adr = mk_address(indices);
    command=sprintf('rotamers %s %s %i',adr,label,T);
    hMain.store_undo=false;
    hMain.dynamic_rotamers=false;
    cmd(hMain,command);
end;

for k0=1:length(model.sites),
    for k1=1:length(model.sites{k0}),
        for k=1:length(model.sites{k0}(k1).residue),
            if sum(abs(indices-model.sites{k0}(k1).residue(k).indices)) == 0,
                id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                if strcmpi(label,label_defs.residues(id).short_name) || strcmpi(label,label_defs.residues(id).tc),
                    if T == model.sites{k0}(k1).residue(k).T,
                        NO_pos=model.sites{k0}(k1).residue(k).NOpos;
                    end;
                end;
            end;
        end;
    end;
end;


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.copy = true;
update_plot(handles);

% --- Executes on button press in pushbutton_model_forth.
function pushbutton_model_forth_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_model_forth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'diagnostics')
    return
end;
if handles.curr_model < handles.diagnostics.success,
    handles.curr_model = handles.curr_model + 1;
    handles.edit_current_model.String = sprintf('%i',handles.curr_model);
    update_plot(handles);
else
    add_msg_board('Warning: No more models.');
end;
guidata(hObject,handles);

% --- Executes on button press in pushbutton_model_back.
function pushbutton_model_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_model_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.curr_model > 1
    handles.curr_model = handles.curr_model - 1;
    handles.edit_current_model.String = sprintf('%i',handles.curr_model);
    update_plot(handles);
else
    add_msg_board('Warning: Model number must be positive.');
end
guidata(hObject,handles);

function update_plot(handles)

if handles.progress < 1 
    return
end

if handles.diagnostics.success < 1 
    return
end

if handles.copy
    figure;
    copy_mode = true;
    handles.copy = false;
else
    axes(handles.axes_multi_plot);
    copy_mode = false;
    cla;
end

cm = handles.curr_model;
cr = handles.curr_restraint;

handles.text_model_probability.String = sprintf('%5.3f',handles.diagnostics.probabilities(cm));
if ~isempty(handles.distributions)
    handles.text_restraint_type.String = handles.distributions(cr).type;
end
if ~copy_mode
    set(gca,'OuterPosition',handles.plot_position);
    set(gca,'XTickLabelRotation','default');
    set(gca,'XTickLabelMode','auto');
    cla
end
if handles.radiobutton_crosslinks.Value
    hold on;
    xlink_fulfill = handles.diagnostics.xlink_fulfill;
    [ml,nm] = size(xlink_fulfill);
    for km = 1:nm
        plot(xlink_fulfill(:,km),'.');
    end
    if ml > 0
        plot([1,ml],[handles.xlink_threshold,handles.xlink_threshold],'Color',[0.25,0.25,0.25]);
        xlink_labels = cell(1,ml);
        for kl = 1:ml
            [~,ctag,~,resnum]=mk_address_parts(handles.restraints.xlinks(kl).indices1);
            xltag = sprintf('%s%i',ctag,resnum);
            [~,ctag,~,resnum]=mk_address_parts(handles.restraints.xlinks(kl).indices2);
            xltag = sprintf('%s/%s%i',xltag,ctag,resnum);
            xlink_labels{kl} = xltag;
        end
    end
    if ~copy_mode
        set(gca,'OuterPosition',[36.8360 2.5 116.7300 17.8]);
    end
    set(gca,'XTick',1:ml);
    set(gca,'XTickLabelRotation',90);
    set(gca,'XTickLabel',xlink_labels);
    title('Crosslink fulfillment');
end
if handles.radiobutton_DEER.Value
    distr = handles.distributions(cr).all_distr;
    info = handles.distributions(cr);
    hold on;
    plot(info.rax,info.sum_distr,'k');
    plot(info.rax, distr(cm,:),'Color',[0,0,0.75]); 
    my_title = sprintf('DEER %s-%s <r> = %4.1f Å, sr = %4.1f Å Overlap: %5.3f',info.adr1,info.adr2,info.meanr,info.sigr,info.overlap);
    if handles.distributions(cr).r_req >= 0 && handles.distributions(cr).sigr_req >= 0
        sim_distr = (info.rax-info.r_req/10).^2/(2*0.01*info.sigr_req^2);
        sim_distr = exp(-sim_distr);
        sim_distr = sim_distr/sum(sim_distr);
        plot(info.rax,sim_distr,'Color',[0.75,0,0]);
        overlap = sum(min([info.sum_distr;sim_distr]));
        my_title = sprintf('%s, Fulfillment: %5.3f',my_title,overlap);
    end
    title(my_title);
end

guidata(handles.pushbutton_run,handles);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

diagnostics = handles.diagnostics;
fname = save_data(diagnostics,[handles.restraints.newID '_diagnostics.mat'],'diagnostics');
add_msg_board(sprintf('Diagnostics information saved to file %s',fname));
fname = save_data(handles.distributions,[handles.restraints.newID '_distributions.mat'],'distance distribution');
add_msg_board(sprintf('Distance distributions saved to file %s',fname));
guidata(handles.pushbutton_save,handles);

function pfname = save_data(data,suggestion,comment)

global general

my_path=pwd;
cd(general.DEER_files);

[fname,pname]=uiputfile('*.mat',sprintf('Save %s information',comment),suggestion);
if isequal(fname,0) || isequal(pname,0),
    add_msg_board(sprintf('Saving of %s information on the modeling run was aborted.',comment));
    return;
end;
reset_user_paths(pname);
general.DEER_files=pname;

pfname = fullfile(pname,fname);
save(pfname,'data');

cd(my_path);


% --- Executes on button press in checkbox_std_seed.
function checkbox_std_seed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_std_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_std_seed

function [overlap,shift] = get_overlap(rax,distr1,distr2)

distr1 = distr1/sum(distr1);
distr2 = distr2/sum(distr2);
shift = sum(rax.*distr1) - sum(rax.*distr2); 
% fprintf(1,'Sim. mean distance: %5.2f\n',sum(rax.*distr1));
% fprintf(1,'Restraint mean distance: %5.2f\n',sum(rax.*distr2));
overlap = sum(min([distr1;distr2]));
    



function edit_SANS_chi2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SANS_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SANS_chi2 as text
%        str2double(get(hObject,'String')) returns contents of edit_SANS_chi2 as a double

[v,handles]=edit_update_MMM(handles,hObject,0.1,1000,15,'%5.1f',0);
handles.SANS_threshold = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_SANS_chi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SANS_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SAXS_chi2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SAXS_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SAXS_chi2 as text
%        str2double(get(hObject,'String')) returns contents of edit_SAXS_chi2 as a double

[v,handles]=edit_update_MMM(handles,hObject,0.1,1000,9,'%5.1f',0);
handles.SAXS_threshold = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_SAXS_chi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SAXS_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xlink_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xlink_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xlink_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_xlink_threshold as a double

[v,handles]=edit_update_MMM(handles,hObject,3,1000,30,'%5.1f',0);
handles.xlink_threshold = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_xlink_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xlink_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_current_model_Callback(hObject, eventdata, handles)
% hObject    handle to edit_current_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_current_model as text
%        str2double(get(hObject,'String')) returns contents of edit_current_model as a double

[v,handles]=edit_update_MMM(handles,hObject,1,handles.diagnostics.success,1,'%i',1);
handles.curr_model = v;
update_plot(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_current_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_current_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_restraint_forth.
function pushbutton_restraint_forth_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restraint_forth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'diagnostics')
    return
end;
if handles.curr_restraint < length(handles.distributions),
    handles.curr_restraint = handles.curr_restraint + 1;
    handles.edit_restraint.String = sprintf('%i',handles.curr_restraint);
    update_plot(handles);
else
    add_msg_board('Warning: No more restraints.');
end;
guidata(hObject,handles);

% --- Executes on button press in pushbutton_restraint_back.
function pushbutton_restraint_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_restraint_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'diagnostics')
    return
end;
if handles.curr_restraint > 1,
    handles.curr_restraint = handles.curr_restraint - 1;
    handles.edit_restraint.String = sprintf('%i',handles.curr_restraint);
    update_plot(handles);
else
    add_msg_board('Warning: Restraint number must be positive.');
end;
guidata(hObject,handles);


function edit_restraint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_restraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_restraint as text
%        str2double(get(hObject,'String')) returns contents of edit_restraint as a double

[v,handles]=edit_update_MMM(handles,hObject,1,length(handles.distributions),1,'%i',1);
handles.curr_restraint= v;
update_plot(handles);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit_restraint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_restraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_SANS_chi2.
function radiobutton_SANS_chi2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_SANS_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_plot(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_SANS_fit.
function radiobutton_SANS_fit_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_SANS_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_plot(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_SAXS_fit.
function radiobutton_SAXS_fit_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_SAXS_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_plot(handles);


% --- Executes on button press in radiobutton_SANS_chi2.
function radiobutton_SANS_chi2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SANS_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SANS_chi2

update_plot(handles);


% --- Executes on button press in radiobutton_SANS_fit.
function radiobutton_SANS_fit_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SANS_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SANS_fit

update_plot(handles);


% --- Executes on button press in radiobutton_SAXS_fit.
function radiobutton_SAXS_fit_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_SAXS_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_SAXS_fit

update_plot(handles);


% --- Executes on button press in radiobutton_crosslinks.
function radiobutton_crosslinks_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_crosslinks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_crosslinks

update_plot(handles);


% --- Executes on button press in radiobutton_DEER.
function radiobutton_DEER_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_DEER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_DEER

update_plot(handles);



function edit_xlink_percentage_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xlink_percentage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xlink_percentage as text
%        str2double(get(hObject,'String')) returns contents of edit_xlink_percentage as a double
[v,handles]=edit_update_MMM(handles,hObject,0,100,80,'%5.1f',0);
handles.xlink_percentage = v;
handles.text_xlinks_required.String = sprintf('(%i restraints)',ceil(v*length(handles.restraints.xlinks)/100));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_xlink_percentage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xlink_percentage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dr_std_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dr_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dr_std as text
%        str2double(get(hObject,'String')) returns contents of edit_dr_std as a double

[v,handles]=edit_update_MMM(handles,hObject,0.05,3,1,'%5.3f',0);
handles.p_model = erf(v/sqrt(2));
handles.edit_p_model.String = sprintf('%5.3f',handles.p_model);
handles.restraints.p_model = handles.p_model;
guidata(hObject,handles);

function [restrain,restraints,monitor,cancelled,number,number_monitor] = process_domain_restraints(restraints,modnum)
% Processes restraints for a single flexible domain

global rotamer_libraries
global model
global residue_defs

cancelled=false;

restrain = [];
monitor = [];
number = 0;
number_monitor = 0;

if ~isfield(restraints,'sequence')
    add_msg_board('ERROR: Sequence specification is missing.'); 
    restrain=[];
    cancelled = true;
    return;
end

for k = 1:length(restraints.sequence)
    restrain(k).secondary = 0;
    restrain(k).label = [];
    restrain(k).cis = 0;
    restrain(k).r_beacon = [];
    restrain(k).r_intern = [];
    restrain(k).oligomer = [];
    restrain(k).depth = [];
end

monitor = restrain;

if ~isfield(restraints,'start') || ~isfield(restraints,'end')
    add_msg_board('ERROR: Domain is not specified in restraint file.'); 
    restrain=[];
    cancelled = true;
    return;
else
    res1 = restraints.start;
    rese = restraints.end;
    if isnan(res1) || isnan(rese)
        add_msg_board('ERROR: Residue numbers of domain could not be recognized.'); 
        restrain=[];
        cancelled = true;
        return;
    end
end

restraints.res1 = res1;
restraints.rese = rese;

if isfield(restraints,'DEER')
    poi = 0;
    llist = cell(0);
    label = cell(0);
    for k = 1:length(restraints.DEER)
        [indices,message] = resolve_address(restraints.DEER(k).adr1);
        if message.error ~= 2 && message.error ~= 13 % this site exists and is thus a beacon
            if message.error
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr1)); 
                address = mk_address(indices,1);
                add_msg_board(address);
                restrain=[];
                cancelled = true;
                return;
            end
            poi = poi + 1;
            restraints.DEER(k).type1 = 1;
            indices(3) = modnum;
            restraints.DEER(k).indices1 = indices;
            llist{poi} = mk_address(indices);
            label{poi} = restraints.DEER(k).label1;
        else % this site is in the loop to be modelled
            for kr = 1:length(rotamer_libraries)
                if strcmpi(rotamer_libraries(kr).label,restraints.DEER(k).label1) || strcmpi(rotamer_libraries(kr).tc,restraints.DEER(k).label1)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel1 = NO;
                    restraints.DEER(k).type1 = 0;
                end
            end
        end
        [indices,message]=resolve_address(restraints.DEER(k).adr2);
        if message.error ~= 2 && message.error ~= 13 % this site exists and is thus a beacon
            if message.error
                add_msg_board(sprintf('ERROR: Residue address %s is invalid.',restraints.DEER(k).adr2)); 
                address = mk_address(indices,1);
                add_msg_board(address);
                restrain=[];
                cancelled = true;
                return;
            end
            poi = poi + 1;
            restraints.DEER(k).type2 = 1;
            indices(3) = modnum;
            restraints.DEER(k).indices2 = indices;
            llist{poi} = mk_address(indices);
            label{poi} = restraints.DEER(k).label2;
        else % this site is in the loop to be modelled
            for kr = 1:length(rotamer_libraries)
                if strcmpi(rotamer_libraries(kr).label,restraints.DEER(k).label2) || strcmpi(rotamer_libraries(kr).tc,restraints.DEER(k).label2)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                    restraints.DEER(k).NO_rel2 = NO;
                    restraints.DEER(k).type2 = 0;
                end
            end
        end
    end
    labels = get_labels(llist,label);
    for k = 1:length(restraints.DEER)
        clabel1 = restraints.DEER(k).label1;
        clabel2 = restraints.DEER(k).label2;
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 0 % internal restraint
            resa = correct_section_address(restraints.DEER(k).adr1);
            resb = correct_section_address(restraints.DEER(k).adr2);
            if isempty(resa) || isempty(resb)
                add_msg_board(sprintf('ERROR: Residue address %s or %s is invalid.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO1 = restraints.DEER(k).NO_rel1;
            NO2 = restraints.DEER(k).NO_rel2;
            if resa < resb
                exch = resa; resa = resb; resb = exch;
                exch = NO1; NO1 = NO2; NO2 = exch;
            end
            if restraints.DEER(k).r ~= 0
                [restrain,number] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,clabel1,clabel2);
            else
                [monitor,number_monitor] = mk_internal_restraint(monitor,NO1,NO2,resa,resb,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,clabel1,clabel2);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 ==0 % beacon restraint, first residue beacon
            indices = restraints.DEER(k).indices1;
            adr1 = mk_address(indices);
            for kl = 1:length(labels)
                if sum(abs(indices-labels(kl).indices)) == 0
                    xyz_beacon = labels(kl).xyz;
                end
            end
            res_loop = correct_section_address(restraints.DEER(k).adr2);
            if isempty(res_loop)
                add_msg_board(sprintf('ERROR: Residue address %s inside domain is invalid.',restraints.DEER(k).adr2)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO = restraints.DEER(k).NO_rel2;
            if restraints.DEER(k).r ~= 0
                [restrain,number] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,clabel1,clabel2,indices,adr1);
            else
                [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,clabel1,clabel2,indices,adr1);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 0 && restraints.DEER(k).type2 == 1 % beacon restraint, second residue beacon
            indices = restraints.DEER(k).indices2;
            adr2 = mk_address(indices);
            for kl = 1:length(labels)
                if sum(abs(indices-labels(kl).indices)) == 0
                    xyz_beacon = labels(kl).xyz;
                end
            end
            res_loop = correct_section_address(restraints.DEER(k).adr1);
            if isempty(res_loop)
                add_msg_board(sprintf('ERROR: Residue address %s inside domain is invalid.',restraints.DEER(k).adr1)); 
                restrain=[];
                cancelled = true;
                return;
            end
            NO = restraints.DEER(k).NO_rel1;
            if restraints.DEER(k).r ~= 0
               [restrain,number] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number,clabel1,clabel2,indices,adr2);
            else
               [monitor,number_monitor] = mk_beacon_restraint(monitor,NO,res_loop,xyz_beacon,restraints.DEER(k).r,restraints.DEER(k).sigr,res1,number_monitor,clabel1,clabel2,indices,adr2);
                number_monitor = number_monitor + 1;
            end
        end
        if restraints.DEER(k).type1 == 1 && restraints.DEER(k).type2 == 1 % nonsensical restraint inside defined structure
            add_msg_board(sprintf('Warning: Restraint between residues %s and %s inside rigid bodies will be ignored.',restraints.DEER(k).adr1,restraints.DEER(k).adr2)); 
        end
    end
end

if isfield(restraints,'depth')
    for k = 1:length(restraints.depth)
        clabel = restraints.depth(k).label;
        if strcmpi(restraints.depth(k).label,'CA')
            NO = [];
        else
            for kr = 1:length(rotamer_libraries)
                if strcmpi(rotamer_libraries(kr).label,restraints.depth(k).label)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end
            end
        end
        res = restraints.depth(k).num;
        if isempty(res) || isnan(res)
            add_msg_board(sprintf('ERROR: Residue number %i inside domain is invalid.',res)); 
            restrain=[];
            cancelled = true;
            return;
        end
        if restraints.depth(k).z ~= 0
            [restrain,number] = mk_depth_restraint(restrain,NO,res,restraints.depth(k).z,restraints.depth(k).sigz,res1,number,clabel);
        else
            [monitor,number_monitor] = mk_depth_restraint(monitor,NO,res,restraints.depth(k).z,restraints.depth(k).sigz,res1,number_monitor,clabel);
            number_monitor = number_monitor + 1;
        end
    end
end

if isfield(restraints,'oligomer')
    for k = 1:length(restraints.oligomer)
        clabel = restraints.oligomer(k).label;
        if strcmpi(restraints.oligomer(k).label,'CA')
            NO = [];
        else
            for kr = 1:length(rotamer_libraries)
                if strcmpi(rotamer_libraries(kr).label,restraints.oligomer(k).label)
                    Tvec = rotamer_libraries(kr).T;
                    [mi,pT] = min(abs(Tvec-298));
                    if mi > eps
                        add_msg_board('Warning: No library exactly fits the specified labelling temperature.');
                    end
                    libname = id2tag(pT,rotamer_libraries(kr).files);
                    NO = get_relative_label(libname);
                end
            end
        end
        res = restraints.oligomer(k).num;
        if isempty(res) || isnan(res)
            add_msg_board(sprintf('ERROR: Residue number %i is invalid in oligomer restraint.',res)); 
            restrain=[];
            cancelled = true;
            return;
        end
        if restraints.oligomer(k).r ~= 0
            [restrain,number] = mk_oligomer_restraint(restrain,NO,res,restraints.oligomer(k).mult,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number,clabel);
        else
            [monitor,number_monitor] = mk_oligomer_restraint(monitor,NO,res,restraints.oligomer(k).mult,restraints.oligomer(k).r,restraints.oligomer(k).sigr,res1,number_monitor,clabel);
            number_monitor = number_monitor + 1;
        end
    end
end

if isfield(restraints,'cispeptides')
    for k = 1:length(restraints.cispeptides)
        kr = restraints.cispeptides(k) - res1 + 1;
        restrain(kr).cis = 1;
    end
end

if isfield(restraints,'aprop')
    for k = 1:length(restraints.aprop)
        res = restraints.aprop(k).num;
        if isempty(res) || isnan(res)
            add_msg_board(sprintf('ERROR: Residue number %i in domain is invalid in APROP.',restraints.aprop(k).num)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).aprop = restraints.aprop(k).prop;
    end
end

if isfield(restraints,'bprop')
    for k = 1:length(restraints.bprop)
        res = restraints.bprop(k).num;
        if isempty(res) || isnan(res)
            add_msg_board(sprintf('ERROR: Residue number %i is invalid in BPROP.',restraints.bprop(k).num)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).bprop = restraints.bprop(k).prop;
    end
end

if isfield(restraints,'pprop')
    for k = 1:length(restraints.pprop)
        res = restraints.pprop(k).num;
        if isempty(res) || isnan(res)
            add_msg_board(sprintf('ERROR: Residue number %i is invalid in PPROP.',restraints.pprop(k).num)); 
            restrain=[];
            cancelled = true;
            return;
        end
        kr = res - res1 + 1;
        restrain(kr).pprop = restraints.pprop(k).prop;
    end
end

if isfield(restraints,'helices')
    for k = 1:length(restraints.helices)
        ha = restraints.helices(k).Nterm;
        he = restraints.helices(k).Cterm;
        if isnan(ha) || isnan(he)
            add_msg_board(sprintf('ERROR: Wrong helix specification %i-%i.',restraints.helices(k).Nterm,restraints.helices(k).Cterm)); 
            restrain=[];
            cancelled = true;
            return;
        end
        for ks = 1:length(restraints.sequence)
            off1 = ks - 1 + res1 - ha;
            off2 = he - (ks - 1 + res1);
            if off1 >=0 && off2 >=0
                restrain(ks).secondary = 1;
            end
            if off1 >=2 && off2 >=2
                restrain(ks).secondary = 3;
            end
        end
    end
end

if isfield(restraints,'strands')
    for k = 1:length(restraints.strands)
        ha = restraints.strands(k).Nterm;
        he = restraints.strands(k).Cterm;
        if isnan(ha) || isnan(he)
            add_msg_board(sprintf('ERROR: Wrong strand specification %I-%i.',restraints.strands(k).Nterm,restraints.strands(k).Cterm)); 
            restrain=[];
            cancelled = true;
            return;
        end
        for ks = 1:length(restraints.sequence)
            off1 = ks - 1 + res1 - ha;
            off2 = he - (ks - 1 + res1);
            if off1 >=0 && off2 >=0
                restrain(ks).secondary = 2;
            end
        end
    end
end

% anchor residue information

Na_indices = [];
if isfield(restraints,'Nanchor') && ~strcmp(restraints.Nanchor,'*')
    [indices,message] = resolve_address(restraints.Nanchor);
    indices(3) = modnum;
    restraints.Nanchor = mk_address(indices);
    Na_indices = indices;
    if message.error
        add_msg_board(sprintf('ERROR: N-terminal anchor residue %s does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    end
    Nname = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
    id = tag2id(Nname,upper(residue_defs.restags));
    N_slc = residue_defs.single_letter_code(id);
    indices_p = indices;
    indices_p(4) = indices_p(4) - 1; % index of previous residue before N terminal anchor
    Npname = model.structures{indices_p(1)}(indices_p(2)).residues{indices_p(3)}.info(indices_p(4)).name;
    id = tag2id(Npname,upper(residue_defs.restags));
    Np_slc = residue_defs.single_letter_code(id);
    restraints.Nanchor_p = mk_address(indices_p);
    anchorNp = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor N atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorNp(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor CA atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorNp(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor C atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorNp(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Nanchor_p),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s before N-terminal anchor O atom does not exist.',restraints.Nanchor_p));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorNp(4,:) = coor;
    end
    anchorN = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Nanchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s N atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorN(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Nanchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s CA atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorN(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Nanchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s C atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorN(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Nanchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For N-terminal anchor residue %s O atom does not exist.',restraints.Nanchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorN(4,:) = coor;
    end
    restraints.anchorN = anchorN;
    restraints.anchorNp = anchorNp;
    restraints.Nseq = [Np_slc N_slc];
else
    restraints.anchorN = [];
    restraints.anchorNp = [];
    restraints.Nseq = '';
end
Ca_indices = [];
if isfield(restraints,'Canchor') && ~strcmp(restraints.Canchor,'*')
    [indices,message] = resolve_address(restraints.Canchor);
    indices(3) = modnum;
    restraints.Canchor = mk_address(indices);
    Ca_indices = indices;
    if message.error
        add_msg_board(sprintf('ERROR: C-terminal anchor residue %s does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    end
    Cname = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
    id = tag2id(Cname,upper(residue_defs.restags));
    C_slc = residue_defs.single_letter_code(id);
    indices_n = indices;
    indices_n(4) = indices_n(4) + 1; % index of next residue after C-terminal anchor
    Cnname = model.structures{indices_n(1)}(indices_n(2)).residues{indices_n(3)}.info(indices_n(4)).name;
    id = tag2id(Cnname,upper(residue_defs.restags));
    Cn_slc = residue_defs.single_letter_code(id);
    restraints.Canchor_n = mk_address(indices_n);
    anchorCn = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Canchor_n),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor N atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorCn(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Canchor_n),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor CA atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorCn(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Canchor_n),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor C atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorCn(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Canchor_n),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For residue %s after C-terminal anchor O atom does not exist.',restraints.Canchor_n));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorCn(4,:) = coor;
    end
    anchorC = zeros(4,3);
    [msg,coor] = get_object(sprintf('%s.N',restraints.Canchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s N atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorC(1,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.CA',restraints.Canchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s CA atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorC(2,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.C',restraints.Canchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s C atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
        return
    else
        anchorC(3,:) = coor;
    end
    [msg,coor] = get_object(sprintf('%s.O',restraints.Canchor),'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: For C-terminal anchor residue %s O atom does not exist.',restraints.Canchor));
        add_msg_board(message.text);
        restrain=[];
        cancelled = true;
       return
    else
        anchorC(4,:) = coor;
    end
    restraints.anchorC = anchorC;
    restraints.anchorCn = anchorCn;
    restraints.Cseq = [C_slc Cn_slc];
else
    restraints.anchorC = [];
    restraints.anchorCn = [];
    restraints.Cseq = '';
end
restraints.Na_indices = Na_indices;
restraints.Ca_indices = Ca_indices;


function labels = get_labels(llist,label)

global model
global hMain

T = 298; 

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

% check whether sites are already labelled and whether all restraint sites
% do exist
lindices=zeros(length(labels),4);
for k=1:length(labels),
    cindices=labels(k).indices;
    if ~isempty(cindices),
        lindices(k,:)=cindices;
    end;
end;
poi=0;
to_do_list{1}=' ';
for k=1:length(llist),
    adr1=llist{k};
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
        add_msg_board(sprintf('This site does not exist in current structure %s',mk_address(1)));
    end;
    found=false;
    for l=1:length(labels),
        diff=ind1-lindices(l,:);
        if sum(abs(diff))==0,
            found=true;
        end;
    end;
    if ~found,
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi} = adr1;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label{k},T);
        hMain.store_undo=false;
        hMain.dynamic_rotamers=false;
        cmd(hMain,command);
    end;
end;

if isfield(model,'sites'),
    labels = label_information(model.sites);
else
    labels = cell(0);
end;

function NO = get_relative_label(libname)

load(libname);
midNO = rot_lib.usefull_atoms.midNO;
pops = rot_lib.calibration.pop;
NO = zeros(1,3);
for k = 1:length(rot_lib.library),
    coor = rot_lib.library(k).ecoor;
    NO = NO + pops(k)*(coor(midNO(1),2:4) + coor(midNO(2),2:4))/2;
end;
NO = NO/(sum(pops));

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

function [rmsd,xyz]=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
xyz = [xmean,ymean,zmean];
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function res = correct_section_address(adr)

if adr(1) == '#'
    resstr = adr(2:end);
else
    resstr = adr;
end
res = str2double(resstr);
if isnan(res)
    res = [];
end
if res - floor(res) > eps
    res = [];
end

function [restrain,number] = mk_beacon_restraint(restrain,NO,res_loop,xyz_beacon,rmean,sigr,res1,number,label1,label2,bindices,resb)

scale_units = 1; % rd_restraints_rigiflex already converts to Å

grace = 0.5; % 5 Å uncertainty of label position

k = res_loop - res1 + 1;
kr = length(restrain(k).r_beacon)+1;
restrain(k).label = NO;
restrain(k).r_beacon(kr).xyz = xyz_beacon;
restrain(k).r_beacon(kr).label1 = label1;
restrain(k).r_beacon(kr).label2 = label2;
restrain(k).r_beacon(kr).bindices = bindices;
restrain(k).r_beacon(kr).resb = resb;
if rmean > 0 && sigr > 0,
    restrain(k).r_beacon(kr).type = 'Gaussian';
    restrain(k).r_beacon(kr).par1 = rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = sqrt(sigr^2 + grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_beacon(kr).type = 'bounds';
    restrain(k).r_beacon(kr).par1 = -rmean*scale_units;
    restrain(k).r_beacon(kr).par2 = -sigr*scale_units;
end;

function [restrain,number] = mk_internal_restraint(restrain,NO1,NO2,resa,resb,rmean,sigr,res1,number,label1,label2)

scale_units = 1; % rd_restraints_rigiflex already converts to Å

grace = 0.5; % 5 Å uncertainty of label position

if resa < resb % restraint must be stored at the later site
    exch = resa;
    resa = resb;
    resb = exch;
    exch = NO1;
    NO1 = NO2;
    NO2 = exch;
end;
k = resa - res1 + 1;
k2 = resb - res1 + 1;
kr = length(restrain(k).r_intern)+1;
restrain(k).label = NO1;
restrain(k2).label = NO2;
restrain(k).r_intern(kr).site = k2;
restrain(k).r_intern(kr).label1 = label1;
restrain(k).r_intern(kr).label2 = label2;
restrain(k).r_intern(kr).resb = resb;
if rmean > 0 && sigr > 0,
    restrain(k).r_intern(kr).type = 'Gaussian';
    restrain(k).r_intern(kr).par1 = rmean*scale_units;
    restrain(k).r_intern(kr).par2 = sqrt(sigr^2 + 2*grace^2)*scale_units;
    number = number + 1;
else
    restrain(k).r_intern(kr).type = 'bounds';
    restrain(k).r_intern(kr).par1 = -rmean*scale_units;
    restrain(k).r_intern(kr).par2 = -sigr*scale_units;
end;

function [restrain,number] = mk_depth_restraint(restrain,NO,res,zmean,sigz,res1,number,label)

scale_units = 1; % rd_restraints_rigiflex already converts to Å

k = res - res1 + 1;
kr = length(restrain(k).depth)+1;
restrain(k).depth(kr).label_type = label;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).depth(kr).site = 'label';
else
    restrain(k).depth(kr).site = 'CA';
end;
if zmean > 0 && sigz > 0,
    restrain(k).depth(kr).type = 'Gaussian';
    restrain(k).depth(kr).par1 = zmean*scale_units;
    restrain(k).depth(kr).par2 = sigz*scale_units;
    number = number + 1;
else
    restrain(k).depth(kr).type = 'bounds';
    restrain(k).depth(kr).par1 = -zmean*scale_units;
    restrain(k).depth(kr).par2 = -sigz*scale_units;
end;

function [restrain,number] = mk_oligomer_restraint(restrain,NO,res,n,rmean,sigr,res1,number,label)

scale_units= 1; % rd_restraints_rigiflex already converts to Å

k = res - res1 + 1;
kr = length(restrain(k).oligomer)+1;
restrain(k).oligomer(kr).label_type = label;
if ~isempty(NO),
    restrain(k).label = NO;
    restrain(k).oligomer(kr).site = 'label';
else
    restrain(k).oligomer(kr).site = 'CA';
end;
if rmean > 0 && sigr > 0,
    restrain(k).oligomer(kr).type = 'Gaussian';
    restrain(k).oligomer(kr).par1 = rmean*scale_units;
    restrain(k).oligomer(kr).par2 = sigr*scale_units;
    number = number + 1;
else
    restrain(k).oligomer(kr).type = 'bounds';
    restrain(k).oligomer(kr).par1 = -rmean*scale_units;
    restrain(k).oligomer(kr).par2 = -sigr*scale_units;
end;
restrain(k).oligomer(kr).n = n;

function handles = set_progress_interface(handles)

handles.text_percent_time.String = '0';
handles.text_max_trials.String = '0';
handles.text_success.String = '0';
handles.text_time_left.String = 'Idle.';
handles.text_dmg_fail.String = '0';
handles.text_auxiliary_fail.String = '0';
handles.text_core_fail.String = '0';
handles.text_linker_fail.String = '0';
handles.text_clash_fail.String = '0';
handles.text_xlink_fail.String = '0';
handles.text_SANS_fail.String = '0';
handles.text_SAXS_fail.String = '0';

switch handles.progress
    case 0
        handles.text_ri2.String = '% of max. trials spent';
        handles.text_ri3.String = 'rigid-body models found';
        handles.text_ri4.String = '% failed matrix geometry';
        handles.text_ri5.String = '% failed auxiliary restraints';
        handles.text_ri6.String = '% failed core restraints';
        handles.text_ri7.String = '% failed linker restraints';
        handles.text_ri8.String = '% rigid-body clashes';
        handles.text_ri9.String = '% failed cross-link restraints';
        handles.text_ri10.String = '% failed SANS fitting';
        handles.text_ri11.String = '% failed SAXS fitting';
    case 1
        handles.text_ri2.String = '% of max. trials spent';
        handles.text_ri3.String = '(rigid body-arrangement)';
        handles.text_ri4.String = '(RNA section)';
        handles.text_ri5.String = '% section trials';
        handles.text_ri6.String = '% section time';
        handles.text_ri7.String = ' ';
        handles.text_ri8.String = ' ';
        handles.text_ri9.String = ' ';
        handles.text_ri10.String = ' ';
        handles.text_ri11.String = ' ';
    case 2
        handles.text_ri2.String = '% of max. trials spent';
        handles.text_ri3.String = 'loop models found';
        handles.text_ri4.String = 'backbone models found';
        handles.text_ri5.String = '% restraint violations';
        handles.text_ri6.String = '% internal clashes';
        handles.text_ri7.String = '% clashes with protein';
        handles.text_ri8.String = '% loop side chain clashes';
        handles.text_ri9.String = '% protein side chain clashes';
        handles.text_ri10.String = '% loop closure failures';
        handles.text_ri11.String = '% Ramachandran mismatch';
    case 2
        handles.text_ri2.String = 'combinations will be tested';
        handles.text_ri3.String = 'models assembled';
        handles.text_ri4.String = 'backbone models found';
        handles.text_ri5.String = '% restraint violations';
        handles.text_ri6.String = '% internal clashes';
        handles.text_ri7.String = '% clashes with protein';
        handles.text_ri8.String = '% loop side chain clashes';
        handles.text_ri9.String = '% protein side chain clashes';
        handles.text_ri10.String = '% loop closure failures';
        handles.text_ri11.String = '% Ramachandran mismatch';
end


% --- Executes on button press in checkbox_exhaustive.
function checkbox_exhaustive_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_exhaustive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_exhaustive

if isfield(handles,'restraints')
    [handles,restraints] = analyze_exhaustive(handles,handles.restraints);
    handles.restraints = restraints;
else
    handles.text_exhaustive_resolution.String = sprintf('n.a.');
    handles.edit_max_trials.String = sprintf('%i',handles.max_trials);
end
guidata(hObject,handles);

function [handles,restraints] = analyze_exhaustive(handles,restraints)

tph = 1000000; % trials per hour, to be replaced with machine-specific value
target_resolution = 3; % lowest resolution to be considered sensible with spin labels
min_approach = 5; % minimal approach of two reference points [Å]
max_extension = 180; % 180 maximum distance between any two reference points [Å]

[m,~] = size(restraints.lb);
% augment lower and upper bounds
for k1 = 1:m-1
    for k2 = k1+1:m
        if restraints.lb(k1,k2) < min_approach
            restraints.lb(k1,k2) = min_approach;
            restraints.lb(k2,k1) = min_approach;
        end
        if restraints.ub(k1,k2) < 0.1 % the unset bounds are zero
            restraints.ub(k1,k2) = max_extension;
            restraints.ub(k2,k1) = max_extension;
        end
    end
end
ntrials = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
res = target_resolution;
if restraints.maxtrials > 0
    while ntrials > restraints.maxtrials
        target_resolution = target_resolution + 0.1;
        [ntrials,res] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
    end
end
trials = tph*handles.max_time;
if trials > handles.max_trials
    trials = handles.max_trials;
end
while ntrials > trials
    target_resolution = target_resolution + 0.1;
    [ntrials,res] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
end
handles.text_exhaustive_resolution.String = sprintf('%4.1f',res);
if handles.checkbox_exhaustive.Value
    handles.edit_max_trials.String = sprintf('%i',ntrials);
    restraints.maxtrials = ntrials;
else
    handles.edit_max_trials.String = sprintf('%i',handles.max_trials);
end

function handles = mk_RNA_report_distributions(fid_report,handles)

sadr = mk_address(handles.diagnostics.snum);

% fid_report = fopen(report_name,'wt');

fprintf(fid_report,'\n\n--- DEER restraint report ---\n');

fprintf(fid_report,'\nRestraint file used: %s\n',handles.restraint_file);

restraints = handles.restraints;

distributions = handles.distributions;
poi = length(distributions);

if isfield(restraints.RNA,'DEER')
    for kr = 1:length(restraints.RNA.DEER)
        ind1 = resolve_address(sprintf('%s%s',sadr,restraints.RNA.DEER(kr).adr1));
        ind2 = resolve_address(sprintf('%s%s',sadr,restraints.RNA.DEER(kr).adr2));
        if ~isempty(ind1) && ~isempty(ind2)
            poi = poi + 1;
            distributions(poi).type = 'RNA link restraint';
            distributions(poi).label1 = restraints.RNA.DEER(kr).label1;
            distributions(poi).label2 = restraints.RNA.DEER(kr).label2;
            distributions(poi).ind1 = ind1;
            distributions(poi).ind2 = ind2;
            ind1(1) = handles.diagnostics.snum;
            ind2(1) = handles.diagnostics.snum;
            distr = zeros(handles.diagnostics.success,301); % assume that get_distribution returns a distribution with 300 points
            sum_distr = zeros(1,301);
            fprob = 0;
            sum_var = 0;
            for km = 1:handles.diagnostics.success
                if handles.RNA_link_success(km)
                    ind1(3) = km;
                    ind2(3) = km;
                    NO_pos1 = get_NO_pos(ind1,distributions(poi).label1,298);
                    NO_pos2 = get_NO_pos(ind2,distributions(poi).label2,298);
                    [rax,sim_distr] = get_distribution(NO_pos1,NO_pos2,0.05);
                    rmean = sum(rax.*sim_distr);
                    cvar = sum(sim_distr.*(rax-rmean).^2);
                    distr(km,:) = sim_distr;
                    sum_distr= sum_distr + sim_distr*handles.diagnostics.probabilities(km);
                    sum_var = sum_var + cvar*handles.diagnostics.probabilities(km);
                    fprob = fprob + handles.diagnostics.probabilities(km);
                end
            end
            sum_distr = sum_distr/fprob;
            mean_var = sum_var/fprob;
            distributions(poi).all_distr = distr;
            distributions(poi).sum_distr = sum_distr;
            distributions(poi).r_req = restraints.RNA.DEER(kr).r;
            distributions(poi).sigr_req = restraints.RNA.DEER(kr).sigr;
            [~,ctag1,~,resnum1] = mk_address_parts(ind1);
            [~,ctag2,~,resnum2] = mk_address_parts(ind2);
            distributions(poi).adr1 = sprintf('(%s)%i',ctag1,resnum1);
            distributions(poi).adr2 = sprintf('(%s)%i',ctag2,resnum2);
            distributions(poi).rax = rax;
            fprintf(fid_report,'\nRNA link restraint (%s)%i to (%s)%i\n',ctag1,resnum1,ctag2,resnum2);
            meanr = 10*sum(sum_distr.*rax);
            sigr = sqrt(sum(sum_distr.*(10*rax-meanr).^2));
            distributions(poi).meanr = meanr;
            distributions(poi).sigr = sigr;
            fprintf(fid_report,'<r> = %4.1f Å, sigr = %4.1f Å for restraint <r> = %4.1f Å, sigr = %4.1f Å)',meanr,sigr,restraints.RNA.DEER(kr).r,restraints.RNA.DEER(kr).sigr);
            fprintf(fid_report,'Rigidity index for ensemble: %5.2f; to restraint: %5.2f\n',sqrt(mean_var/sigr^2),sqrt(mean_var/restraints.RNA.DEER(kr).sigr^2));
            sim_distr = (rax-restraints.DEER(kr).r/10).^2/(2*0.01*restraints.DEER(kr).sigr^2);
            sim_distr = exp(-sim_distr);
            sim_distr = sim_distr/sum(sim_distr);
            overlap = sum(min([sum_distr;sim_distr]));
            fprintf(fid_report,' (overlap: %6.4f).\n',overlap);
            mean_overlap = 0;
            overlap_counter = 0;
            for km1 = 1:handles.diagnostics.success-1
                for km2 = km1+1:handles.diagnostics.success
                    if handles.RNA_link_success(km1) && handles.RNA_link_success(km2)
                        overlap_counter = overlap_counter + 1;
                        mean_overlap = mean_overlap + sum(min([distr(km1,:);distr(km2,:)]));
                    end
                end
            end
            mean_overlap = mean_overlap/overlap_counter;
            distributions(poi).overlap = mean_overlap;
            fprintf(fid_report,'Mean overlap of distributions for the ensemble: %6.4f.\n',mean_overlap);
        end
    end
    
    handles.distributions = distributions;
end

function curve = load_SAXS_curve(fname)

fid = fopen(fname);
if fid==-1
    curve = [];
    add_msg_board('Warning. Loading of SAXS curve failed');
    return;
end
nl=0;
curve = zeros(10000,4);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    %         fprintf(1,'%s\n',tline); % echo for debugging
    if nl > 0 % skip first line
        dataset = str2num(tline);
        ncol = length(dataset);
        curve(nl,1:ncol) = dataset;
    end
    nl = nl + 1;
end
curve = curve(1:nl-1,:);
%curve(:,1) = curve(:,1);
fclose(fid);


function info = set_info(chain,snum)

info.idCode = sprintf('X%i',100+snum);
info.class = 'MMM model';
info.depDate = date;
info.title = 'RigiFlex model';
info.chain_tags = ':A:';
info.chain_ids = 1;
info.center = [0,0,0];
info.B_range = [0,20];
info.authors = 'MMM';
info.molecule = 'chimera';
info.organism = 'Dragon musculus';
info.remarks = [4 5 6];
info.SSbonds = [];
info.Modeller_obj = [];
info.Modeller_sid = [];
info.rotamers = false;
info.missing = {};
info.site_tags = ':';
info.sites = [];
info.keywords = '';
info.metal = [];
info.resolution = [];
info.references = [];
info.insertions = [];
info.alternate = false;
info.atoms = chain.atoms{1};
info.residues = length(chain.residues{1});
