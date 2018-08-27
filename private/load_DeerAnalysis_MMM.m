function handles=load_DeerAnalysis_MMM(handles,stand_alone)
% Loads output files of DeerAnalysis (original data, form factor, distance
% distribution)
%
% G. Jeschke, 2009

global general

if ~exist('stand_alone','var'),
    stand_alone = false;
end;

handles.updated=0;
if ~stand_alone,
    handles=set_defaults_MMM(handles); % set all status variables to defaults
end;

my_path=pwd;
cd(general.DEER_files);

ext='*.txt;';
[fname,pname]=uigetfile(ext,'Load DeerAnalysis result');
if isequal(fname,0)||isequal(pname,0),
    return; 
end;
reset_user_paths(pname);
general.DEER_files=pname;
cd(pname);

% separate filename from extension
dots=findstr('.',fname);
ext_pos=length(fname);
if ~isempty(dots), ext_pos=dots(length(dots))-1; end;
bas_name=fname(1:ext_pos-4); % name without extension and result file extension
figname=['Deer simulation - ' fname]; % tell user, which file is current
if ~stand_alone,
    set(handles.DEER,'Name',figname);
end;

handles.project_dir=pname;
handles.bas_name=bas_name;
handles.source_file=[pname fname];


orig_name=[bas_name '_bckg.dat'];
data=load(orig_name);
handles.texp=data(:,1);
handles.vexp=data(:,2);
handles.v_orig=data(:,2);
handles.t_orig=1000*(data(:,1)-data(1,1));

ff_name=[bas_name '_fit.dat'];
data=load(ff_name);
handles.tdip=data(:,1);
handles.cluster=data(:,2);
handles.exp_depth=1-data(end,3);
if ~stand_alone,
    set(handles.checkbox_form_factor,'Enable','on');
end;

distr_name=[bas_name '_distr.dat'];
data=load(distr_name);
handles.rexp=data(:,1);
handles.dexp=data(:,2);

cd(my_path);
