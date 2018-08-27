function handback=get_dataset_MMM(handles)
%
% Load data set in the selected format and save it in data structure handles
% simplified version for MMM
%
% G. Jeschke, 2009

global general

handles.updated=0;
handles=set_defaults(handles); % set all status variables to defaults

my_path=pwd;
cd(general.exp_files);

ext='*.DTA;';
[fname,pname]=uigetfile(ext,'Load experimental dataset');
if isequal(fname,0)||isequal(pname,0),
    handback=handles;
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
figname=['Deer simulation - ' fname]; % tell user, which file is current
set(handles.DEER,'Name',figname);
handles.project_dir=pname;
handles.bas_name=bas_name;
handles.source_file=[pname fname];

[x,~,z,vb]=get_elexsys_MMM(bas_name);


[m,~]=size(z); % catches cases where 2D data without explicit summation over nuclear modulation are offered
if m > 1
    z = sum(z);
end

% x=x-x(1)*ones(size(x));

handles.t_orig=x;
handles.v_orig=z;

[texp,vexp,zt,dt]=pre_process_MMM(x,z);

handles.dt=dt;
handles.zero_time=zt;
set(handles.edit_zero_time,'String',sprintf('%5i',handles.zero_time));

handles.bas_name=bas_name;

handles.texp=texp'/1000;
handles.vexp=vexp';
handles.vb=vb;

handles.rexp=1.5:0.05:10;
handles.dexp=[];

handback=handles;

cd(my_path);

