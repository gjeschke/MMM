function check_for_Modeller_completion(obj, event, string_arg)

global general
global hMain

job=get(obj,'UserData');
add_msg_board(sprintf('Modeller job with log file %s',job.log));

key=job.finished;

fname=fullfile(general.tmp_files, job.log);

is_completed=false;

fid=fopen(fname,'r');

if fid~=-1,
    isheader=1;
    while isheader,
        tline = fgetl(fid);
        if ~ischar(tline), break, end;
        if strfind(tline,key),
            isheader=true;
            is_completed=true;
        end;
    end;
end;

fclose(fid);

if is_completed,
    runtime=etime(clock,job.started);
    hours=floor(runtime/3600);
    minutes=floor((runtime-3600*hours)/60);
    seconds=round(runtime-3600*hours-60*minutes);
    add_msg_board(sprintf('is completed after %i h %i min % i s.',hours,minutes,seconds));
    stop(obj);
    button = questdlg('Do you want to import results? (otherwise, results can only be imported after closing and reopening MMM)','Modeller job completed.','Yes','No','Yes');
    if strcmp(button,'Yes'),
        set(hMain.figure,'Pointer','watch');
        drawnow
        import(job);
        add_msg_board('Imported structure is not automatically displayed.');
        add_msg_board('Use command "show [:] coil" or any other means of display.');
        set(hMain.figure,'Pointer','arrow');
        delete(fullfile(general.tmp_files,[job.name '.mat']));
    end;
else
    add_msg_board('is still running...');
end;

drawnow


function import(job)

global general
global hMain
global graph_settings
global model

pname=general.tmp_files;
idCode=job.targ_ID;
logname=job.log;
ensemble=assess_ensemble(logname,'',job.ensemble_size,job.min_GA341);
if isempty(ensemble),
    msg='ERROR: No ensemble model could be generated.';
    add_msg_board(msg);
    return
end;
ensemble_size=length(ensemble);
if ensemble_size < job.ensemble_size,
    add_msg_board(sprintf('Warning: only %i models have sufficient GA341 score, while %i models were requested.',ensemble_size,job.ensemble_size));
    add_msg_board('Importing smaller ensemble. Consider increasing the number of models.');
end;
[message,snum,stag,models]=rd_pdb_ensemble(ensemble,idCode,pname);
if message.error,
    add_msg_board(message.text);
else
    correct_Modeller_ensemble(snum,job.ares,job.chain_tags);
    add_msg_board(sprintf('Ensemble read into structure %i with tag [%s].',snum,stag));
    add_msg_board(sprintf('This ensemble has %i chain models',models));
end;

id=[stag 'r'];
exists=tag2id(id,model.structure_tags);
while ~isempty(exists),
    id=[id 'a'];
    exists=tag2id(id,model.structure_tags);
end;
add_msg_board('Now removing spin labels...');
[repnum,msg]=replace(snum,':R1A:IA1:',false,false);
add_msg_board(sprintf('%i spin labels were removed.',repnum));
add_msg_board('Now repacking side chains...');
for modnum=1:models,
    message=repack(snum,modnum,0,job.clean_targ_seq);
end;

handles=guidata(hMain.figure);

if hMain.virgin,
    hMain.virgin=0;
    % initialize display
    axes(handles.axes_model);
    cla(handles.axes_model,'reset');
    axis(handles.axes_model,'equal');
    axis(handles.axes_model,'off');
    set(handles.axes_model,'Clipping','off');
    set(hMain.figure,'Renderer','opengl');
    hold(handles.axes_model,'on');
    guidata(handles.axes_model,handles);
end;


% initialize display
axes(handles.axes_model);
axis(handles.axes_model,'equal');
axis(handles.axes_model,'off');
set(handles.axes_model,'Clipping','off');
set(hMain.figure,'Renderer','opengl');
hold(handles.axes_model,'on');

view(handles.axes_model,graph_settings.az,graph_settings.el);
lighting(handles.axes_model,'gouraud'); 
material(handles.axes_model,'shiny');

axes(handles.axes_frame);
cla(handles.axes_frame,'reset');
axis(handles.axes_frame,'equal');
axis(handles.axes_frame,'off');
hold(handles.axes_frame,'on');

plot3(handles.axes_frame,[0 1],[0 0],[0 0],'r','LineWidth',2); % x axis
plot3(handles.axes_frame,[0 0],[0 1],[0 0],'g','LineWidth',2); % y axis
plot3(handles.axes_frame,[0 0],[0 0],[0 1],'b','LineWidth',2); % z axis
axis(handles.axes_frame,[-0.1,1.1,-0.1,1.1,-0.1,1.1]);

guidata(handles.axes_frame,handles);

set_view([],false);

guidata(handles.axes_model,handles);
