function tmp_check
% function tmp_check
% checks whether the temporary directory should be cleaned up, if so,
% whether there are old files and deletes any old files (after asking the
% user)
%
% G. Jeschke, 2009

global general

time_stamp=strcat(general.tmp_files,'time_stamp.mat');
current=now;
if exist(time_stamp,'file'), % check whether directory was cleaned up and if so, load the corresponding time stamp
    load(time_stamp);
    add_msg_board(sprintf('Temporary files were last cleaned %s',datestr(cleaned)));
    if current-cleaned>general.tmp_check,
        cleaned=current;
        save(time_stamp,'cleaned');
    else
        return; % do nothing if the last check was recent
    end;
end;

% disp(general.tmp_files);
tmp_dir=dir(general.tmp_files); % load temporary directory
files=numel(tmp_dir);
to_delete=zeros(1,files);
for k=1:files, % check for and mark old files
    if current-tmp_dir(k).datenum>general.old_file && ~tmp_dir(k).isdir,
        to_delete(k)=1;
        disp(sprintf('File %s is %6.2f days old.',tmp_dir(k).name,current-tmp_dir(k).datenum));
    end;
end;

if sum(to_delete)>0, % check whether there are any old files
    add_msg_board(sprintf('Found %i old temporary files %i',sum(to_delete)));
    button = questdlg(sprintf('Delete the temporary files older than %i days now?',general.old_file),sprintf('Found %i old temporary files %i',sum(to_delete)),'Yes','No','Yes');
    if strcmp(button,'Yes'),
        oldstate=recycle;
        recycle('off');
        for k=1:files,
            if to_delete(k),
                fname=strcat(general.tmp_files,tmp_dir(k).name);
                delete(fname);
            end;
        end;
        recycle(oldstate);
        add_msg_board('Old temporary files deleted.');
    else
        add_msg_board('Deleting old temporary files refused by user.');
    end;
else
    add_msg_board('There were no old temporary files to be deleted.');
end;

cleaned=current;
save(time_stamp,'cleaned');
