function job_check
% function job_check
% checks whether the temporary directory contains information on background
% jobs, for instance Modeller jobs and reinstalls periodic test for job
% completion and result import call
%
% G. Jeschke, 2011

global general
tmp_dir=dir(general.tmp_files); % load temporary directory
files=numel(tmp_dir);
for k=1:files, % check for and mark old files
    if ~tmp_dir(k).isdir,
        if strcmpi(tmp_dir(k).name(1:4),'job-'), % this is a job file
            add_msg_board(sprintf('Starting periodic check for background job %s',tmp_dir(k).name));
            load(tmp_dir(k).name);
            switch job.type
                case 'Modeller'
                    t = timer('TimerFcn',@check_for_Modeller_completion, 'StartDelay', 2.0, 'Period', 120.0, 'ExecutionMode', 'fixedSpacing' ,'UserData', job);
                    start(t);
                otherwise
                    add_msg_board('Job file indicates unknown background job type.');
            end;
        end;
    end;
end;

