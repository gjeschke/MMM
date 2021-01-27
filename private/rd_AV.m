% function [rax,distr] = rd_AV(fname);
%
% read in two AV .xyz files in fname
% and calculates distance distribution [rax, distr]
% fname is 1x2 cell from fullfile()
function [rax,distr] = rd_AV(fname,volume_res)

if nargin < 2
    volume_res = 3;
end

if numel(fname) ~= 2
    add_msg_board('Warning: Two AV files are required for a distance distribution.');
    return
end

dx_filename = cell(1,2);
filename = cell(1,2);
for i = 1:2
    dx_filename{i} = [fname{i}(1:end-4) '_res' num2str(volume_res) '.dx']; % standard dx_filename (fullfile() format)
    [path,filename{i},~] = fileparts(dx_filename{i});
    if isfile(dx_filename{i})
        add_msg_board('Found dx grid data file');
    else
        add_msg_board('Conversion to dx grid data...');
        success = conv_AVcloud2dx(fname{i},dx_filename{i},volume_res);
        if ~success
            add_msg_board('Warning: Conversion of AV file to .dx grid failed.');
            return
        end
        add_msg_board('... done.');
    end
end

distr_filename = fullfile(path,[filename{1} '_' filename{2} '_distr.dat']);
if isfile(distr_filename)
    add_msg_board('Found AV distance distribution file');
    data = importdata(distr_filename);
    rax = data(:,1);
    distr = data(:,2);
else
    add_msg_board('Calculating distance distribution...');
    [rax, distr] = conv_dx2distr(dx_filename{1}, dx_filename{2});
    dlmwrite(distr_filename,[rax distr],'delimiter','\t','precision',16);
    add_msg_board('... done.');
end
    


