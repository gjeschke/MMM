function [file_list,pop] = rd_ensemble_description(fname)

file_list = cell(1,1000);
pop = zeros(1,1000);
poi = 0;

fid=fopen(fname);
if fid == -1
    add_msg_board('Warning: File list does not exist');
    file_list = {};
    return;
end

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        if ~contains(tline,'%')
            myline = textscan(tline,'%s');
            args = myline{1};
            arg1 = char(args(1));
            if ~isempty(arg1)
                poi = poi + 1;
                file_list{poi} = arg1;
                pop(poi) = str2double(char(args(2)));
            end
        end
    end
end

fclose(fid);

file_list = file_list(1:poi);
pop = pop(1:poi);
