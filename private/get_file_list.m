function [my_list,comments] = get_file_list(fname)
% Reads a list of file names from a list file, % is the comment character,
% comments are skipped
% the first argument on a line is considered as a file name, other
% arguments are skipped
% output is stored in a cell vector
% comments are returned in a cell string


my_list = cell(1,10000);
comments = cell(1,100);
poi = 0;

fid=fopen(fname);
if fid == -1
    add_msg_board('Warning: File list does not exist');
    my_list = {};
    return;
end

compoi = 0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            compoi = compoi + 1;
            comments{compoi} = tline(k(1)+1:end);
        else
            myline = textscan(tline,'%s');
            args=myline{1};
            arg1 = char(args(1));
            if ~isempty(arg1)
                poi = poi + 1;
                my_list{poi} = arg1;
            end
        end
    end
end

fclose(fid);

my_list = my_list(1:poi);
comments = comments(1:compoi);
