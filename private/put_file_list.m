function fname = put_file_list(fname,my_list,selection,append,aux_data,comments)
% Writes a list of file names to a list file
%
% fname     file name of the list file
% my_list   cell vector of file names comprising the list
% selection optional vector of selected entries in my_list, if present,
%           only the selected entries are written out, if not present, the
%           whole list is written out, if empty, the output file is empty,
%           use NaN if you want to skip this argument (whole list will be
%           written then)
% append    optional flag, if true and file fname already exists, the
%           entries are appended to the existing file, defaults to false
% aux_data  optional additional numerical data, they are written in the
%           form %12.6f if in size [m,n] of the array, either matches the 
%           length of the selection or n matches it and m = 1,
% comments  optional comment lines, cell string, added at the beginning of
%           the file
%
% fname     the input argument is echoed on success, as output argument,
%           fname is empty when an error occured
%
% G. Jeschke, 2.7.2019

if ~exist('selection','var') || sum(isnan(selection))
    selection = 1:length(my_list);
end

if ~exist('append','var')
    append = false;
end

if append 
    fid=fopen(fname,'at');
else
    fid=fopen(fname,'wt');
end

if ~exist('aux_data','var')
    aux_data = [];
else
    [m,n] = size(aux_data);
    if n == length(selection) && m == 1
        aux_data = aux_data.';
        n = 1;
    elseif m ~= length(selection)
        aux_data = [];
        add_msg_board(sprintf('Warning [put_file_list]: Auxiliary data matrix [%i,%i] does not fit length of file list selection (%i). Auxiliary data is skipped.',m,n,length(selection)));
    end
end

if ~exist('comments','var')
    comments = {};
end

if fid == -1
    add_msg_board('ERROR: File list cannot be written');
    fname = '';
    return;
end

for k = 1:length(comments)
    fprintf(fid,'%% %s\n',comments{k});
end

for k = 1:length(selection)
    fprintf(fid,'%s',my_list{selection(k)});
    if ~isempty(aux_data)
        for k2 = 1:n
            fprintf(fid,'%12.6f',aux_data(k,k2));
        end
    end
    fprintf(fid,'\n');
end

fclose(fid);

