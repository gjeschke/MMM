% function conv_AVcloud2dx(in_filename, out_filename)
% to turn fluorescence AV .xyz output file into .dx volmap data as for MMM & VMD
% 
% DX file header
% # Data calculated by the VMD volmap function
% # object 1 class gridpositions counts 155 123 153
% # origin 4.55 -14.95 -19.55
% # delta 0.1 0 0
% # delta 0 0.1 0
% # delta 0 0 0.1
% # object 2 class gridconnections counts 155 123 153
% # object 3 class array type double rank 0 items 2916945 data follows
% data: z fastest, then y, then x
%
% Example:
% in_filename = 'Q388A488(D).xyz';
% out_filename = 'Q388A488(D)_res1.dx';
% 
% Daniel Klose, 2020
%
function success = conv_AVcloud2dx(in_filename, out_filename, volume_res, plot_flag)


volume_box_size = 81; % Size of Voxel Cube in Angstrom

cloud_probability1 = 0.995; % 0.995 is MMM default (currently unused)
cloud_probability2 = 0.5; % 0.995 is MMM default (currently unused)

% currently no interpolation: ntimes = 1; % recursive interpolation degree, 1 or 2 seems best
num_upscale = 1e5; % may be 1e5, normalization upscale to evade 1e-6 limit in VMD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
success = 0;

if nargin < 4
    plot_flag = 0;
end
if nargin < 3
    plot_flag = 0;
    volume_res = 3; % grid resolution in Ang.
end
if nargin < 2
    plot_flag = 0;
    volume_res = 3; % grid resolution in Ang.
    out_filename = [in_filename(1:end-4) '_res' num2str(volume_res) '.dx'];
end
if nargin == 0
    add_msg_board('conv_AVcloud2dx: Please provide input filename!');
    return;
end

fid = fopen(in_filename);
if fid == -1
    add_msg_board('conv_AVcloud2dx: Could not open input file!');
    return;
end
tline = fgetl(fid);
Npoints = sscanf(tline,'%f');
description = fgetl(fid);
dye_type = description(1);

% read AV data
tline = fgetl(fid);
data = [];
while strcmp(tline(1:2),[dye_type ' ']) % indicates last line
    data = vertcat(data, sscanf(tline(3:end),'%f %f %f')');
    tline = fgetl(fid);
end
middle = sscanf(tline(5:end),'%f %f %f');
fclose(fid);

% calculate data points on grid
box_center = middle;
box_origin = middle - volume_box_size/2;
% currently defined manually:
% volume_box_size = ceil(max([abs(max(data(:,1))-min(data(:,1))), abs(max(data(:,2))-min(data(:,2))), abs(max(data(:,3))-min(data(:,3)))]) + 2*volume_res);
box = zeros(volume_box_size/volume_res+1,volume_box_size/volume_res+1,volume_box_size/volume_res+1);
center_int = round((volume_box_size/volume_res+1)/2);
for i = 1:length(data(:,1)) % loop over all frames (numbering acc. to array!!)
            coord_diff = data(i,:)'-box_center;
            coord_index = round(coord_diff./volume_res) + center_int;
            % build up density (though values might be unique
            box(coord_index(1),coord_index(2),coord_index(3)) = box(coord_index(1),coord_index(2),coord_index(3)) + 1;
end
total_volume = (volume_res^3)*length(find(box~=0)); % in Ang^3

% normalize to probability
cartesian_grid = permute(box,[3 2 1]);
cartesian_grid = cartesian_grid./sum(sum(sum(cartesian_grid)));

dx_data = zeros([3 ceil(numel(box)/3)]);
dx_data(:,end) = NaN;
dx_data(1:length(cartesian_grid(:))) = cartesian_grid(:).*num_upscale;
dim_change_offset_new = sum(isnan(dx_data(:)));
dx_data = dx_data';

try
    fid = fopen(out_filename,'w');
    fprintf(fid,'%s\n',['# file calculated by AVcloud2dx (DAKL, 2019) from ' description]);
    fprintf(fid,'%s %i %i %i\n','object 1 class gridpositions counts', length(box(1,1,:)), length(box(1,:,1)),length(box(:,1,1)));
    fprintf(fid,'%s %.3f %.3f %.3f\n','origin', box_origin);
    fprintf(fid,'%s %.1f %i %i\n','delta', volume_res, 0, 0);
    fprintf(fid,'%s %i %.1f %i\n','delta', 0, volume_res, 0);
    fprintf(fid,'%s %i %i %.1f\n','delta', 0, 0, volume_res);
    fprintf(fid,'%s %i %i %i\n','object 2 class gridconnections counts', length(box(1,1,:)), length(box(1,:,1)),length(box(:,1,1)));
    fprintf(fid,'%s %i %s\n','object 3 class array type double rank 0 items', length(box(:)), 'data follows');
    fclose(fid);

    dlmwrite(out_filename,dx_data,'-append','delimiter',' ');
    fid = fopen(out_filename,'a'); % append writing
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','object "AV density via Matlab [A^-3]" class field');
catch
    add_msg_board('conv_AVcloud2dx: Could not open/write to output file!');
    return;
end
% calc isoval1
cartesian_grid_sorted = sort(cartesian_grid(:),'descend');
cutoff_index = find(cumsum(cartesian_grid_sorted)>=cloud_probability1,1,'first');
isoval1 = cartesian_grid_sorted(cutoff_index);

% calc percentage1
cartesian_grid_sorted = sort(cartesian_grid(:));
cartesian_grid_logic = cartesian_grid_sorted >= isoval1;
probability_here1 = sum(cartesian_grid_sorted.*cartesian_grid_logic);

% calc isoval2
cartesian_grid_sorted = sort(cartesian_grid(:),'descend');
cutoff_index = find(cumsum(cartesian_grid_sorted)>=cloud_probability2,1,'first');
isoval2 = cartesian_grid_sorted(cutoff_index);

% calc percentage2
cartesian_grid_sorted = sort(cartesian_grid(:));
cartesian_grid_logic = cartesian_grid_sorted >= isoval2;
probability_here2 = sum(cartesian_grid_sorted.*cartesian_grid_logic);

if plot_flag
    figure();
    isosurface(permute(box,[2 3 1]),1);
    title('org. equal probabilities, isovalue = 1');
    figure();
    isosurface(permute(cartesian_grid,[2 1 3]),isoval1);
    title(['isovalue1 = ' num2str(isoval1)]);
    figure();
    isosurface(permute(cartesian_grid,[2 1 3]),isoval2);
    title(['isovalue2 = ' num2str(isoval2)]);
end

try
    fprintf(fid,'%s\n',['# Isovalue1 of ' num2str(isoval1*num_upscale) ' => probability of ' num2str(probability_here1)]);
    fprintf(fid,'%s\n',['# Isovalue2 of ' num2str(isoval2*num_upscale) ' => probability of ' num2str(probability_here2)]);
    fprintf(fid,'%s\n',['# delete trailing NaNs: ' num2str(dim_change_offset_new)]);
    fprintf(fid,'%s\n',['# Total volume: ' num2str(total_volume)]);
    
    fclose(fid);
catch
    add_msg_board('conv_AVcloud2dx: Could not open/write to output file!');
    return;
end
success = 1;

end





