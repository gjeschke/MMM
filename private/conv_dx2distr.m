% function [rax, distr] = conv_dx2distr(in_file1, in_file2)
% calc distance distribution between 2 .dx volmap datasets
% use .dx files from conv_AVcloud2dx()
% or from VMD
% 
% Output: Distance distribution [rax, distr] (in nm)
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
% in_file1 = 'Q388A488(D)_res1.dx';
% in_file2 = 'S475A594(A)_res1.dx';
%
% Daniel Klose, 2020
%
function [rax, distr] = conv_dx2distr(in_file1, in_file2)

% load dataset 1
data_struct1 = importdata(in_file1,' ',8);
data1 = data_struct1.data';
% description1 = data_struct1.textdata{1};
grid_counts1 = (sscanf(data_struct1.textdata{2},'%*s %*i %*s %*s %*s %i %i %i'))'; % x, y, z convention
origin1 = sscanf(data_struct1.textdata{3},'%*s %f %f %f'); % get origin vector
delta1 = sscanf(data_struct1.textdata{4},'%*s %f %*f %*f'); % read grid spacing
% delta_after_interp = delta_old/(2^ntimes);

% dim_after_interp = arrayfun(@(x) 2^ntimes*x-(2^ntimes-1) ,grid_counts);
cartesian_grid1 = zeros([grid_counts1(3) grid_counts1(2) grid_counts1(1)]);

if size(data1(:,1)) ~= 3
    error('Wrong dimension of input data!');
end
dim_change_offset1 = sum(isnan(data1(:)));

cartesian_grid1(:) = data1(1:(end-dim_change_offset1)); % data must be 3 x n
cartesian_grid1 = permute(cartesian_grid1,[3 2 1]);
cartesian_grid_xyz1 = cartesian_grid1./sum(sum(sum(cartesian_grid1)));
% add_msg_board(['Grid1-Sum = ' num2str(sum(sum(sum(cartesian_grid_xyz1))))]);

% load dataset 2
data_struct2 = importdata(in_file2,' ',8);
data2 = data_struct2.data';
% description2 = data_struct2.textdata{1};
grid_counts2 = (sscanf(data_struct2.textdata{2},'%*s %*i %*s %*s %*s %i %i %i'))'; % x, y, z convention
origin2 = sscanf(data_struct2.textdata{3},'%*s %f %f %f'); % get origin vector
delta2 = sscanf(data_struct2.textdata{4},'%*s %f %*f %*f'); % read grid spacing
% delta_after_interp = delta_old/(2^ntimes);

% dim_after_interp = arrayfun(@(x) 2^ntimes*x-(2^ntimes-1) ,grid_counts);
cartesian_grid2 = zeros([grid_counts2(3) grid_counts2(2) grid_counts2(1)]);

if size(data2(:,1)) ~= 3
    error('Wrong dimension of input data!');
end
dim_change_offset2 = sum(isnan(data2(:)));

cartesian_grid2(:) = data2(1:(end-dim_change_offset2)); % data must be 3 x n
cartesian_grid2 = permute(cartesian_grid2,[3 2 1]);
cartesian_grid_xyz2 = cartesian_grid2./sum(sum(sum(cartesian_grid2)));
% add_msg_board(['Grid2-Sum = ' num2str(sum(sum(sum(cartesian_grid_xyz2))))]);


% calc distance distribution
non_zero_ind1 = find(cartesian_grid_xyz1 ~= 0);
non_zero_ind2 = find(cartesian_grid_xyz2 ~= 0);
size1 = size(cartesian_grid_xyz1);
size2 = size(cartesian_grid_xyz2);

Rori = origin2 - origin1;
r_hist = linspace(0,200,401)';
y_hist = zeros(size(r_hist));

% tic;
for i = non_zero_ind1'
    [h, k, l] = ind2sub(size1,i);
    r1_vec = delta1.*([h; k; l]-0.5);
    prob1 = cartesian_grid_xyz1(h,k,l);
    for j = non_zero_ind2'
        [m, n, o] = ind2sub(size2, j);
        r2_vec = delta2.*([m; n; o]-0.5);
        dist_here = norm(Rori - r1_vec + r2_vec); % for hist. with r= 0..200 Ang. with 401 points
        hist_ind = round(2*dist_here);
        prob_here = prob1 * cartesian_grid_xyz2(m,n,o);
        y_hist(hist_ind) = y_hist(hist_ind) + prob_here;
    end
end
% toc

rax = r_hist./10; % conv to nm
distr = y_hist;


