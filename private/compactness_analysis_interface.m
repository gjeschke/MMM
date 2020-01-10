function handles = compactness_analysis_interface(handles)

[fname,pname]=uigetfile('*.dat','Load ensemble description from file');

poi = strfind(fname,'.dat');
basname = fname(1:poi-1);
ensemble_pdb = fullfile(pname,strcat(basname,'.pdb'));

hfig = gcf;

if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Ensemble loading cancelled by user');
    return
end

[file_list,p] = rd_ensemble_description(fullfile(pname,fname));

if isempty(file_list)
    add_msg_board('ERROR: Reading of ensemble file failed.');
    return
end

set(hfig,'Pointer','watch');

if ~exist(ensemble_pdb,'file')
    add_msg_board('PDB file of ensemble not found. Assembling from individual files');
    N = length(file_list);
    for k = 1:length(file_list)
        [ccoor,cN,n,resaxis,sequence] = rd_CA_trace(file_list{k});
        if cN > 1
            add_msg_board(sprintf('ERROR: File %s is supposed to be a single model, but contains %i models.',file_list{k},cN));
            return
        end
        if k == 1
            coor = zeros(N*n,3);
            bas = 0;
            n1 = n;
        else
            if n ~= n1
                add_msg_board(sprintf('ERROR: Number %i of CA atoms in file %s does not match number of CA atoms in first file (%i).',n,file_list{k},n1));
                return
            end
        end
        coor(bas+1:bas+n,:) = ccoor;
        bas = bas + n;
    end
else
    [coor,N,n,resaxis,sequence] = rd_CA_trace(ensemble_pdb);
end
set(hfig,'Pointer','arrow');

options.verbose = true;
options.resaxis = resaxis;
[C,R0,nu,raxis,Rg_distr,P,R0_ee,nu_ee] = local_compactness(coor,N,n,p,options);

add_msg_board(sprintf('Radius of gyration fit parameters: R0 = %4.2f Å, nu = %5.3f',R0,nu));
add_msg_board(sprintf('Root mean square segment length fit parameters: R0ee = %4.2f Å, nuee = %5.3f',R0_ee,nu_ee));

x = [resaxis(1) resaxis(end)];
[tags,colors] = mk_color_list('aa_colors.def');
rgb = zeros(length(sequence),3);
for k = 1:length(sequence)
    svgid = tag2id(sequence(k),tags);
    svgcolor = id2tag(svgid,colors);
    rgb(k,:) = svg2rgb(svgcolor);
end

figure;
image(x,x,C,'CDataMapping','scaled');
hold on;
for k = 1:length(sequence)
    plot(resaxis(k),resaxis(1)-2,'.','MarkerSize',10,'Color',rgb(k,:));
    plot(resaxis(1)-2,resaxis(k),'.','MarkerSize',10,'Color',rgb(k,:));
end
curr_axis = gca;
set(curr_axis,'CLim',[-1,1]);
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
axis equal
mymap = ones(51,3);
for k = 1:25
    mymap(k,2:3) = k/26*[1,1];
    mymap(k+26,1:2) = (25-k)/25*[1,1];
end
mymap = flipud(mymap);
colormap(mymap);
title('Relative deviation from random coil radius of gyration');
figure;
image(x,x,P,'CDataMapping','scaled');
hold on;
for k = 1:length(sequence)
    plot(resaxis(k),resaxis(1)-2,'.','MarkerSize',10,'Color',rgb(k,:));
    plot(resaxis(1)-2,resaxis(k),'.','MarkerSize',10,'Color',rgb(k,:));
end
curr_axis = gca;
set(curr_axis,'CLim',[-1,1]);
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
axis equal
colormap(mymap);
title('Relative deviation of root mean square distance');
figure;
[n1,n2] = size(Rg_distr);
plot(1,1,'k.');
plot(n1,n2,'k.');
image([1,n2],[raxis(1),raxis(end)],Rg_distr,'CDataMapping','scaled');
curr_axis = gca;
curr_axis.YDir = 'normal';
colorbar;
axis tight
xlabel('Segment length');
ylabel('Distribution of radius of gyration');
colormap bone
mymap = colormap;
mymap = flipud(mymap);
colormap(mymap);
title('Variation of the radius of gyration');


function [file_list,pop] = rd_ensemble_description(fname)

file_list = cell(1,10000);
pop = zeros(1,10000);
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
                if length(args) > 1
                    pop(poi) = str2double(char(args(2)));
                end
            end
        end
    end
end

fclose(fid);

file_list = file_list(1:poi);
pop = pop(1:poi);
if sum(pop) == 0
    pop = ones(1,poi)/poi;
end