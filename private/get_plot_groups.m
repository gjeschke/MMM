function plot_groups = get_plot_groups(fname)
% Extracts plot groups from a restraint file

plot_groups = [];
fid=fopen(fname);
if fid == -1
    add_msg_board('ERROR: Restraint file does not exist');
    return;
end

poi = 0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            if k(1)>1
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end
        end
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#')
            switch upper(char(args(2)))
                case 'PLOTGROUP'
                    poi = poi+1;
                    plot_groups(poi).color = svg2rgb(lower(char(args(3))));
                    plot_groups(poi).conformers = zeros(1,1000);
                    poi2 = 0;
                    for karg = 4:length(args)
                        nums = str2num(char(args(karg)));
                        plot_groups(poi).conformers(1+poi2:length(nums)+poi2) = nums;
                        poi2 = poi2 + length(nums);
                    end
                    plot_groups(poi).conformers=  plot_groups(poi).conformers(1:poi2);
            end
        end
    end
end

fclose(fid);