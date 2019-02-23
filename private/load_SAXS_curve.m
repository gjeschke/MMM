function curve = load_SAXS_curve(fname)

fid = fopen(fname);
if fid==-1
    curve = [];
    add_msg_board('Warning. Loading of SAXS curve failed');
    return;
end
nl=0;
curve = zeros(10000,4);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    ncol = 0;
    %         fprintf(1,'%s\n',tline); % echo for debugging
    if nl > 0 % skip first line
        dataset = str2num(tline);
        ncol = length(dataset);
        curve(nl,1:ncol) = dataset;
    end
    if ncol > 0 || nl == 0
        nl = nl + 1;
    end
end
curve = curve(1:nl-1,:);
% curve(:,1) = 10*curve(:,1);
fclose(fid);
