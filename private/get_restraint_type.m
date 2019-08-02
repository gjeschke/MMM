function type = get_restraint_type(fname)
% Determines the type of a restraint file

fid=fopen(fname);
if fid == -1
    add_msg_board('ERROR: Restraint file does not exist');
    type = '';
    return;
end

type = 'GENERIC';

mode = 0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline) || mode<0, break, end
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
                case 'RIGIFLEX'
                    type = 'RIGIFLEX';
                case 'FLEX'
                    type = 'FLEX';
                case 'TINKER'
                    type = 'TINKER';
                case 'LOCATE'
                    type = 'LOCATE';
                case 'NETWORK'
                    type = 'NETWORK';
                case 'DOCKING'
                    type = 'DOCKING';
                case 'END'
                    mode = -1;
            end
        end
    end
end

fclose(fid);