function convert_psipred(iname)
% function convert_psipred(iname)
%
% converts PSIPRED primary output to secondary structure propensity
% restraints in MMM format
%
% iname     full name of the input file, including extension
%           the output name is the input name with .MMM.dat appended
%
% G. Jeschke, 29.8.2017

maxres = 5000;
seq = char(double(' ')*ones(1,maxres));
helix_prop = zeros(1,maxres);
strand_prop = zeros(1,maxres);

fid=fopen(iname);
if fid == -1
    add_msg_board('ERROR: PSIPRED file does not exist');
    return;
end;

maxres = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) break, end
    if ~isempty(tline)
        k = strfind(tline,'#'); % remove comments
        if isempty(k)
            myline = textscan(tline,'%s');
            args=myline{1};
            if length(args) == 6
                resnum = str2double(char(args(1)));
                if resnum > maxres
                    maxres = resnum;
                end
                seq(resnum) = char(args(2));
                helix_prop(resnum) = round(5*str2double(char(args(5))))/5;
                strand_prop(resnum) = round(5*str2double(char(args(6))))/5;
            end
        end
    end
end

fclose(fid);
seq = seq(1:maxres);
helix_prop = helix_prop(1:maxres);
strand_prop = strand_prop(1:maxres);

oname = strcat(iname,'.MMM.dat');

fid = fopen(oname,'wt');
if fid == -1
    add_msg_board('ERROR: Output file could not be created.');
    return;
end;

fprintf(fid,'# SEQUENCE %s\n',seq);
if sum(helix_prop(helix_prop > strand_prop)) > 0
    fprintf(fid,'# APROP\n');
    for k = 1:maxres
        if helix_prop(k) > 0 && helix_prop(k) > strand_prop(k)
            fprintf(fid,'%i %4.1f\n',k,helix_prop(k));
        end
    end
end

if sum(strand_prop(strand_prop > helix_prop)) > 0
    fprintf(fid,'# BPROP\n');
    for k = 1:maxres
        if strand_prop(k) > 0 && strand_prop(k) > helix_prop(k)
            fprintf(fid,'%i %4.1f\n',k,strand_prop(k));
        end
    end
end

fprintf(fid,'# END\n');

fclose(fid);

