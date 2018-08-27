function [chi2,outname,status,result,fit] = fit_SAXS_by_crysol(datafile,pdbfile,sm)
%

poi = strfind(pdbfile,'.pdb');
if isempty(poi)
    outname = strcat(pdbfile,'00.fit');
    pdbfile = strcat(pdbfile,'.pdb');
else
    outname = [pdbfile(1:poi-1) '00.fit'];
end;

s=which('crysol.exe');
cmd=[s ' ' pdbfile ' ' datafile ' -cst'];
if exist('sm','var') && ~isempty(sm)
    cmd = [cmd ' -sm ' sprintf('%4.2f',sm)];
end
[status,result]=dos(cmd);

% chi = [];
% poi = strfind(result,'Chi:');
% if ~isempty(poi),
%     rem = textscan(result(poi+4:end),'%s');
%     args = rem{1};
%     chi = str2double(char(args(1)));
% end;
% 
% poi = strfind(result,'Data fit       saved to file');
% if ~isempty(poi),
%     rem = textscan(result(poi+length('Data fit       saved to file'):end),'%s');
%     args = rem{1};
%     outname = char(args(1));
% end;

fit = zeros(10000,4);

chi2 = [];

fid = fopen(outname);
if fid==-1
    fit = [];
    return;
end;
nl=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    %         fprintf(1,'%s\n',tline); % echo for debugging
    if nl > 0 % skip first line
        dataset = str2num(tline);
        ncol = length(dataset);
        fit(nl,1:ncol) = dataset;
    else
        poi = strfind(tline,'Chi^2:');
        if ~isempty(poi)
            rem = textscan(tline(poi+6:end),'%s');
            args = rem{1};
            chi2 = str2double(char(args(1)));
        else
            poi = strfind(tline,'Chi:');
            rem = textscan(tline(poi+4:end),'%s');
            args = rem{1};
            chi = str2double(char(args(1)));
            chi2 = chi^2;
        end;
    end
    nl=nl+1;
end
fit = fit(1:nl-1,:);
if ncol == 4 % remove error column from output of newer CRYSON
    fit = fit(:,[1 2 4]);
end
fclose(fid);

