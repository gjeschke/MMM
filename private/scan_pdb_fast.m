function maxat = scan_pdb_fast(fname)
% function protein=scan_pdb_fast(fname)
%
% Finds highest atom number maxat in PDB file
%
% G. Jeschke, 2014

maxat = 10000;
% maxres = 0;

fid=fopen(fname);
if fid==-1,
    return;
end;

maxat2 = 0;
nl=0;
while 1
    tline = fgetl(fid);
%             fprintf(1,'%s\n',tline); % echo for debugging

    nl=nl+1;
    if ~ischar(tline), break, end
    if length(tline)>=6, % catches too short end line of MolProbity files
        record=tline(1:6);
    elseif length(tline) >=3 && strcmpi(tline(1:3),'END');
        record='END   ';
    end;
    switch record
        case 'MODEL '
            if length(tline)>=11,
                curr_model=sscanf(tline(11:end),'%i');
            else
                curr_model=1;
            end;                
            if curr_model>1, 
                break; 
            end;
        case {'ATOM  ','HETATM'}
            maxat2 = maxat2+1;
            atnum=str2double(tline(7:11));
            if atnum > maxat, maxat = atnum; end; 
%             resnum=str2double(tline(23:26));
%             if resnum > maxres, maxres = resnum; end; 
        case 'END   '
            break;
    end;
    record = '';
end

fclose(fid);

if maxat2 > maxat,
    maxat = maxat2;
end;
maxat = 1000*(ceil(maxat/1000)+1);
% maxres = 100*(ceil(maxres/100)+1),

