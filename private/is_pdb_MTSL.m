function labeled = is_pdb_MTSL(pdbid)
% function labeled = is_pdb_MTSL(pdbid)
%
% Determines whether a PDB file contains an MTSL residue (residue code R1A)
% only the SEQRES records are scanned for this
%
% pdbid         PDB identifier, file is downloaded to the temp directory
%               and deleted there after analysis
% labeled       flag, true if MTSL is contained, false otherwise, empty, if
%               no PDB filed with this pdbid was found
%
% G. Jeschke, 2012


labeled = false;

fname = get_pdb_file(pdbid);

if isempty(fname),
    labeled=[];
    return
end;

fid=fopen(fname);
if fid==-1,
    labeled=[];
    return;
end;

nl=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    if length(tline)>=6, % catches too short end line of MolProbity files
        record=tline(1:6);
    else
        record='END   ';
    end;
    switch record
        case 'SEQRES'
            k = strfind(upper(tline),'R1A');
            if ~isempty(k),
                labeled = true;
                break; 
            end
        case 'HETNAM'
            if strcmpi(tline(12:14),'MTN'),
                labeled = true;
                break; 
            end
    end;
    if strcmpi(record,'END   '),
            break;
    end;
end

fclose(fid);

delete(fname);