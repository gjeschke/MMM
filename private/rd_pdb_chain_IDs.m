function chain_IDs = rd_pdb_chain_IDs(fname)
% chain_IDs = rd_pdb_chain_IDs(fname)
%
% Finds all chain IDs in a PDB file that contain ATOM records
%
% G. Jeschke, 12.06.2020

chain_IDs = '$$$$$$$$$$$$$$$$$$$$$$$$$$';

if ~contains(fname,'.')
    fname = strcat(fname,'.pdb');
end

ifid=fopen(fname);
if ifid==-1
    warning('rd_pdb_chain_IDs:Input file could not be opened');
    return;
end

nc = 0;
while 1
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if length(tline) >= 54 
        if strcmpi(tline(1:4),'ATOM')
            cid = tline(22);
            if ~contains(chain_IDs,cid)
                nc = nc + 1;
                chain_IDs(nc) = cid;
            end
        end
    end
end
fclose(ifid);

chain_IDs = chain_IDs(1:nc);