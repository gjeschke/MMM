function [protein,chain_tags,chain_ids]=test_rd_pdb_coarse(pdbid)

[protein,chain_tags,chain_ids]=rd_pdb_coarse(pdbid);

fid=fopen('test_output.txt','wt');
for k=1:length(chain_ids),
    tag=id2tag(k,chain_tags);
    fprintf(fid,'[%s](%s) Database reference: %s\n',pdbid,tag,protein(k).dbref);
    for kk=1:length(protein(k).resnum),
        fprintf(fid,'%s%i CA %8.3f %8.3f %8.3f\n',protein(k).restype(kk),protein(k).resnum(kk),protein(k).Ca(kk,:));
        if ~strcmpi(protein(k).restype(kk),'G'),
            fprintf(fid,'%s%i CB %8.3f %8.3f %8.3f\n',protein(k).restype(kk),protein(k).resnum(kk),protein(k).Cb(kk,:));
        end;
    end;
end;
fclose(fid);