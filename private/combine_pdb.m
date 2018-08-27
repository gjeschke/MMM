function message = combine_pdb(fname,chainlist)
% message = combine_pdb(fname,chainlist)
%
% Combines chain sections from multiple PDB files into a single PDB file
% the PDB file is minimal (only atom records)
% all files must have the extension .pdb
%
% fname     name of the output PDB file
% chainlist array of structures specifying the chain sections with fields
%           .file   name of the PDB file containing the chain section
%           .chain  chain to be extracted
%           .ochain chain identifier in the output file
%           example .file = 'RNA_sec3_mod1';
%                   .chain = 'A';
%                   .ochain = 'D';
%           inserts chain A from file RNA_sec3_mod1 at the current position
%           into the output file as a section of chain D
%
% G. Jeschke, 16.4.2018

message.error = 0;
message.text = 'OK';

fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='File could not be written';
    return;
end;

for k = 1:length(chainlist)
    iname = strcat(chainlist(k).file,'.pdb');
    ifid=fopen(iname);
    if ifid==-1,
        message.error=2;
        message.text= sprintf('Input file %s could not be opened',iname);
        fclose(fid);
        return;
    end;
    nat = 0;
    while 1
        tline = fgetl(ifid);
        if ~ischar(tline), break, end
        if length(tline)>=6, % catches too short end line of MolProbity files
            record=tline(1:6);
        else
            record='END   ';
        end;
        switch record
            case {'ATOM  ','HETATM'} % these are the only records that are used
                chain_tag=tline(22);
                if chain_tag == chainlist(k).chain
                    nat = nat + 1;
                    tline(22) = chainlist(k).ochain;
                    fprintf(fid,'%s\n',tline);
                end
        end
    end
    fclose(ifid);
    fprintf(1,'%i atoms inserted from file %s chain %s into chain %s\n',chainlist(k).file,chainlist(k).chain,chainlist(k).ochain);
end

fprintf(fid,'END   \n');
fclose(fid);