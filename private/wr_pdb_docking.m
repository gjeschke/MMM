function message=wr_pdb_docking(fname,idCode)
% function wr_pdb(fname,idCode)
%
% Writes a PDB file in format V 3.20 of the current structure
% only the 
%   header, (HEADER)
%   title, (TITLE)
%   remarks on format and the originating program (REMARK 4, REMARK 5)
%   sequence information, (SEQRES)
%   sequence modifications, (MODRES)
%   information on nonstandard residues without their names (HET)
%   secondary structure information, this comes from an original PDB
%   structure (if any) or from HELIX, SHEET, LOOP definitions within MMM,
%   sheets are stored only in single-strand format (no registration)
%   definitions within MMM take precedence (HELIX, SHEET)
%   unit cell information (CRYST1), if present in the original file
%   atom coordinates, (ATOM, HETATM)
%   connection information (CONECT)
%   model numbers (MODEL, ENDMDL)
%   master record (MASTER)
% are written
%
% fname     file name for output, extension .pdb is appended, if no
%           extension is present
% idCode    four-letter PDB code or pseudo-PDB code
%
% Ye. Polyhach, G. Jeschke, 2013

global model

% initialize checksum values

% append proper extension, if extension is missing
if isempty(strfind(fname,'.')), fname=strcat(fname,'.pdb'); end;

% yepo adjustment 300813---------------------------------------------------
% snum=model.current_structure;
all_structures_str=textscan(model.structure_tags,'%s','delimiter',':','MultipleDelimsAsOne',1);
all_struct=all_structures_str{1};
for ii=1:numel(all_struct)
    if strcmp(idCode,all_struct{ii})
        snum=ii;
    end
end
%--------------------------------------------------------------------------

snumold = model.current_structure;
model.current_structure = snum;

message = wr_pdb(fname,idCode);

model.current_structure = snumold;
