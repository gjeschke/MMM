function rgb=color_by_charge(slc)
% RGB color triple corresponding to the charge of a residue encoded by
% single-letter code slc
%
% G. Jeschke, 2009

global residue_defs

rgb=[0.75,0.75,0.75]; % default for no charge or unknown charge, light grey

resnum=findstr(slc,residue_defs.single_letter_code);
if ~isempty(resnum) && resnum<=length(residue_defs.residues),
    info=residue_defs.residues(resnum);
    switch info.charge,
        case 2, rgb=[0,0,0.5];
        case 1, rgb=[0.0,0.0,1];
        case -1, rgb=[1,0,0];
        case -2, rgb=[0.5,0,0];
    end;
end;
