function coor=atom_coor(indices,tag)
%
%
% Returns the coordinates of an atom in a residue, the residue ist
% addressed by indices, the atom by a tag, if the atom does not exist in
% this residue, an empty array is returned
%
% indices   internal indices, only 1:4 are used
% tag       atom tag, as defined in PDB templates
%
% coor      Cartesian coordinates of the atom
%
% G. Jeschke, 2009

global model

coor=[];

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
id=tag2id(tag,info.atom_tags);
if ~isempty(id),
    num=info.atom_numbers{id};
    [mm,nn]=size(num);
    if nn==1, num=[num 1]; end;
    coor=zeros(1,3);
    pop=0;
    for kk=1:mm,
        coor=coor+num(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(num(kk,1),:);
        pop=pop+num(kk,2);
    end;
    coor=coor/pop;
end;