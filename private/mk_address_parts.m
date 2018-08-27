function [stag,ctag,modelnum,resnum,icode,atag,ltag,resname]=mk_address_parts(indices)
% Makes the address parts for given indices of a single object
% index arrays are not allowed as input
% 
% all unspecidied parts are empty, for wrong input indices, all parts are
% empty, check isempty(stag) to recognize invalid inputs
%
% indices   indices
% 

global model

stag='';
ctag='';
modelnum=[];
resnum=[];
icode='';
atag='';
ltag='';

if isempty(indices),
    return;
end;

[m,n]=size(indices);

if m>1 || n<1 || n>6,
    return;
end;

cindices=indices;
cindices=cindices(cindices>0);

% make structure tag
ids=model.structure_ids;
sid=[];
for kk=1:length(ids),
    if ids(kk)==cindices(1),
        sid=kk;
    end;
end;
stag=id2tag(sid,model.structure_tags);
if length(cindices)>1
    ctag=id2tag(cindices(2),model.chain_tags{cindices(1)});
    if length(cindices)>2,
        modelnum=cindices(3);
        if length(cindices)>3,
            info=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4));
            if isfield(model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)),'insertion_code'),
                icode=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).insertion_code;
            else
                icode=' ';
            end;
            resnum=info.number;
            resname=info.name;
            if length(cindices)>4,
                atag=id2tag(cindices(5),info.atom_tags);
                if length(cindices)>5,
                    ltag=id2tag(cindices(6),info.location_tags);
                end;
            end;
        end;
    end;
end;

