function address=mk_address(indices,name)
% Makes the address string for given indices,
% for index arrays, a cell string is returned
%
% indices   indices
% name      optional flag that determines whether the residue name should
%           be appended as a comment, default 0 (no residue name), 1 name
%           is appended
% address   address string
% 
% If all indices are zero, an empty string is returned

global model

if nargin<2,
    name=0;
end;

address='';
if isempty(indices),
    return;
end;

[m,n]=size(indices);

for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    
    if isempty(cindices),
        address = '';
        return;
    end;

    % make structure tag
    ids=model.structure_ids;
    sid=[];
    for kk=1:length(ids),
        if ids(kk)==cindices(1),
            sid=kk;
        end;
    end;
    stag=id2tag(sid,model.structure_tags);
    address=sprintf('[%s]',stag);
    if length(cindices)>1
        ctag=id2tag(cindices(2),model.chain_tags{cindices(1)});
        address=sprintf('%s(%s)',address,ctag);
        if length(cindices)>2,
            address=sprintf('%s{%i}',address,cindices(3));
            if length(cindices)>3,
                info=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4));
                if isfield(model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)),'insertion_code'),
                    icode=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).insertion_code;
                else
                    icode='';
                end;
                address=strtrim(sprintf('%s%i%s',address,info.number,icode));
                if length(cindices)>4,
                    atag=id2tag(cindices(5),info.atom_tags);
                    address=sprintf('%s.%s',address,atag);
                    if length(cindices)>5,
                        ltag=id2tag(cindices(6),info.location_tags);
                        address=sprintf('%s:%s',address,ltag);
                    end;
                end;
            end;
        end;
    end;
    if length(cindices)>=4 && name,
        resname=info.name;
        address=[address '; ' resname];
    end;
    address0{k}=address;
end;

if m==1,
    address=address0{1};
else
    address=address0;
end;