function address=description(indices)
% function address=description(indices)
% 
% clear text description of the object addressed by indices
%
% G. Jeschke, 2009

global model

address='';
if isempty(indices), return; end;

indices=indices(indices>0);

% make structure tag
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==indices(1),
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
address=sprintf('Struct. %s',stag);
if length(indices)>1
    ctag=id2tag(indices(2),model.chain_tags{indices(1)});
    address=sprintf('%s, chain %s',address,ctag);
    if length(indices)>2,
        address=sprintf('%s{%i}',address,indices(3));
        if length(indices)>3,
            info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
            address=sprintf('%s, res. %s-%i',address,info.name,info.number);
            if length(indices)>4,
                atag=id2tag(indices(5),info.atom_tags);
                address=sprintf('%s at. %s',address,atag);
                if length(indices)>5,
                    ltag=id2tag(indices(6),info.location_tags);
                    if ~isempty(ltag) && ~strcmp(ltag,' '),
                        address=sprintf('%s loc. %s',address,ltag);
                    end;
                end;
            end;
        end;
    end;
end;   