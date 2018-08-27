function [library,NOpos]=find_library_NOpos(adr)
% returns the rotamer library for a labeled residue at a given address
% empty string is returned, if residue is not found or library was not
% stored

global model

library='';
NOpos=[];
indices=resolve_address(adr);
if isempty(indices) || length(indices)~=4,
    return;
end;


nscan=length(model.sites);
for k=1:nscan,
    if isfield(model.sites{k},'library'),
        currlib=model.sites{k}.library;
    else
        currlib='';
    end;    
    for kk=1:length(model.sites{k}.residue),
        diff=sum(abs(indices-model.sites{k}.residue(kk).indices));
        if diff==0,
            library=currlib;
            NOpos=model.sites{k}.residue(kk).NOpos;
        end;
    end;
end;