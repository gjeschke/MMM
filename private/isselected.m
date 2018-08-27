function flag=isselected(identifier)
% function flag=isselected(identifier)
%
% Determines whether an object is selected
%
% identifier    address of index vector to identify object, if identifier
%               has more than one row, output flag is true only if  a l l
%               objects are selected
% flag          true (1) if object is selected, false (0) if not
%
% G. Jeschke, 2010

flag=false;

selection=resolve_address('*');

if isempty(selection),
    return;
end;

if isa(identifier,'char'),
    [indices,message]=resolve_address(identifier);
    if message.error,
        return;
    end;
else
    indices=identifier;
end;

[m,n]=size(indices);
[ms,ns]=size(selection);

flag=true;
for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    found=false;
    for kk=1:ms,
        sindices=selection(kk,:);
        sindices=sindices(sindices>0);
        if length(sindices)==length(cindices),
            if sum(abs(sindices-cindices))==0,
                found=true;
            end;
        end;
    end;
    if ~found, flag=false; break; end;
end;
