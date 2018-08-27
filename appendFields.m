function existStr=appendFields(existStr,sourceStr,keepsource)

% adds all fields from the source structure sourceStr to an existing
% structure existStr without removing any field frim the existing structrue

% keepsource - flag: if missing or 0, sourceStr will be erased; if 1 - will be
% kept.

% !!! tested only for a simple structure so far (no nested structrue for example)
% !!! (280813)

%--------------------------------------------------------------------------
if nargin<2
    keep=0;
else
    keep=keepsource;
end

sourceFields=fieldnames(sourceStr);
[ne,~]=size(sourceFields);
for ii=1:ne
    existStr.(sourceFields{ii,1})=sourceStr.(sourceFields{ii,1});
end

if keep==0
    clear sourceStr;
end

    
    