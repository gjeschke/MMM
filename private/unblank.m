function newstr=unblank(str)
% remove all blanks from an input string

newstr='';
for k=1:length(str),
    if ~strcmp(str(k),' '),
        newstr=[newstr str(k)];
    end;
end;
