function newstr=compact(str)
% condense all sequences of blanks in an input string to a single blank
% and remove all leading and trailing blanks

newstr='';
str=strtrim(str);
blanks=0;
for k=1:length(str),
    if strcmp(str(k),' '),
        blanks=blanks+1;
    else
        blanks=0;
    end;
    if blanks<2,
        newstr=[newstr str(k)];
    end;
end;