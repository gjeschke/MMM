function name=basname(FileName)
% Strips the last extension (after a .) from a file name, if there is one

points=strfind(FileName,'.');
if isempty(points)
    name=FileName;
else
    name=FileName(1:points(end)-1);
end;
