function newtag=mk_structure_tag
% newtag=mk_structure_tag
%    Makes a new structure tag with four characters
%
% newtag    new structure tag, empty, if now new tag can be found
%
% G. Jeschke, 2010

global model

newtag='+AAA';
c1=0;
while ~isempty(strfind(model.structure_tags,newtag)) && c1<25,
    c1=c1+1;
    newtag(4)=char(double('A')+c1);
end;
c1=0;
while ~isempty(strfind(model.structure_tags,newtag)) && c1<25,
    c1=c1+1;
    newtag(3)=char(double('A')+c1);
end;
c1=0;
while ~isempty(strfind(model.structure_tags,newtag)) && c1<25,
    c1=c1+1;
    newtag(2)=char(double('A')+c1);
end;
if ~isempty(strfind(model.structure_tags,newtag)),
    newtag='';
end;
