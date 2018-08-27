function ids = tag2ids(tag,tags)
% function ids = tag2ids(tag,tags)
%
% Returns the identification codes corresponding to a string tag by
% comparison with a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if tag is not found in tags or if the position is larger than the number
% of elements in codelist, an empty id is given back
%
% tag       (string) tag to be found, example: 'Be'
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:Be:'
%
% ids       ids telling the positions in the tag list, for example above:
%           [4,8]
%
% G. Jeschke, 24.3.2017

delimiter = ':'; % colon as default delimiter

etag = [delimiter tag delimiter];
positions = strfind(tags,etag);
ids = zeros(1,length(positions));
for k = 1:length(positions)
    ids(k) = 1+sum(find(tags==delimiter) < positions(k));
end;

