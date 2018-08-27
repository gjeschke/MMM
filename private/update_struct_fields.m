function struct1=update_struct_fields(struct1,struct2)
% Updates all fields in struct1 that also exist in struct2 with the
% values they have in struct2, 
% thus allows assignments between partially mismatching struct variables,
% for instance in arays of struct variables
% fields that exist in struct1, but not struct2, retain their old values
% fields that exist in struct2, but not struct1, are ignored
% field name matching is case-sensitive
%
% Assignment of the result to struct1 is mandatory
%
% G. Jeschke, 2010

names1=fieldnames(struct1);
names2=fieldnames(struct2);

for k=1:length(names2),
    tf=strcmp(names2{k},names1);
    if sum(tf),
        value=getfield(struct2,names2{k});
        struct1=setfield(struct1,names2{k},value);
    end;
end;