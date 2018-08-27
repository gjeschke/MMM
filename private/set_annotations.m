function message=set_annotations(indices,annotations)
% assigns new annotations to an object with given indices
%
% annotations are a structure with mandatory fields
%    .privacy    privacy level 0: world, 1: group, 2: owner, vector of privacy levels
%                corresponding to annotation pages (elements of cellstring text)
%    .keywords   integer id list of keywords for this 
%	             annotation, can be empty, but must be present
%    .references reference list for this annotation, can be empty, but must be
%   	         present (integer id's to reference list)
%    .text       cellstring with annotation text  
%
% G. Jeschke, 2009

global model

message.error=0;
message.text='No error';

[m,n]=size(indices);
if m>1,
    message.error=9;
    message.text='Annotations can be set only for a single object.';
    return
end;

m=length(indices);
n=length(model.annotations);

if m==0,
    message.error=5;
    message.text='No object indexed.';
    return
end;

aindices=zeros(1,6);

aindices(1:m)=indices;
found=0;
poi=n+1;
k=1;
if n>0,
    while ~found && k<=n,
        if sum(abs(indices-model.annotations(k).indices(1:m)))==0,
            found=1;
            poi=k;
        end;
        k=k+1;
    end;
end;

model.annotations(poi).indices=aindices;
model.annotations(poi).info=annotations;

