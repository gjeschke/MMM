function [message,annotations]=get_annotations(indices)
% returns the annotations for an object with given indices
%
% G. Jeschke, 2009

global model

message.error=0;
message.text='No error';

if isempty(indices),
    annotations=[];
    message.error=5;
    message.text='No object indexed.';
    return
end;

[m,n]=size(indices);
if m>1,
    annotations=[];
    message.error=9;
    message.text='Annotations can be returned only for a single object.';
    return
end;

n=length(model.annotations);
m=length(indices);

if m==0,
    annotations=[];
    message.error=5;
    message.text='No object indexed.';
    return
end;

if n==0,
    annotations=[];
    message.error=6;
    message.text='No annotations.';
    return
end;


found=0;
poi=n+1;
k=1;
if n>0,
    while ~found && k<=n,
        if length(model.annotations(k).indices)>=m,
            if sum(abs(indices-model.annotations(k).indices(1:m)))==0,
                found=1;
                poi=k;
            end;
        end;
        k=k+1;
    end;
end;

if found,
    annotations=model.annotations(poi).info;
else
    annotations=[];
    message.error=4;
    message.text='No annotation for this object';	
end;


