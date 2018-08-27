function message=add_annotation(indices,type,text,keywords)
% function add_annotation(indices,type,text,keywords)
%
% Add an automatic annotation to the object addressed by indices
% type      predefined type of automatic annotation, can be
%           Alternate       alternate locations defined for this residue
%           Binding         binding sites
%           Coordinating    atom that coordinates a metal atom
%           General         any unspecified automatic annotation
%           Inserted        inserted residue (PDB crap)
%           Metal           metal atom coordinated by...
%           Missing         missing atom coordinates
%           Mutation        mutated residues (MODRES records)
%           Spin            spin label
% text      the annotation text
% keywords  optional cellstring of keywords

global model

mesage.error=0;
message.text='No error.';

if nargin<4,
    keywords={};
end;

[message,annotations]=get_annotations(indices);

if message.error>0 && ~isempty(annotations),
    add_msg_board(sprintf('%s annotation could not be set for %s',type,mk_address(indices)));
    add_msg_board(message.text);
    return
end;

if isempty(annotations),
    annotations.privacy=0;
    annotations.keywords=[];
    annotations.references=[];
    annotations.text{1}='';
end;

if isempty(annotations.privacy),
    annotations.privacy=0;
end;

if ~isempty(keywords),
    register_keyword(keywords,indices);
    for k=1:length(keywords),
       id=tag2id(keywords{k},model.keywords); % find keyword
       if ~isempty(id) && ~sum(annotations.keywords==id),
            n=length(annotations.keywords);
            annotations.keywords(n+1)=id;
       end;
    end;
end;

% Check, whether there is already an automatic annotation page of this type
page=0;
test=strcat('*',type);
for k=1:length(annotations.text),
    old_text=annotations.text{k};
    if ~isempty(old_text),
        old_text=old_text(1,:);
        token=strtok(old_text);
        if strcmp(token,test),
            page=k;
            break;
        end;
    end;
end;

if page==0,
    if isempty(annotations.text{1}),
        page=1;
    else
        page=length(annotations.text)+1;
        annotations.privacy=[annotations.privacy,0];
    end;
    switch type
        case 'Alternate'
            annotations.text{page}='*Alternate atom locations exist for this residue or electron density partially missing:';
        case 'Binding'
            annotations.text{page}='*Binding sites in which this residue is involved:';
        case 'Coordinating'
            annotations.text{page}='*Coordinating the metal atom:';
        case 'General'
            annotations.text{page}='*General automatic annotation:';
        case 'Inserted'
            annotations.text{page}='*Inserted residue (PDB structure with residue number inconsistency):';
        case 'Metal'
            annotations.text{page}='*Metal atom coordinated by atoms:';
        case 'Missing'
            annotations.text{page}='*Missing atom coordinates:';
        case 'Mutation'
            annotations.text{page}='*Mutation at this residue:';
        case 'Spin'
            annotations.text{page}='*Spin labeled residue:';
        otherwise % if the type is not known, make a general annotation
            annotations.text{page}='*General automatic annotation:';
    end;
end;

annotations.text{page}=strvcat(annotations.text{page},text);

message=set_annotations(indices,annotations);

if message.error>0,
    add_msg_board(sprintf('%s annotation could not be set for %s',type,mk_address(indices)));
    add_msg_board(message.text);
    return
end;

function register_keyword(key,indices)

global model

if isempty(key),
    return;
end;

if isempty(indices),
    return;
end;

cindices=zeros(1,6);
cindices(1:length(indices))=indices;

if isfield(model,'keywords'),
   for k=1:length(key),
       id=tag2id(key{k},model.keywords); % find keyword
       if isempty(id), % keyword is new
           model.keywords=sprintf('%s%s:',model.keywords,key{k}); % extend keyword list
           poi=length(model.keys);
           model.keys(poi+1).indices=cindices;
       else
           found=0;
           [m,n]=size(model.keys(id).indices);
           for kk=1:m,
               if sum(abs(model.keys(id).indices(kk,:)-cindices))==0,
                   found=1;
                   break;
               end;
           end;
           if ~found,
               model.keys(id).indices=[model.keys(id).indices; cindices];
           end;
       end;
   end;
else
    model.keywords=sprintf(':%s:',key{1});
    model.keys(1).indices=cindices;
    if length(key)>1,
        for k=2:length(key),
           model.keywords=sprintf('%s%s:',model.keywords,key{k});
           model.keys(k).indices=cindices;
        end;
    end;
end;
