function [cindices,full_info]=context(address,radius)
% Indices of objects within a context radius of the input object, operates
% on mean coordinates of objects
%
% address   address (if string) or indices (otherwise) of the input object
% radius    context radius
%
% cindices  index matrix of objects within the context radius, objects are
%           of analogous type as input object (chain, residue, atom)
% full_info plain text info (cell string array) on objects within context
%           radius
%
% G. Jeschke, 2009

cindices=[];
full_info={};

if isa(address,'char'),
    [indices,message]=resolve_address(address);
else
    indices=address;
end;

indices=indices(indices>0);

if length(indices)==2, % expand index array to all coordinate sets if object is a chain
    adr=mk_address(indices);
    adr=sprintf('%s{:}',adr);
    indices=resolve_address(adr);
end;

[m,n]=size(indices);

cindices=zeros(5000,6);
distance=zeros(1,5000);
maxn=0;
poi=0;
for k=1:m,
    aindices=indices(k,:);
    [msg,coor1]=get_object(aindices,'xyz');
    coor1=mean(coor1,1);
    adr=mk_address(aindices(1));
    adr=sprintf('%s(:){%i}',adr,aindices(3));
    switch length(aindices)
        case 3,
            iindices=resolve_address(adr);
            comp=3;
        case 4,
            [msg,iindices]=get_object(adr,'children');
            comp=4;
        case {5,6},
            [msg,iindices]=get_object(adr,'descendants');
            comp=6;
    end;
    if iscell(iindices),
        iindices=cat(1,iindices{:});
    end;
    [mm,nn]=size(iindices);
    for kk=1:mm,
        aiindices=iindices(kk,:);
        aiindices=aiindices(aiindices>0);
        if length(aiindices)==comp && sum(aiindices(1:length(aindices))~=aindices),
            [msg,coor2]=get_object(aiindices,'xyz');
            dist=norm(coor1-mean(coor2,1));
            if dist<=radius,
                poi=poi+1;
                nn=length(aiindices);
                cindices(poi,1:nn)=aiindices;
                distance(poi)=dist;
                if nn>maxn, maxn=nn; end;
            end;
        end;
    end;
end;

if poi>0,
    cindices=cindices(1:poi,1:maxn);
else
    cindices=[];
end;

if ~isempty(cindices),
    full_info{1}=sprintf('Objects within %4.1f Å of selected object:',radius);
    for k=1:poi,
        full_info{k+1}=sprintf('%s at %4.1f Å',description(cindices(k,:)),distance(k));
    end;
else
    full_info={'No object on same hierarchy level within context radius'};
end;


