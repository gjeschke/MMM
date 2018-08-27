function address=mk_residue_address(rnum,cset)
% Makes the address string for a residue in the current structure and chain
%
% rnum      residue number
% cset      (optional) coordinate set number, defaults to 1
% address   address string

global model

if nargin<2, cset=1; end;

snum=model.current_structure;
ids=model.structure_ids;
sid=[];
for k=1:length(ids),
    if ids(k)==snum,
        sid=k;
    end;
end;
stag=id2tag(sid,model.structure_tags);
ctag=model.current_chain;
address=sprintf('[%s](%s){%i}%i',stag,ctag,cset,rnum);