function repair_sidechains(snum)
% function repair_sidechains(snum)
% Tries to repair sidechains of incompletely defined residues via SCWRL4
% this is successful for all residues whose backbone atoms are defined
%

global model

indices=model.info{snum}.incomplete;
[m,n]=size(indices);
errors=zeros(1,m);

for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    if length(cindices)==4,
        tlc=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).name;
        message=mutate(cindices,tlc);
        errors(k)=message.error;
    else
        errors(k)=10;
    end;
end;
still_incomplete=indices;
failed=0;
for k=1:m,
    adr=mk_address(indices(k,:),true);
    if errors(k),
        failed=failed+1;
        still_incomplete(failed,:)=indices(k,:);
        add_msg_board(sprintf('Residue %s could not be repaired.',adr));
    else
        add_msg_board(sprintf('Sidechain of residue %s was repaired.',adr));
    end;
end;
still_incomplete=still_incomplete(1:failed,:);
model.info{snum}.incomplete=still_incomplete;
if failed==0,
    add_msg_board('All residues could be repaired.');
else
    add_msg_board(sprintf('Warning: %i residues are still incomplete.',failed));
end;
