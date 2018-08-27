function msg=reassign_secondary(indices)
% function reassign_secondary(indices)
% Reassigns secondary structure for a given chain by DSSP information, if
% this is present, this can be restricted to a chain, chain model, or a
% stretch of residues or can encompass a whole structure
%
% indices   indices of the objects that should be reassigned, [1xn] vector,
%           or [mxn] index array
%           length 1    whole structure is reassigned
%           length 2    all models of one chain are reassigned
%           length 3    one models of one chain is reassigned
%           [mxn], n=4  all addressed residues are reassigned, in this case
%                       helix and strand ranges of all affected chains are
%                       redetermined
%
% G. Jeschke, 2010

global model

msg.error=0;
msg.text='No error.';

if model.locked,
    msg.error=1;
    msg.text='Secondary structure is locked.';
    add_msg_board('ERROR: Secondary structure is locked. Use unlock button first.');
    return
end;

num_dssp=0;
num_changed=0;
[m,n]=size(indices);
chains=zeros(1000,2);
numc=0;
for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    n=length(cindices);
    if n<4,
        [msg,indices1]=get_object(cindices,'descendants');
    else
        indices1=cindices;
    end;
    [m1,n1]=size(indices1);
    for kk=1:m1,
        cindices1=(indices1(kk,:));
        cindices1=cindices1(cindices1>0);
        % Determine which chains are affected
        newchain=true;
        for kkk=1:numc,
            if length(cindices1)>1 && sum(abs(cindices1(1:2)-chains(kkk,:)))==0,
                newchain=false;
                break;
            end;
        end;
        if newchain,
            numc=numc+1;
            chains(numc,:)=cindices1(1:2);
        end;
        if length(cindices1)==4, % this is a residue
            info=model.structures{cindices1(1)}(cindices1(2)).residues{cindices1(3)}.info(cindices1(4));
            oldsec=info.secondary;
            if isfield(info,'dssp'), % DSSP information present
                num_dssp=num_dssp+1;
                dssp_type=char(info.dssp);
                switch dssp_type,
                    case {'H','G','I'}
                        info.secondary=1;
                    case 'E'
                        info.secondary=2;
                    otherwise
                        info.secondary=0;
                end;
                if info.secondary~=oldsec,
                    num_changed=num_changed+1;
                end;
                model.structures{cindices1(1)}(cindices1(2)).residues{cindices1(3)}.info(cindices1(4))=info;
            end;
        end;
    end;
end;
add_msg_board(sprintf('For %i residues DSSP information was present.',num_dssp));
add_msg_board(sprintf('For %i residues secondary structure information changed.',num_changed));

chains=chains(1:numc,:);
add_msg_board(sprintf('Helices and strands are now redetermined for %i chain(s).',numc));

for k=1:numc,
    snum=chains(k,1);
    cnum=chains(k,2);
    message=pdb_secondary(snum,cnum);
end;
