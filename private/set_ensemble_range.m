function set_ensemble_range(snum)
% function set_ensemble_range(snum)
%
% determines the range of C_alpha r.m.s.d. of individual residues in an
% ensemble of structures (coordinate sets on hierarchy level 3, i.e.
% different models of the same macromolecular complex)
% 
% the result is stored in field model.info{snum}.ensemble_range in the form
% of a vector [min_rmsd,max_rmsd]
% if only a single model exists (for all chains)
% model.info{snum}.ensemble_range=[] is set, use
% isempty(model.info{snum}.ensemble_range) to test for this condition
%
% snum      number of the structure

global model

min_rmsd=1e6;
max_rmsd=-1e6;
valid=false;
chains=length(model.structures{snum}(:));
for cnum=1:chains,
    models=length(model.structures{snum}(cnum).residues);
    if models>1,
        valid=true;
        coor=zeros(models,3);
        residues=length(model.structures{snum}(cnum).residues{1}.info);
        for rnum=1:residues,
            present=0;
            for mnum=1:models,
                adr=mk_address([snum,cnum,mnum,rnum]);
                indices=resolve_address([adr '.CA']);
                [msg,ccoor]=get_object(indices,'coor');
                if ~isempty(ccoor),
                    coor(mnum,:)=ccoor;
                    present=present+1;
                end;
            end;
            if present>0,
                mcoor=sum(coor,1)/present;
                diff=zeros(1,models);
                for mnum=1:models,
                    diff(mnum)=norm(coor(mnum,:)-mcoor);
                end;
                rmsd=sqrt(sum(diff.^2/(models-1)));
                if rmsd<min_rmsd, min_rmsd=rmsd; end;
                if rmsd>max_rmsd, max_rmsd=rmsd; end;
            end;
        end;
    end;
end;
if valid,
    model.info{snum}.ensemble_range=[min_rmsd,max_rmsd];
    add_msg_board(sprintf('Ensemble Calpha r.ms.d. ranges between %5.2f and %5.2f Å',min_rmsd,max_rmsd));
else
    model.info{snum}.ensemble_range=[];
end;
