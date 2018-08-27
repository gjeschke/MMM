function message=consolidate_chain(indices)
% function consolidate_chain(indices)
%
% consolidates the coordinates, Bfactors, and connection records for a
% chain by removing empty coordinate records
%
% indices   1x2 vector of integer, MMM indices addressing the chain
%
% message   error diagnostics, message.error 0: No error. 1: Isotope list
%           inconsistent (nothing changed) 2:: warning that in a chain with
%           several models, some models were inconsistent with model 1,
%           these inconsistent models were removed
%
% G. Jeschke, 2011

global model

message.error=0;
message.text='No error.';

snum=indices(1);
cnum=indices(2);

% check consistency of chain models and create translation tables
models=length(model.structures{snum}(cnum).xyz);
[mi,ni]=size(model.structures{snum}(cnum).isotopes);
for mnum=1:models,
    [ma,na]=size(model.structures{snum}(cnum).xyz{mnum});
    if mi~=ma,
        message.error=1;
        message.text='ERROR: Isotope list is inconsistent with coordinate list.';
        add_msg_board(message.text);
        return
    end;
end;
mnum=1;
translation=zeros(1,ma);
info=model.structures{snum}(cnum).residues{mnum}.info;
newpoi=0;
newxyz=zeros(ma,3);
newBfac=zeros(1,ma);
newBtens=zeros(ma,6,'int32');
maxconn0=model.structures{snum}(cnum).maxconn;
newconn=zeros(ma,maxconn0,'int32');
newiso=zeros(ma,2,'single');
maxconn=0;
for rnum=1:length(info), % loop over residues
    pointers=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers;
    for anum=1:length(pointers), % loop over atoms
        pointer=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum};
        [loc,n]=size(pointer);
        for lnum=1:loc, % loop over locations
            poi=pointer(lnum,1); % actual coordinate set number
            newpoi=newpoi+1;
            pointer(lnum,1)=newpoi;
            translation(poi)=newpoi;
            newiso(newpoi,:)=model.structures{snum}(cnum).isotopes(poi,:);
            newxyz(newpoi,:)=model.structures{snum}(cnum).xyz{mnum}(poi,:);
            newBfac(newpoi)=model.structures{snum}(cnum).Bfactor{mnum}(poi);
            newBtens(newpoi,:)=model.structures{snum}(cnum).Btensor{mnum}(poi,:);
            newconn(newpoi,:)=model.structures{snum}(cnum).conn(poi,:);
            mb=sum(newconn(newpoi,:)>0);
            if mb>maxconn, maxconn=mb; end;
        end;
        model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum}=pointer;
    end;
end;
newiso=newiso(1:newpoi,:);
newxyz=newxyz(1:newpoi,:);
newBfac=newBfac(1:newpoi);
newBtens=newBtens(1:newpoi,:);
newconn=newconn(1:newpoi,1:maxconn);
for anum=1:newpoi,
    for bnum=1:maxconn,
        if newconn(anum,bnum)>0,
            newconn(anum,bnum)=translation(newconn(anum,bnum));
        end;
    end;
end;
model.structures{snum}(cnum).isotopes=newiso;
model.structures{snum}(cnum).conn=newconn;
model.structures{snum}(cnum).xyz{mnum}=newxyz;
model.structures{snum}(cnum).Bfactor{mnum}=newBfac;
model.structures{snum}(cnum).Btensor{mnum}=newBtens;
maxconn1=maxconn;
newpoi1=newpoi;
newiso1=newiso;
newconn1=newconn;
mpoi=1;
if models>1,
    for mnum=2:models,
        info=model.structures{snum}(cnum).residues{mnum}.info;
        translation=zeros(1,ma);
        newpoi=0;
        newxyz=zeros(ma,3);
        newBfac=zeros(1,ma);
        newBtens=zeros(ma,6,'int32');
        newconn=zeros(ma,maxconn0,'int32');
        newiso=zeros(ma,2,'single');
        maxconn=0;
        for rnum=1:length(info), % loop over residues
            pointers=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers;
            for anum=1:length(pointers), % loop over atoms
                pointer=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum};
                [loc,n]=size(pointer);
                for lnum=1:loc, % loop over locations
                    poi=pointer(lnum,1); % actual coordinate set number
                    newpoi=newpoi+1;
                    pointer(lnum,1)=newpoi;
                    translation(poi)=newpoi;
                    newxyz(newpoi,:)=model.structures{snum}(cnum).xyz{mnum}(poi,:);
                    newBfac(newpoi)=model.structures{snum}(cnum).Bfactor{mnum}(poi);
                    newBtens(newpoi,:)=model.structures{snum}(cnum).Btensor{mnum}(poi,:);
                    mb=sum(newconn(newpoi,:)>0);
                    if mb>maxconn, maxconn=mb; end;
                end;
                model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum}=pointer;
            end;
        end;
        if newpoi~=newpoi1,
            add_msg_board(sprintf('Warning: Atom number in chain model %i inconsistent with chain model 1.',mnum));
            add_msg_board(sprintf('Removing chain model %i.',mnum));
            message.error=2;
            message.text='Warning: Inconsistent model(s) removed.';
            continue;
        end;
        newxyz=newxyz(1:newpoi,:);
        newBfac=newBfac(1:newpoi);
        newBtens=newBtens(1:newpoi,:);
        model.structures{snum}(cnum).xyz{mnum}=newxyz;
        model.structures{snum}(cnum).Bfactor{mnum}=newBfac;
        model.structures{snum}(cnum).Btensor{mnum}=newBtens;        
    end;
end;


