function insert_peptide_section(res1,rese,structure,template_indices,Nanchor_res,Canchor_res)
% inserts a peptide chain section section of peptide chain

% Bfactor is set to twice the maximum Bfactor of existing residues

global model
global residue_defs


snum = template_indices(1);
cnum = template_indices(2);
mnum = template_indices(3);

if ~isempty(Nanchor_res)
    model.structures{snum}(cnum).residues{mnum}.info(Nanchor_res).terminal = 0;
end

Bmax=max(model.info{snum}.B_range);

newresidues=res1:rese;
[ma,~]=size(model.structures{snum}(cnum).xyz{mnum});
[ma0,maxconn]=size(model.structures{snum}(cnum).conn);
newiso=[model.structures{snum}(cnum).isotopes; zeros(20*length(newresidues),2,'single')];
newxyz=[model.structures{snum}(cnum).xyz{mnum}; zeros(20*length(newresidues),3)];
newBfac=[model.structures{snum}(cnum).Bfactor{mnum}, zeros(1,20*length(newresidues))];
newBtens=[model.structures{snum}(cnum).Btensor{mnum}; zeros(20*length(newresidues),6,'int32')];
sequence=model.structures{snum}(cnum).sequence;
restags=model.structures{snum}(cnum).residues{mnum}.residue_tags;

newrnum=zeros(1,length(newresidues));
nr=length(model.structures{snum}(cnum).residues{mnum}.info);
rpoi=0;
for k=1:length(newresidues),
    rnum0=newresidues(k);
    rnum=0;
    for kk=1:length(structure(1).residues{1}.info),
        if structure(1).residues{1}.info(kk).number == rnum0,
            rnum=kk;
            break;
        end;
    end;
    if rnum==0,
        add_msg_board(sprintf('Warning: Residue %i not found in modelled loop',rnum0));
        continue;
    end;
    tag=structure(1).residues{1}.info(rnum).name;
    id=tag2id(tag,upper(residue_defs.restags),residue_defs.single_letter_code);
    sequence(rnum0)=id;
    if mnum==1 && isfield(model.structures{snum}(cnum),'seqexist'),
        model.structures{snum}(cnum).seqexist(rnum0)=1;
    end;
    nr=nr+1;
    rpoi=rpoi+1;
    newrnum(rpoi)=nr;
    % model.structures{snum}(cnum).residues{mnum}.info(nr)=structure(1).residues{1}.info(rnum);
    restag=sprintf('%i:',rnum0);
    restags=strcat(restags,restag);
    pointers=structure(1).residues{1}.info(rnum).atom_numbers;
    if isempty(Canchor_res) && k == length(newresidues),
        model.structures{snum}(cnum).residues{mnum}.info(nr).terminal = 1;
    end;
    for anum=1:length(pointers), % loop over atoms
        pointer=structure(1).residues{1}.info(rnum).atom_numbers{anum};
        [loc,n]=size(pointer);
        for lnum=1:loc, % loop over locations
            poi=pointer(lnum,1); % actual coordinate set number
            ma=ma+1;
            pointer(lnum,1)=ma;
            newiso(ma,:)=structure(1).isotopes(poi,:);
            newxyz(ma,:)=structure(1).xyz{1}(poi,:);
            newBfac(ma)=2*Bmax;
            newBtens(ma,:)=structure(1).Btensor{1}(poi,:);
        end;
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_numbers{anum}=pointer;
        model.structures{snum}(cnum).residues{mnum}.info(nr).name=tag;
        model.structures{snum}(cnum).residues{mnum}.info(nr).type=structure(1).residues{1}.info(rnum).type;
        model.structures{snum}(cnum).residues{mnum}.info(nr).secondary=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).hetflag=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).connected=0;
        model.structures{snum}(cnum).residues{mnum}.info(nr).number=structure(1).residues{1}.info(rnum).number;
        model.structures{snum}(cnum).residues{mnum}.info(nr).atom_tags=structure(1).residues{1}.info(rnum).atom_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).elements=structure(1).residues{1}.info(rnum).elements;
        model.structures{snum}(cnum).residues{mnum}.info(nr).location_tags=structure(1).residues{1}.info(rnum).location_tags;
        model.structures{snum}(cnum).residues{mnum}.info(nr).insertion_code=structure(1).residues{1}.info(rnum).insertion_code;
    end;
end;
newiso=newiso(1:ma,:);
newxyz=newxyz(1:ma,:);
newBfac=newBfac(1:ma);
newBtens=newBtens(1:ma,:);

if mnum==1,
    model.structures{snum}(cnum).sequence=sequence;
    model.structures{snum}(cnum).isotopes=newiso;
    model.structures{snum}(cnum).conn=[model.structures{snum}(cnum).conn; zeros(ma-ma0,maxconn)];
end;

model.structures{snum}(cnum).xyz{mnum}=newxyz;
model.structures{snum}(cnum).Bfactor{mnum}=newBfac;
model.structures{snum}(cnum).Btensor{mnum}=newBtens;
model.structures{snum}(cnum).residues{mnum}.residue_tags=restags;

% make internal bonds in new residues
for k=1:length(newrnum),
    if newrnum(k) > 0,
        model.structures{snum}(cnum)=mk_internal_bonds(model.structures{snum}(cnum),newrnum(k),residue_defs);
    else
        disp('Aber Hallo!');
    end;
end;

% sort residues by number
info=model.structures{snum}(cnum).residues{mnum}.info;
numbers=length(info);
for k=1:length(info),
    numbers(k)=info(k).number;
end;
info0=info;
[sorted,oldnumbers]=sort(numbers);
tags=':';
for k=1:length(oldnumbers),
    info(k)=info0(oldnumbers(k));
    tag=id2tag(oldnumbers(k),model.structures{snum}(cnum).residues{mnum}.residue_tags);
    tags=[tags tag ':'];
end;
model.structures{snum}(cnum).residues{mnum}.residue_tags=tags;
model.structures{snum}(cnum).residues{mnum}.info=info;
