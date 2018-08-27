function replace_sidechain(structure,indices,mutated,modnum,first)
% Replaces a residue in an MMM-internal structure addressed by indices by
% the corresponding residue from a structure description read in by
% rd_pdb.m
%
% structure     structure with the new sidechain read in by rd_prb.m, must
%               the first model (coordinate set) is used
% indices       indices of the residue to be replaced in the MMM-internal
%               structure
% mutated       residue index, coressponding to indices(4) in the structure
%               read in by rd_pdb (chain index matches the one in indices)
% modnum        model in the internal structure in which the residue is to
%               be replaced, substitutes for indices(3), optional, if
%               missing, indices(3) is used
% first         optional flag that tells that this is done in the first of
%               several models, defaults to true
%
% G. Jeschke, 2010

global model
global chemistry

if nargin<4,
    modnum=indices(3);
end;

if nargin<5,
    first=true;
end;

% disp(indices);
info=model.structures{indices(1)}(indices(2)).residues{modnum}.info(indices(4));
N=tag2id('N',info.atom_tags);
atnum_N=info.atom_numbers{N};
atnum_N=atnum_N(1,1);
C=tag2id('C',info.atom_tags);
atnum_C=info.atom_numbers{C};
atnum_C=atnum_C(1,1);
Bfactor=model.structures{indices(1)}(indices(2)).Bfactor{modnum};
Bfac0=Bfactor(atnum_N);
Btensor=model.structures{indices(1)}(indices(2)).Btensor{modnum};
Btens0=Btensor(atnum_N,:);
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
atoms=model.structures{indices(1)}(indices(2)).atoms{modnum};
conn=model.structures{indices(1)}(indices(2)).conn;
[mconn,maxconn]=size(conn);
xyz=model.structures{indices(1)}(indices(2)).xyz{modnum};
[mxyz,nxyz]=size(xyz);conn_N=conn(atnum_N,:);
conn_C=conn(atnum_C,:);
res_atnums=zeros(1,200);
poi=0;
for k=1:length(info.atom_numbers), % extract all atom numbers of old residue
    catom=info.atom_numbers{k};
    catnum=catom(:,1);
    res_atnums(poi+1:poi+length(catnum))=catnum;
    poi=poi+length(catnum);
end;
res_atnums=res_atnums(1:poi);
for k=1:length(conn_N), % remove all residue-internal bonds
    if sum(find(conn_N(k)==res_atnums)),
        conn_N(k)=0;
    end;
end;
conn_N=conn_N(conn_N~=0);
for k=1:length(conn_C), % remove all residue-internal bonds
    if sum(find(conn_C(k)==res_atnums)),
        conn_C(k)=0;
    end;
end;
conn_C=conn_C(conn_C~=0);


newxyz=structure(indices(2)).xyz{1};
newinfo=structure(indices(2)).residues{1}.info(mutated);
res_atnums=zeros(1,200);
poi=0;
for k=1:length(newinfo.atom_numbers), % extract all atom numbers of old residue
    catom=newinfo.atom_numbers{k};
    catnum=catom(:,1);
    res_atnums(poi+1:poi+length(catnum))=catnum;
    poi=poi+length(catnum);
end;
res_atnums=res_atnums(1:poi);
for k=1:length(newinfo.atom_numbers), % extract all atom numbers of old residue
    newinfo.atom_numbers{k}=k+mxyz;
end;
newinfo.number=info.number;
newconn=structure(indices(2)).conn;
newbonds=newconn(res_atnums,:);
[mb,newmaxconn]=size(newbonds);
if maxconn>newmaxconn,
    newbonds=[newbonds,zeros(mb,maxconn-newmaxconn)];
end;
if newmaxconn>maxconn,
    conn=[conn,zeros(mconn,newmaxconn-maxconn)];
end;
N=tag2id('N',newinfo.atom_tags);
newconn_N=newbonds(N,:);
C=tag2id('C',newinfo.atom_tags);
newconn_C=newbonds(C,:);
for k=1:length(newconn_N), % remove all residue-external bonds
    if newconn_N(k)~=0,
        if sum(find(newconn_N(k)==res_atnums))==0,
            if ~isempty(conn_N),
                newconn_N(k)=conn_N(1);
            else
                newconn_N(k)=0;
            end;
        else
            whichone=find(newconn_N(k)==res_atnums,1);
            newconn_N(k)=whichone+mxyz;
        end;
    end;
end;
pp=0;
for k=1:length(newconn_N),
    if newconn_N(k)~=0,
        pp=pp+1;
        newconn_N(pp)=newconn_N(k);
    end;
end;
for k=pp+1:length(newconn_N),
    newconn_N(k)=0;
end;
for k=1:length(newconn_C), % remove all residue-external bonds
    if newconn_C(k)~=0
        if sum(find(newconn_C(k)==res_atnums))==0,
            if ~isempty(conn_C),
                newconn_C(k)=conn_C(1);
            else
                newconn_C(k)=0;
            end;
        else
            whichone=find(newconn_C(k)==res_atnums,1);
            newconn_C(k)=whichone+mxyz;
        end;
    end;
end;
pp=0;
for k=1:length(newconn_C),
    if newconn_C(k)~=0,
        pp=pp+1;
        newconn_C(pp)=newconn_C(k);
    end;
end;
for k=pp+1:length(newconn_C),
    newconn_C(k)=0;
end;
for k=1:mb,
    if k~=N && k~=C,
        for kk=1:newmaxconn,
            if newbonds(k,kk)~=0
                whichone=find(res_atnums==newbonds(k,kk),1);
                newbonds(k,kk)=whichone+mxyz;
            end;
        end;
    end;
end;
newbonds(N,:)=newconn_N;
newbonds(C,:)=newconn_C;
newcoor=newxyz(res_atnums,:);
xyz=[xyz;newcoor];
conn=[conn;newbonds];
atoms=atoms+length(res_atnums);
Bfactor=[Bfactor,Bfac0*ones(1,length(res_atnums))];
extBtens0=repmat(Btens0,length(res_atnums),1);
Btensor=[Btensor;extBtens0];
newinfo.secondary=info.secondary;
ma=length(newinfo.elements);
newisotopes=zeros(ma,2);
newisotopes(:,1)=newinfo.elements';
for k=1:ma,
    newisotopes(k,2)=chemistry.pse(newisotopes(k,1)).mass;
end;
model.structures{indices(1)}(indices(2)).Bfactor{modnum}=Bfactor;
model.structures{indices(1)}(indices(2)).Btensor{modnum}=Btensor;
model.structures{indices(1)}(indices(2)).atoms{modnum}=atoms;
if first,
    model.structures{indices(1)}(indices(2)).conn=conn;
    model.structures{indices(1)}(indices(2)).isotopes=[isotopes;newisotopes];
end;
model.structures{indices(1)}(indices(2)).xyz{modnum}=xyz;
model.structures{indices(1)}(indices(2)).residues{modnum}.info(indices(4))=update_struct_fields(model.structures{indices(1)}(indices(2)).residues{modnum}.info(indices(4)),newinfo);

