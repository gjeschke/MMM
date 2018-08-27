function [coor,elements,conn,NO]=get_virtual_label(labeled_site,rotamer,get_H)
% function [coor,elements,conn]=get_virtual_label(labeled_site,rotamer,get_H)
%
% provides coordinates, element identifiers, and connection information for
% the rotamer with number rotamer in the site scan result labeled_site
% the coordinates conform to the current backbone coordinates of the
% labeled residue
%
% Input:
% labeled_site  rotamer ensemble descrition for the labeled site, 
%               corresponds to model.sites{ks}(kc).residue(kr)
% rotamer       number of the rotamer to be returned
% get_H         optional flag that determines whether hydrogen coordinates
%               and connections are to be returned, defaults to true
%
% Output:
% coor          Cartesian coordinates of the rotamer
% elements      element vector of the atoms
% conn          connection table, relative to coordinate matrix
% NO            indices into coordinate array coor for the NO group
%
% G. Jeschke, 2010

global label_defs
global model

if nargin<3,
    get_H=true;
end;

label=labeled_site.label;
rotamers=labeled_site.rotamers;
indices=labeled_site.indices;
adr=mk_address(indices);

coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};

% get information on the label type
labnum=tag2id(label,label_defs.restags);
atags=label_defs.residues(labnum).atoms;
info.name=label_defs.residues(labnum).tc;
NR_tag=id2tag(1,label_defs.residues(labnum).frame);
OR_tag=id2tag(2,label_defs.residues(labnum).frame);
N_label=tag2id(NR_tag,atags);
O_label=tag2id(OR_tag,atags);

atags0=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_tags;
N=tag2id('N',atags0);
CA=tag2id('CA',atags0);
C=tag2id('C',atags0);

% get further information on the current residue
anum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{N};
apoi=anum(:,1)';
Ncoor=mean(coors(apoi,:),1);
anum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{CA};
apoi=anum(:,1)';
CAcoor=mean(coors(apoi,:),1);
anum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{C};
apoi=anum(:,1)';
Ccoor=mean(coors(apoi,:),1);


N=tag2id('N',atags);
CA=tag2id('CA',atags);
C=tag2id('C',atags);

coor1=[Ncoor;CAcoor;Ccoor];

conn=label_defs.residues(labnum).conn;

coor=rotamers(rotamer).coor;
coor2=[coor(N,:);coor(CA,:);coor(C,:)];
x=coor1(1,:)-coor1(2,:); % x axis is along C_alpha-N bond
x=x/norm(x);    % unit vector along x
yp=coor1(3,:)-coor1(2,:); % y axis is in the plane spanned by x axis and C-Ca bond
yp=yp/norm(yp);
z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z=z/norm(z);
y=cross_rowvec(z,x); % real (corrected) y axis 
dircos=[x;y;z];
Rp=dircos; % rotation matrix for conversion to standard frame

newatoms=0;
[m,n]=size(coor);
newpoi=zeros(1,m);
elements=label_defs.residues(labnum).elements;
for k=1:m,
    if ~isnan(coor(k,1)),
        if get_H || elements(k)~=1,
            newatoms=newatoms+1;
            newpoi(newatoms)=k;
            if k==N_label,
                N_new=newatoms;
            end;
            if k==O_label,
                O_new=newatoms;
            end;
        end;
    end;
end;
[mc,nc]=size(conn);
nconny=zeros(newatoms,nc);
ncoor=zeros(newatoms,3);
elements=zeros(1,newatoms);
for k=1:newatoms,
    ncoor(k,:)=coor(newpoi(k),:);
    elements(k)=label_defs.residues(labnum).elements(newpoi(k));
    c1=[];
    for kk=1:nc,
        if conn(newpoi(k),kk)~=0,
            c0=conn(newpoi(k),kk);
            for kkk=1:length(newpoi),
                if newpoi(kkk)==c0, c1=[c1 kkk]; end;
            end;
        end;
        nconny(k,1:length(c1))=c1;
    end;
end;
conn=nconny;

[m,n] = size(coor);
for k = 1:m,
    newvec = coor(k,:)*Rp + CAcoor;
    coor(k,:) = newvec;
end;


NO=[N_new,O_new];