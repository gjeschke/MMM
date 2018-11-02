function [message,argout]=set_residue(indices,property,argin)
% function [message,argout]=set_residue(indices,property,argin)
%
% Sets properties of a residue in MMM
%
% indices   index vector that identifies the residue in the model
% property  property to be set, given as a string, available properties and
%           corresponding set functions are defined at the beginning of the
%           function source code in the switch board
% argin     further input arguments given as a structure argin{:}
%
% message   error message structure with fields .error and .text,
%           .error=0 indicates no error
% argout    further output arguments given as structure argout{n}
%
% G. Jeschke, 2009-2013

argout={};

message.error=0;
message.text='';

if nargin<3
    argin{1}='';
end

switch property
    case 'annotations'
        message=set_annotations(indices,argin{1});
    case 'color'
        message=color_residue(indices,argin{1},0);
    case 'colorscheme'
        message=colorscheme_residue(indices,argin);
    case 'hide'
        message=hide_residue(indices);
    case 'label'
        message=label_residue(indices,argin{1},argin{2},argin{3},argin{4});
    case 'ribboncolor'
        message=color_residue(indices,argin{1},1);
    case 'secondary'
        message=secondary_residue(indices,argin{1});
    case 'show'
        message=show_residue(indices,argin);
    case 'transform'
        message=transform_residue(indices,argin);
    case 'transparency'
        message=transparency_residue(indices,argin{1});
    case 'uncolor'
        message=uncolor_residue(indices);
%     case 'unlabel'
%         message=unlabel_residue(indices);
    case 'untransparency'
        message=untransparency_residue(indices);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function message=unlabel_residue(indices)

global model

message.error=0;
message.text='';

if isfield(model,'native'),
    found=0;
    for k=1:length(model.native),
        if sum(model.native(k).indices-indices)==0,
            found=k;
        end;
    end;
end;

if ~found,
    message.error=1;
    message.text='No information on native sidegroup found.';
    adr=mk_address(indices);
    add_msg_board(sprintf('### Warning ### Unlabeling impossible.'));
    add_msg_board(sprintf('No information on native sidegroup of residue %s found.',adr));
    return
end;

model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4))=model.native(found).info;

function message=label_residue(indices,label,rotamers,attach_frame,exclude)

global label_defs

message.error=0;
message.text='';

labnum=tag2id(label,label_defs.restags);
if isfield(label_defs.residues(labnum),'class') && ~isempty(label_defs.residues(labnum).class),
    class = label_defs.residues(labnum).class; 
else
    class = 'nitroxide';
end;
if isfield(label_defs.residues(labnum),'attachment') && ~isempty(label_defs.residues(labnum).attachment)
    attachment = label_defs.residues(labnum).attachment;
else
    attachment = 'peptide';
end;
if strcmp(attachment,'peptide'), 
    % if strcmp(class,'nitroxide')
    message=label_peptide_residue(indices,label,rotamers);
else
    % if strcmp(class,'nitroxide')
    message=label_nucleotide_residue(indices,label,rotamers,attach_frame,exclude);
end;

function message=label_peptide_residue(indices,label,rotamers)
% Attaches a specified label to a residue, whereas a right-handed
% attachment frame is superimposed
%
% label         label type, three-letter identifier of PDB template or
%               pseudo-PDB template
% rotamers      list of rotamers and their populations
%               rotamers(i).coor  coordinates
%               rotamers(i).pop   populations

global model
global label_defs
global chemistry

loctags0=':A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:';
loctags=':';

rr=length(rotamers);

message.error=0;
message.text='';


% get information on the current residue
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));

if strcmpi(info.name,label),
    adr=mk_address(indices);
    add_msg_board(sprintf('Residue %s is already labeled with label %s. Nothing changed',adr,label));
    message.error=1;
    message.text='Warning! Same label was already attached.';
    return
end;

if isfield(model,'native'),
    poi=length(model.native);
    model.native(poi+1).indices=indices;
    model.native(poi+1).info=info;
else
    model.native(1).indices=indices;
    model.native(1).info=info;
end;

isotopes=model.structures{indices(1)}(indices(2)).isotopes;
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
Bfactors0=model.structures{indices(1)}(indices(2)).Bfactor{indices(3)};
Btensors0=model.structures{indices(1)}(indices(2)).Btensor{indices(3)};
conny=model.structures{indices(1)}(indices(2)).conn;

[atoms,n]=size(coors);
[isonum,n]=size(isotopes);
if atoms<isonum, % workaround for chains with several coordinate sets
    new_model=0;
else
    new_model=1;
end;


% get information on the label type
labnum=tag2id(label,label_defs.restags);
atags=label_defs.residues(labnum).atoms;
NR_tag=id2tag(1,label_defs.residues(labnum).frame);
OR_tag=id2tag(2,label_defs.residues(labnum).frame);
info.name=label_defs.residues(labnum).tc;
mylabel.adr=mk_address(indices);
mylabel.name=label_defs.residues(labnum).short_name;
mylabel.tc=label_defs.residues(labnum).tc;
mylabel.NOpos=zeros(length(rotamers),4);
conn=label_defs.residues(labnum).conn;
[atx,maxconn0]=size(conny);
model.structures{indices(1)}(indices(2)).maxconn=maxconn0;
[mconn,nc]=size(conn);
if nc>model.structures{indices(1)}(indices(2)).maxconn,
    model.structures{indices(1)}(indices(2)).maxconn=nc;
end;



% get further information on the current residue
N=tag2id('N',atags);
adr1=[mylabel.adr '.N'];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
Ncoor=coors(atom(1,1),:);
CA=tag2id('CA',atags);
adr1=[mylabel.adr '.CA'];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
CAcoor=coors(atom(1,1),:);
C=tag2id('C',atags);
adr1=[mylabel.adr '.C'];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
Ccoor=coors(atom(1,1),:);
O=tag2id('O',atags);
adr1=[mylabel.adr '.O'];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
Ocoor=coors(atom(1,1),:);

coor1=[Ncoor;CAcoor;Ccoor];


atom_numbers={};
for rotamer=1:length(rotamers),
    locid=mod(rotamer,26);
    if locid==0, locid=26; end;
    locmid=floor((rotamer-1)/26);
    if locmid>0,
        loctag=strcat(id2tag(locmid,loctags0),id2tag(locid,loctags0));
    else
        loctag=id2tag(locid,loctags0);
    end;
    loctags=sprintf('%s%s:',loctags,loctag);
    coor=rotamers(rotamer).coor;
    [m,n]=size(coor);
    ncoor=coor;
    % made consistent with get_rotamers 3.9.2013
    x=coor1(1,:)-coor1(2,:); % x axis is along C_alpha-N bond
    x=x/norm(x);    % unit vector along x
    yp=coor1(3,:)-coor1(2,:); % y axis is in the plane spanned by x axis and C-Ca bond
    yp=yp/norm(yp);
    z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
    z=z/norm(z);
    y=cross_rowvec(z,x); % real (corrected) y axis 
    dircos=[x;y;z];
    Rp=dircos; % rotation matrix for conversion to standard frame
    newtags=':';
    newatoms=0;
    newpoi=zeros(1,m);
    for k=1:m,
        atag=id2tag(k,atags);
        if ~strcmpi(atag,'H2') && ~strcmpi(atag,'OXT') && ~strcmpi(atag,'HXT') && ~strcmpi(atag,'H'), % these atoms are not attached
            ncoor(k,:)=coor(k,:)*Rp + CAcoor; % 3.9.2013
            if ~isnan(ncoor(k,1)),
                newtags=[newtags atag ':'];
                newatoms=newatoms+1;
                newpoi(newatoms)=k;
                if k==N, N_new=newatoms; end;
                if k==C, C_new=newatoms; end;
            end;
        end;
        if strcmpi(atag,NR_tag) || strcmpi(atag,OR_tag),
            mylabel.NOpos(rotamer,1:3)=mylabel.NOpos(rotamer,1:3)+ncoor(k,:)/2;
            mylabel.NOpos(rotamer,4)=rotamers(rotamer).pop;
        end;
    end;
    ncoor(N,1:3)=coor1(1,1:3);
    ncoor(CA,1:3)=coor1(2,1:3);
    ncoor(C,1:3)=coor1(3,1:3);
    ncoor(O,1:3)=Ocoor;
%     disp('New coordinates');
%     disp(ncoor);
    nisotopes=zeros(atoms+newatoms,2);
    ncoors=zeros(atoms+newatoms,3);
    Bfactors=zeros(1,atoms+newatoms);
    Btensors=zeros(atoms+newatoms,6,'int32');
    nconny=zeros(atoms+newatoms,model.structures{indices(1)}(indices(2)).maxconn);
    if new_model,
        nisotopes(1:atoms,:)=isotopes;
        nconny(1:atoms,1:maxconn0)=conny;
    else
        nisotopes=isotopes;
        nconny=conny;
    end;
    ncoors(1:atoms,:)=coors;
    Bfactors(1:atoms)=Bfactors0;
    Btensors(1:atoms,:)=Btensors0;
    elements=zeros(1,newatoms);
    for k=1:newatoms,
        if rotamer==1, % initialize atom number table
            atnum=zeros(length(rotamers),2);
            atom_numbers{k}=atnum;
        end;
        atnum=atom_numbers{k};
        ncoors(atoms+k,:)=ncoor(newpoi(k),:);
        elements(k)=label_defs.residues(labnum).elements(newpoi(k));
        mass=chemistry.pse(elements(k)).mass;
        nisotopes(atoms+k,1)=elements(k);
        nisotopes(atoms+k,2)=mass;
        c1=[];
        for kk=1:nc,
            if conn(newpoi(k),kk)~=0,
                c0=conn(newpoi(k),kk);
                for kkk=1:length(newpoi),
                    if newpoi(kkk)==c0, c1=[c1 kkk+atoms]; end;
                end;
            end;
            nconny(atoms+k,1:length(c1))=c1;
        end;
        atnum(rotamer,1)=atoms+k;
        atnum(rotamer,2)=rotamers(rotamer).pop;
        atnum(rotamer,3)=rotamer;
        if length(rotamers)==1,
            atom_numbers{k}=atnum(1,1);
        else
            atom_numbers{k}=atnum;
        end;
    end;
    
    if rotamer==1,
        % Make backbone bonds to previous and next residue
        N_new=N_new+atoms;
        C_new=C_new+atoms;
        basadr=mk_address(indices(1:3)); % address of chain model
        mynum=info.number;
        if mynum>1,
            prev_adr=[basadr sprintf('%i',mynum-1) '.C']; % address of peptide carbonyl of previous residue
            [pindices,message]=resolve_address(prev_adr);
            if ~isempty(pindices),
                [msg,p_atnum]=get_atom(pindices,'number');
                if ~isempty(p_atnum),
                    newconn=nconny(N_new,:);
                    newconn=newconn(newconn>0);
                    newconn=[newconn p_atnum(1,1)]; % first location is connected
                    if new_model,
                        nconny(N_new,1:length(newconn))=newconn;
                    end;
                end;
            end;
            next_adr=[basadr sprintf('%i',mynum+1) '.N']; % address of peptide nitrogen of next residue
            [nindices,message]=resolve_address(next_adr);
            if ~isempty(nindices),
                [msg,n_atnum]=get_atom(nindices,'number');
                if ~isempty(n_atnum),
                    newconn=nconny(C_new,:);
                    newconn=newconn(newconn>0);
                    newconn=[newconn n_atnum(1,1)]; % first location is connected
                    if new_model,
                        nconny(C_new,1:length(newconn))=newconn;
                    end;
                end;
            end;
        end;
    end;
    
    % update total number of atoms
    atoms=atoms+newatoms;
    if new_model,
        isotopes=nisotopes(1:atoms,:);
        conny=nconny(1:atoms,1:model.structures{indices(1)}(indices(2)).maxconn);
    end;
    coors=ncoors(1:atoms,:);
    Bfactors0=Bfactors(1:atoms);
    Btensors0=Btensors(1:atoms,:);
    maxconn0=model.structures{indices(1)}(indices(2)).maxconn;
end;

% store information on the new residue definition
info.hetflag=1;
info.connected=0;
info.type=1;
info.atom_tags=newtags;
info.atom_numbers=atom_numbers;
info.elements=elements;
info.location_tags=loctags;

model.structures{indices(1)}(indices(2)).isotopes=nisotopes;
model.structures{indices(1)}(indices(2)).xyz{indices(3)}=ncoors;
model.structures{indices(1)}(indices(2)).Bfactor{indices(3)}=Bfactors;
model.structures{indices(1)}(indices(2)).Btensor{indices(3)}=Btensors;
model.structures{indices(1)}(indices(2)).conn=nconny;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4))=info;

if ~isfield(model,'labels')
    model.labels(1)=mylabel;
else
    nlabels=length(model.labels);
    model.labels(nlabels+1)=mylabel;
end;

rnum2=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
seq=model.structures{indices(1)}(indices(2)).sequence;
if indices(4)>=1 && indices(4)<=length(seq),
    seq(rnum2)='?';
end;
model.structures{indices(1)}(indices(2)).sequence=seq;

% show_residue(indices,'label_frame');
% highlight_selection;

function message=label_nucleotide_residue(indices,label,rotamers,attach_frame,attach_exclude)
% Attaches a specified label to a residue, whereas a right-handed
% attachment frame is superimposed
%
% label         label type, three-letter identifier of PDB template or
%               pseudo-PDB template
% rotamers      list of rotamers and their populations
%               rotamers(i).coor  coordinates
%               rotamers(i).pop   populations
% attach_frame  PDB atom tags (MMM tag string) for atoms on the yp axis, x
%               axis, and the origin of the attachment frame
% exclude       atoms to be excluded from attachment

global model
global label_defs
global chemistry
global rotamer_libraries

loctags0=':A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:';
loctags=':';

rr=length(rotamers);

message.error=0;
message.text='';

% get attachment frame atom tags
ypax_name = id2tag(1,attach_frame);
xax_name = id2tag(2,attach_frame);
orig_name = id2tag(3,attach_frame);


% get information on the current residue
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));

if strcmpi(info.name,label),
    adr=mk_address(indices);
    add_msg_board(sprintf('Residue %s is already labeled with label %s. Nothing changed',adr,label));
    message.error=1;
    message.text='Warning! Same label was already attached.';
    return
end;

if isfield(model,'native'),
    poi=length(model.native);
    model.native(poi+1).indices=indices;
    model.native(poi+1).info=info;
else
    model.native(1).indices=indices;
    model.native(1).info=info;
end;

isotopes=model.structures{indices(1)}(indices(2)).isotopes;
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
Bfactors0=model.structures{indices(1)}(indices(2)).Bfactor{indices(3)};
Btensors0=model.structures{indices(1)}(indices(2)).Btensor{indices(3)};
conny=model.structures{indices(1)}(indices(2)).conn;

[atoms,n]=size(coors);
[isonum,n]=size(isotopes);
if atoms<isonum, % workaround for chains with several coordinate sets
    new_model=0;
else
    new_model=1;
end;


% get information on the label type
labnum=tag2id(label,label_defs.restags);
atags=label_defs.residues(labnum).atoms;
NR_tag=id2tag(1,label_defs.residues(labnum).frame);
OR_tag=id2tag(2,label_defs.residues(labnum).frame);
info.name=label_defs.residues(labnum).tc;
mylabel.adr=mk_address(indices);
mylabel.name=label_defs.residues(labnum).short_name;
mylabel.tc=label_defs.residues(labnum).tc;
mylabel.NOpos=zeros(length(rotamers),4);
conn=label_defs.residues(labnum).conn;
[atx,maxconn0]=size(conny);
model.structures{indices(1)}(indices(2)).maxconn=maxconn0;
[mconn,nc]=size(conn);
if nc>model.structures{indices(1)}(indices(2)).maxconn,
    model.structures{indices(1)}(indices(2)).maxconn=nc;
end;

% make attachment frame
adr1=[mylabel.adr '.' xax_name];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
xax_coor=coors(atom(1,1),:);
adr1=[mylabel.adr '.' orig_name];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
orig_coor=coors(atom(1,1),:);
adr1=[mylabel.adr '.' ypax_name];
cindices=resolve_address(adr1);
atom=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers{cindices(5)};
ypax_coor=coors(atom(1,1),:);

exclude = ':';
% determine atoms that to be removed
if ~isempty(label_defs.residues(labnum).remove),
    exclude = label_defs.residues(labnum).remove;
end;

% determine the atoms of the existing residue, which must be included
old_tags = model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_tags;
old_atoms = 0;
old_numbers = model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).atom_numbers;
elements = model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.info(cindices(4)).elements;
old_elements = zeros(1,length(old_numbers));
atom_numbers = old_numbers;
newtags=':';
for k = 1:length(old_numbers),
    ctag = id2tag(k,old_tags); % atom tag of existing atom
    id = tag2id(ctag,exclude); % is this atom among the ones to be excluded?
    if isempty(id), % if not, put it into atom number list
        old_atoms = old_atoms + 1;
        atom_numbers{old_atoms} = old_numbers{k};
        old_elements(old_atoms) = elements(k);
        newtags = [newtags ctag ':'];
    end;
end;
atom_numbers = atom_numbers(1:old_atoms);
old_elements = old_elements(1:old_atoms);

for rotamer=1:length(rotamers),
    locid=mod(rotamer,26);
    if locid==0, locid=26; end;
    locmid=floor((rotamer-1)/26);
    if locmid>0,
        loctag=strcat(id2tag(locmid,loctags0),id2tag(locid,loctags0));
    else
        loctag=id2tag(locid,loctags0);
    end;
    loctags=sprintf('%s%s:',loctags,loctag);
    coor=rotamers(rotamer).coor;
    [m,n]=size(coor);
    ncoor=coor;
    % made consistent with get_rotamers 3.9.2013
    x=xax_coor - orig_coor; % x axis is along C_alpha-N bond
    x=x/norm(x);    % unit vector along x
    yp=ypax_coor - orig_coor; % y axis is in the plane spanned by x axis and C-Ca bond
    yp=yp/norm(yp);
    z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
    z=z/norm(z);
    y=cross_rowvec(z,x); % real (corrected) y axis 
    dircos=[x;y;z];
    Rp=dircos; % rotation matrix for conversion to standard frame
    newatoms=0;
    newpoi=zeros(1,m);
    for k=1:m,
        atag=id2tag(k,atags);
        exid = tag2id(atag,attach_exclude);
        if isempty(exid), % the atom should be attached
            ncoor(k,:)=coor(k,:)*Rp + orig_coor; % 3.9.2013
            if ~isnan(ncoor(k,1)),
                newtags=[newtags atag ':'];
                newatoms=newatoms+1;
                newpoi(newatoms)=k;
            end;
            if strcmpi(atag,NR_tag) || strcmpi(atag,OR_tag),
                mylabel.NOpos(rotamer,1:3)=mylabel.NOpos(rotamer,1:3)+ncoor(k,:)/2;
                mylabel.NOpos(rotamer,4)=rotamers(rotamer).pop;
            end;
        end;
    end;
    elements=zeros(1,newatoms + old_atoms);
    elements(1:old_atoms) = old_elements;
    nisotopes=zeros(atoms+newatoms,2);
    ncoors=zeros(atoms+newatoms,3);
    Bfactors=zeros(1,atoms+newatoms);
    Btensors=zeros(atoms+newatoms,6,'int32');
    nconny=zeros(atoms+newatoms,model.structures{indices(1)}(indices(2)).maxconn);
    if new_model,
        nisotopes(1:atoms,:)=isotopes;
        nconny(1:atoms,1:maxconn0)=conny;
    else
        nisotopes=isotopes;
        nconny=conny;
    end;
    ncoors(1:atoms,:)=coors;
    Bfactors(1:atoms)=Bfactors0;
    Btensors(1:atoms,:)=Btensors0;
    for k=1:newatoms,
        if rotamer==1, % initialize atom number table
            atnum=zeros(length(rotamers),2);
            atom_numbers{k + old_atoms} = atnum;
        end;
        atnum=atom_numbers{k + old_atoms};
        ncoors(atoms+k,:)=ncoor(newpoi(k),:);
        elements(k + old_atoms)=label_defs.residues(labnum).elements(newpoi(k));
        mass=chemistry.pse(elements(k + old_atoms)).mass;
        nisotopes(atoms+k,1)=elements(k + old_atoms);
        nisotopes(atoms+k,2)=mass;
        c1=[];
        for kk=1:nc,
            if conn(newpoi(k),kk)~=0,
                c0=conn(newpoi(k),kk);
                for kkk=1:length(newpoi),
                    if newpoi(kkk)==c0, c1=[c1 kkk+atoms]; end;
                end;
            end;
            nconny(atoms+k,1:length(c1))=c1;
        end;
        atnum(rotamer,1)=atoms+k;
        atnum(rotamer,2)=rotamers(rotamer).pop;
        atnum(rotamer,3)=rotamer;
        if length(rotamers)==1,
            atom_numbers{k + old_atoms}=atnum(1,1);
        else
            atom_numbers{k + old_atoms}=atnum;
        end;
    end;
    
    % update total number of atoms
    atoms=atoms+newatoms;
    if new_model,
        isotopes=nisotopes(1:atoms,:);
        conny=nconny(1:atoms,1:model.structures{indices(1)}(indices(2)).maxconn);
    end;
    coors=ncoors(1:atoms,:);
    Bfactors0=Bfactors(1:atoms);
    Btensors0=Btensors(1:atoms,:);
    maxconn0=model.structures{indices(1)}(indices(2)).maxconn;
end;

% store information on the new residue definition
info.hetflag=1;
info.connected=0;
info.type=2;
info.atom_tags = newtags;
info.atom_numbers=atom_numbers;
info.elements=elements;
info.location_tags=loctags;

model.structures{indices(1)}(indices(2)).isotopes=nisotopes;
model.structures{indices(1)}(indices(2)).xyz{indices(3)}=ncoors;
model.structures{indices(1)}(indices(2)).Bfactor{indices(3)}=Bfactors;
model.structures{indices(1)}(indices(2)).Btensor{indices(3)}=Btensors;
model.structures{indices(1)}(indices(2)).conn=nconny;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4))=info;

if ~isfield(model,'labels')
    model.labels(1)=mylabel;
else
    nlabels=length(model.labels);
    model.labels(nlabels+1)=mylabel;
end;

rnum2=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
seq=model.structures{indices(1)}(indices(2)).sequence;
if indices(4)>=1 && indices(4)<=length(seq),
    seq(rnum2)='?';
end;
model.structures{indices(1)}(indices(2)).sequence=seq;

function message=secondary_residue(indices,type)
% Hides a residue by deleting graphics objects
%

global model

message.error=0;
message.text='';

switch type
    case 'loop'
        sec=0;
    case 'helix'
        sec=1;
    case 'strand'
        sec=2;
end;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).secondary=sec;



function message=color_residue(indices,rgb,onlyribbon)
% Sets color of residue graphics
%

global model

message.error=0;
message.text='';

if nargin<3,
    onlyribbon=0;
end;

graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;

if ~onlyribbon,
    label_gobjects=[];
    if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'label_graphics'),
        label_graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics;
        if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics),
            gobjects=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.objects;
            if ~isempty(gobjects),
                label_gobjects=[label_gobjects gobjects];
            end;
        end;
    end;
end;

if ~isempty(graphics),
    if ~isempty(graphics.objects), 

        if graphics.mode==1,
            prop='Color';
        else
            prop='FaceColor';
        end;

        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                set(graphics.objects(k),prop,rgb);
            end;
        end;
        graphics.color=[rgb;graphics.color];
        model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=graphics;
    end;
end;
if ~onlyribbon,
    if ~isempty(label_gobjects),
        if isfield(label_graphics,'color'),
            label_graphics.color=[rgb;label_graphics.color];
        else
            label_graphics.color=rgb;
        end;
        for k=1:length(label_gobjects),
            if ishandle(label_gobjects(k)) && label_gobjects(k)~=0,
                h = set(label_gobjects(k));
                if isfield(h,'FaceColor'),
                    set(label_gobjects(k),'FaceColor',rgb);
                end;
            end;
        end;
    end;
    model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics=label_graphics;
end;

if ~onlyribbon,
    if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics') && length(model.structures{indices(1)}(indices(2)).atom_graphics)>=indices(3),
        for anum=1:length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers)
            atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{anum};
            [m,n]=size(atom);
            for k=1:m,
                if ~isempty(model.structures{indices(1)}(indices(2)).atom_graphics),
                    if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
                        graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
                        if isempty(graphics), continue; end;
                        if isempty(graphics.objects), continue; end;
                        for kk=1:length(graphics.objects),
                            if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                                if isprop(graphics.objects(kk),'Color'),
                                    set(graphics.objects(kk),'Color',rgb);
                                elseif isprop(graphics.objects(kk),'FaceColor'),
                                    set(graphics.objects(kk),'FaceColor',rgb);
                                end;
                            end;
                        end;
                        graphics.color=[rgb;graphics.color];
                        model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=graphics;
                    end;
                end;
            end;
        end;
    end;
end;

function message=uncolor_residue(indices)
% Resets color of residue graphics to previous color
%

global model

message.error=0;
message.text='';

[msg,info_text]=get_residue(indices,'info');
atom_adress=[info_text{1} '.:'];
set_object(atom_adress,'uncolor');

graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;


if ~isempty(graphics),
    if ~isempty(graphics.objects), 

        [m,n]=size(graphics.color);
        if m<2,
            message.error=3;
            message.text='No previous color available';
            return
        end;

        graphics.color=graphics.color(2:m,:);
        rgb=graphics.color(1,:);

        if graphics.mode==1,
            prop='Color';
        else
            prop='FaceColor';
        end;


        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                set(graphics.objects(k),prop,rgb);
            end;
        end;
        model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=graphics;
    end;
end;

[msg,info_text]=get_residue(indices,'info');
atom_adress=[info_text{1} '.:'];
set_object(atom_adress,'uncolor');

function message=colorscheme_residue(indices,argin)
% Sets color of residue graphics
%

global model
global graph_settings
global residue_defs

% parameters for colorscheme difference
mindiff=0.2; % minimum per residue r.m.s.d. (Å) for reporting a difference in the message board
maxdiff=3; % per residue r.m.s.d. (Å) corresponding to upper end of color scale (red)
diffgrades=20; % number of shadings (color grades) in the color scale
report_missing=true; % flag that indicates whether residues missing in the 
                     % reference structure are reported in message board
report_nonmatch=true; % flag that indicates whether residues with non-matching atom number in
                      % the reference structure are reported in message board

message.error=0;
message.text='';

scheme=argin{1};

if length(argin)>1, % treat address argument for difference color scheme
    [stag,ctag,modelnum,resnum]=mk_address_parts(indices);
    adr=argin{2};
    dindices0=resolve_address(adr);
    [stag2,ctag2,modelnum2,resnum2]=mk_address_parts(dindices0);
    if ~isempty(resnum2),
        dindices=resolve_address(sprintf('[%s](%s){%i}%i',stag2,ctag2,modelnum2,resnum2));
    elseif ~isempty(modelnum2),
        dindices=resolve_address(sprintf('[%s](%s){%i}%i',stag2,ctag2,modelnum2,resnum));
    elseif ~isempty(ctag2),
        dindices=resolve_address(sprintf('[%s](%s)%i',stag2,ctag2,resnum));
    elseif ~isempty(stag2),
        dindices=resolve_address(sprintf('[%s]%i',stag2,resnum));
    else
        dindices=[];
    end;
else
    dindices=[];
end;

if length(argin)>2,
    arg3=str2double(argin{3});
    if ~isempty(arg3) && ~isnan(arg3) && arg3>0,
        maxdiff=arg3;
    end;
end;

graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;

label_graphics=[];
if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'label_graphics'),
    if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics),
        gobjects=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.objects;
        if ~isempty(gobjects),
            label_graphics=[label_graphics gobjects];
        end;
    end;
end;

known=0; % flag for known scheme

switch scheme
    case 'chain'
        chains=length(model.structures{indices(1)});
        rgb=color_grade(indices(2),chains);
        known=1;
    case 'secondary'
        sec_type=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).secondary;
        rgb=graph_settings.loop_color;
        if sec_type==1,
            rgb=graph_settings.helix_color;
        end;
        if sec_type==2,
            rgb=graph_settings.sheet_color;
        end;
        known=1;
    case 'sequence'
        num1=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(1).number;
        nume=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(end).number;
        residues=nume-num1+1;
        residues2=length(model.structures{indices(1)}(indices(2)).sequence);
        if residues2<residues,
            residues=residues2-num1+1;
        end;
        rgb=color_grade(indices(4)-num1+1,residues);
        known=1;
    case 'Bfactor'
        [message,Bfactors]=get_residue(indices,'Bfactor');
        Bfactor=sqrt(mean(Bfactors));
        rgb=color_grade(Bfactor-graph_settings.Bmin+1,graph_settings.Bmax-graph_settings.Bmin+1);
        known=1;
    case 'Bfactor_tight'
        [message,Bfactors]=get_residue(indices,'Bfactor');
        Bfactor=sqrt(mean(Bfactors));
        Bmin=floor(sqrt(model.info{indices(1)}.B_range(1)));
        Bmax=ceil(sqrt(model.info{indices(1)}.B_range(2)));
        rgb=color_grade(Bfactor-Bmin+1,Bmax-Bmin+1);
        known=1;
    case 'charge'
        num=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
        if num>0 && num<=length(model.structures{indices(1)}(indices(2)).sequence),
            slc=model.structures{indices(1)}(indices(2)).sequence(num);
        else
            slc='?';
        end;
        rgb=color_by_charge(slc);
        known=1;
    case 'hydropathy'
        sequence=model.structures{indices(1)}(indices(2)).sequence;
        residue=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
        if residue<=length(sequence),
            resnum=findstr(sequence(residue),residue_defs.single_letter_code);
        else
            resnum=[];
        end;
        if ~isempty(resnum),
            hydropathy=residue_defs.residues(resnum).hydropathy;
            hmin=10*graph_settings.min_hydropathy;
            hmax=10*graph_settings.max_hydropathy;
            rgb=color_grade(10*hydropathy-hmin+1,hmax-hmin+1);
        else
            rgb=[0.5,0.5,0.5];
        end;
        known=1;
    case 'helix_propensity'
        sequence=model.structures{indices(1)}(indices(2)).sequence;
        residue=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
        if residue<=length(sequence),
            resnum=findstr(sequence(residue),residue_defs.single_letter_code);
        else
            resnum=[];
        end;
        if ~isempty(resnum),
            hp=residue_defs.residues(resnum).helix_propensity;
            hmin=floor(10*sqrt(graph_settings.min_helix_propensity));
            hmax=ceil(10*sqrt(graph_settings.max_helix_propensity));
            rgb=color_grade(10*sqrt(hp)-hmin+1,hmax-hmin+1);
        else
            rgb=[0.5,0.5,0.5];
        end;
        known=1;
    case 'ensemble'
        models=length(model.structures{indices(1)}(indices(2)).residues);
        if models>1,
            coor=zeros(models,3);
            present=0;
            for mnum=1:models,
                cindices=indices;
                cindices(3)=mnum;
                adr=mk_address(cindices);
                nindices=resolve_address([adr '.CA']);
                [msg,ccoor]=get_object(nindices,'coor');
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
                rmsd_min=floor(10*model.info{indices(1)}.ensemble_range(1));
                rmsd_max=ceil(10*model.info{indices(1)}.ensemble_range(2));
                rgb=color_grade(10*rmsd-rmsd_min+1,rmsd_max-rmsd_min+1);
            else
                rgb=[0.5,0.5,0.5];
            end;
        else
            rgb=[0.5,0.5,0.5];
        end;
        known=1;
    case 'difference'
        if ~isempty(dindices),
            anum1=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
            atnums1=zeros(1,length(anum1));
            for k=1:length(anum1),
                atn=anum1{k};
                atnums1(k)=atn(1,1);
            end;
            xyz1=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
            coor1=xyz1(atnums1,:);
            [m1,n]=size(coor1);
            anum2=model.structures{dindices(1)}(dindices(2)).residues{dindices(3)}.info(dindices(4)).atom_numbers;
            atnums2=zeros(1,length(anum2));
            for k=1:length(anum2),
                atn=anum2{k};
                atnums2(k)=atn(1,1);
            end;
            xyz2=model.structures{dindices(1)}(dindices(2)).xyz{dindices(3)};
            coor2=xyz2(atnums2,:);
            [m2,n]=size(coor2);
            if m1==m2,
                rms=sqrt(sum(sum((coor1-coor2).^2))/m1);
                adr1=mk_address(indices,true);
                adr2=mk_address(dindices,true);
                if rms>mindiff,
                    add_msg_board(sprintf('Residue %s coordinate difference to residue %s is %5.1f Å',adr1,adr2,rms));
                end;
                dgrade=1+round((diffgrades-1)*rms/maxdiff);
                if dgrade>diffgrades, dgrade=diffgrades; end;
                rgb=color_grade(dgrade,diffgrades);
            else
                if report_nonmatch,
                    adr1=mk_address(indices,true);
                    adr2=mk_address(dindices,true);
                    add_msg_board(sprintf('Residuse %s and %s have different numbers of atoms.',adr1,adr2));
                end;
                rgb=[0.5,0.5,0.5];
            end;
        else
            if report_missing,
                adr1=mk_address(indices,true);
                add_msg_board(sprintf('For residue %s, the corresponding residue is missing in the reference structure.',adr1));
            end;
            rgb=[0.5,0.5,0.5];
        end;
        known=1;
end;

if known,
    if ~isempty(graphics),
        if ~isempty(graphics.objects), 

            if graphics.mode==1,
                prop='Color';
            else
                prop='FaceColor';
            end;

            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    set(graphics.objects(k),prop,rgb);
                end;
            end;
            graphics.color=[rgb;graphics.color];
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=graphics;
        end;
    end;
    if ~isempty(label_graphics), 
        for k=1:length(label_graphics),
            if ishandle(label_graphics(k)) && label_graphics(k)~=0,
                set(label_graphics(k),'FaceColor',rgb);
            end;
        end;
    end;
else
    msg.error=1;
    msg.text='Unknown color scheme';
end;

% [msg,info_text]=get_residue(indices,'info');
% atom_adress=[info_text{1} '.:'];
% set_object(atom_adress,'colorscheme',{scheme});

if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics')
    for anum=1:length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers)
        atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{anum};
        [m,n]=size(atom);
        for k=1:m,
            if iscell(model.structures{indices(1)}(indices(2)).atom_graphics) && length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
                graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
                if isempty(graphics), continue; end;
                if isempty(graphics.objects), continue; end;
                if graphics.mode==1,
                    prop='Color';
                else
                    prop='FaceColor';
                end;
                for kk=1:length(graphics.objects),
                    if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                        set(graphics.objects(kk),prop,rgb);
                    end;
                end;
                graphics.color=[rgb;graphics.color];
                model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1))=graphics;
            end;
        end;
    end;
end;

function message=transform_residue(indices,matrices)
% Coordinate transformation defined by an affine 4x4 transformation matrix
% or a cell array of such matrices

[message,newindices]=get_residue(indices,'children');
[m,n]=size(newindices);
for k=1:m
    message=set_atom(newindices(k,:),'transform',matrices);
end

function message=transparency_residue(indices,alpha)
% Sets transparency of residue graphics
%

global model

message.error=0;
message.text='';

if alpha<0, 
    ralpha=1;
else
    ralpha=alpha;
end; % default transparency for residues is opaque


graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;

if ~isempty(graphics),
    if ~isempty(graphics.objects), 

        if graphics.mode==1,
            message.error=3;
            message.text='No transparency for line objects.';
            return
        end;

        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                set(graphics.objects(k),'FaceAlpha',ralpha);
            end;
        end;
        graphics.opaque=[ralpha;graphics.opaque];
        model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=graphics;
    end;
end;

[msg,info_text]=get_residue(indices,'info');
atom_adress=[info_text{1} '.:'];
set_object(atom_adress,'transparency',{alpha});

function message=untransparency_residue(indices)
% Resets transparency of residue graphics to previous value
%

global model

message.error=0;
message.text='';

[msg,info_text]=get_residue(indices,'info');
atom_adress=[info_text{1} '.:'];
set_object(atom_adress,'untransparency');

graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;


if ~isempty(graphics),
    if ~isempty(graphics.objects), 

        m=length(graphics.opaque);
        if m<2,
            message.error=3;
            message.text='No previous transparency value available';
            return
        end;

        graphics.opaque=graphics.opaque(2:m);
        alpha=graphics.opaque(1);

        if graphics.mode==1,
            message.error=3;
            message.text='No transparency for line objects.';
            return
        end;


        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                set(graphics.objects(k),'FaceAlpha',alpha);
            end;
        end;
        model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics=graphics;
    end;
end;

[msg,info_text]=get_residue(indices,'info');
atom_adress=[info_text{1} '.:'];
set_object(atom_adress,'untransparency');


function message=hide_residue(indices)
% Hides a residue by deleting graphics objects
%

global model
global hMain

message.error=0;
message.text='';

graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics.objects=[];
model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics.mode=0;

if ~isempty(graphics),
    if ~isempty(graphics.objects),
        for k=1:length(graphics.objects),
            if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                delete(graphics.objects(k));
                model.unrecord=[model.unrecord graphics.objects(k)];
            end;
        end;
    end;
end;

if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'label_graphics'),
    if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics),
        gobjects=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.objects;
        if ~isempty(gobjects),
            for k=1:length(gobjects),
                if ishandle(gobjects(k)), delete(gobjects(k)); end;
            end;
        end;
        model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.objects=[];
        model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.mode=0;
    end;
end;

% [msg,info_text]=get_residue(indices,'info');
% atom_adress=[info_text{1} '.:'];
% set_object(atom_adress,'hide');

if isfield(model.structures{indices(1)}(indices(2)),'atom_graphics') && ~isempty(model.structures{indices(1)}(indices(2)).atom_graphics),
    for anum=1:length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers)
        atom=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{anum};
        [m,n]=size(atom);
        for k=1:m,
            if length(model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)})>=atom(k,1)
                graphics=model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1));
                model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1)).objects=[];
                model.structures{indices(1)}(indices(2)).atom_graphics{indices(3)}(atom(k,1)).mode=0;
                if isempty(graphics), continue; end;
                if isempty(graphics.objects), continue; end;
                for kk=1:length(graphics.objects),
                    if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                        delete(graphics.objects(kk));
                        if isfield(hMain,'unrecord'),
                            model.unrecord=[hMain.unrecord graphics.objects(kk)];
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

function message=show_residue(indices,argin)
% Plots a residue by calling plot routines for all atoms
%
% mode  string that determines the appearance of the plot

global model
mode = argin{1};

message.error=0;
message.text='';

plotted=0;

if isempty(mode),
    mode='ribbon';
end;

resname=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
water=strcmpi(resname,'HOH'); % check whether this is a water 

gobjects=[];

switch mode
    case {'CaWire','CaStick','coil','ribbon'}
        if ~water,
            gobjects=plot_backbone_segment(indices,mode);
        end;
        if ~isempty(gobjects), plotted=1; end;
    case {'wire','stick','ball&stick','space-filling'}
        if ~water,
            atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
            elnums=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).elements;  
            for k=1:length(atoms),
                gobjects=plot_atom([indices k],elnums(k),mode);
                plotted=plotted+length(gobjects);
            end;
        end;
    case {'wire!','stick!','ball&stick!','space-filling!'} % explicit water modes
        atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
        elnums=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).elements;  
        for k=1:length(atoms),
            gobjects=plot_atom([indices k],elnums(k),mode(1:end-1));
            plotted=plotted+length(gobjects);
        end;
    case 'water'
        if water,
            atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
            elnums=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).elements;  
            for k=1:length(atoms),
                gobjects=plot_atom([indices k],elnums(k),'stick');
                plotted=plotted+length(gobjects);
            end;
        end;
    case 'label'
        gobjects = plot_label(indices,0);
        plotted=plotted+length(gobjects);
    case 'mushroom'
        if length(argin) == 2
            gobjects = plot_mushroom(indices,0,argin{2});
        else
            gobjects = plot_mushroom(indices,argin{3},argin{2});
        end;
        plotted=plotted+length(gobjects);
    case 'label_frame'
        gobjects = plot_label(indices,1);        
        plotted = plotted+length(gobjects);
    case 'label_hide'
        if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics),
            gobjects=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.objects;
            if ~isempty(gobjects),
                for k=1:length(gobjects),
                    if ishandle(gobjects(k)), delete(gobjects(k)); end;
                end;
            end;
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.objects=[];
            model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).label_graphics.mode=0;
            plotted=plotted+length(gobjects);
        end;
    otherwise
        message.error=1;
        message.text='Unknown display type.';
        return
end;

if plotted==0,
    message.error=2;
    message.text='Nothing to plot.';
end;

