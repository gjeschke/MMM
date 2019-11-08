function [message,argout]=get_residue(indices,property,varargin)
% function [message,argout]=get_residue(indices,property,varargin)
%
% Gets properties of a residue in MMM
%
% indices   index vector that identifies the atom location in the model
% property  property to be set, given as a string, available properties and
%           corresponding set functions are defined at the beginning of the
%           function source code in the switch board
% varargin  further input arguments given as a cell array varargin{:}
%
% message   error message structure with fields .error and .text,
%           .error=0 indicates no error
% argout    further output arguments given as structure argout{n}
%
% G. Jeschke, 2009-2018

argout=[];
if nargin>=3
    argin=varargin{1};
else
    argin='';
end

message.error=0;
message.text='';

switch property
    case 'absresnum'
        [message,argout]=get_absresnum(indices);
    case 'annotations'
        [message,argout]=get_annotations(indices);
    case 'Bfactor'
        [message,argout]=Bfactor_residue(indices);
    case 'children'
        [message,argout]=children_residue(indices);
    case 'Calpha'
        [message,argout]=Calpha_residue(indices);
    case 'descendants'
        [message,argout]=descendants_residue(indices);
    case 'dssp'
        [message,argout]=dssp_info(indices);
    case 'graphics'
        [message,argout]=graphics_residue(indices);
    case 'info'
        [message,argout]=info_residue(indices);
    case 'mass'
        [message,argout]=mass_residue(indices);
    case 'elements'
        [message,argout]=elements_residue(indices);
    case 'elements_heavy'
        [message,argout]=elements_heavy_residue(indices);
    case 'elements_paradigm'
        [message,argout]=elements_paradigm_residue(indices); % single copy (first alternate location of all atoms)
    case 'label'
        [message,argout] = label_residue(indices,argin);
    case 'xyz'
        [message,argout]=xyz_residue(indices);
    case 'xyz_paradigm' % single copy (first alternate location of all atoms)
        [message,argout]=xyz_paradigm_residue(indices);
    case 'xyz_heavy' % single copy of only heavy (non-hydrogen) atoms
        [message,argout]=xyz_heavy_residue(indices);
    case 'xyz_base'
        [message,argout]=xyz_base(indices);
    case 'xyz_backbone'
        [message,argout]=xyz_backbone(indices);
    case 'name'
        [message,argout]=name_residue(indices);
    case 'slc' % single-letter code
        [message,argout]=slc_residue(indices);
    otherwise
        message.error=3;
        message.text='Property does not exist';
end;

function [message,xyz]=Calpha_residue(indices)
% returns Calpha atom coordinates for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atags=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_tags;
aid=tag2id('CA',atags);
if isempty(aid),
    message.error=1;
    message.text='No Calpha atom for this residue/cofactor.';
    xyz=[];
    return;
end;
anum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers{aid};
apoi=anum(:,1)';
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
xyz=mean(coors(apoi,:),1);

function [message,newindices]=children_residue(indices)
% returns children indices for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
if isempty(atoms)
    newindices=indices;
else
    newindices=zeros(length(atoms),5);
    for k=1:length(atoms),
        newindices(k,1:4)=indices;
        newindices(k,5)=k;
    end;
end;

function [message,name]=name_residue(indices)
% returns children indices for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

name=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;

function [message,slc] = slc_residue(indices)
% returns single_letter code for a residue
%
% indices  index vector, length at least 4

global model
global residue_defs

message.error=0;
message.text='';

slc = '';
name=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;
rid = tag2id(name,residue_defs.restags);
if isempty(rid)
    message.error = 1;
    message.text = 'Unknown residue';
    return
end
slc = residue_defs.single_letter_code(rid);

function [message,absresnum]=get_absresnum(indices)
% returns the absolute residue number (for PRONOX)
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

absresnum=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).absresnum;

function [message,newindices]=descendants_residue(indices)
% returns descendants indices for a residue
%
% indices  index vector, length at least 4

message.error=0;
message.text='';

[message,mindices]=children_residue(indices);
[m,n]=size(mindices);
if n>4, % there are atom children, we have to check for children of their children
    newindices=zeros(1000,6);
    newindices(1:m,1:n)=mindices;
    poi=m;
    locflag=0;
    for k=1:m,
        [message,nindices]=get_atom(mindices(k,:),'descendants');
        [mm,nn]=size(nindices);
        if nn>5,
            locflag=1;
            newindices(poi+1:poi+mm,:)=nindices;
        else
            newindices(poi+1:poi+mm,1:5)=nindices;
        end;
        poi=poi+mm;
    end;
    newindices=newindices(1:poi,:);
    if ~locflag, % none of the atoms had alternate locations, reduce dimensionality of index array
        newindices=newindices(:,1:5);
    end;
else
    newindices=indices;
end;

function [message,info_text]=info_residue(indices)
% returns info string for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

stag=id2tag(indices(1),model.structure_tags,model.structure_ids);
ctag=id2tag(indices(2),model.chain_tags{indices(1)},model.chain_ids{indices(1)});
info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
rtag=id2tag(indices(4),model.structures{indices(1)}(indices(2)).residues{indices(3)}.residue_tags);
adr_string=mk_address(indices);
info_text={adr_string sprintf('Structure: %s',stag) sprintf('Chain      : %s',ctag) sprintf('Model      : %i',indices(3))...
    sprintf('Residue  : %s %s',info.name,rtag)};

function [message,graphics]=graphics_residue(indices)
% returns graphics handles for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

graphics=[];
if isfield(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)),'residue_graphics'),
    graphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).residue_graphics;
end;

function [message,mass]=mass_residue(indices)
% returns mass for a residue, if hydrogens are missing in the model, they
% are currently not counted!
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
na=length(atoms);
mass=0;
for k=1:na,
    locs=atoms{k};
    mass=mass+isotopes(locs(1,1),2);
end;

function [message,xyz]=xyz_residue(indices)
% returns xyz coordinates for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
na=length(atoms);
apoi=zeros(1,3*na);
poi=0;
for k=1:na,
    locs=atoms{k};
    [m,n]=size(locs);
    apoi(poi+1:poi+m)=locs(:,1)';
    poi=poi+m;
end;
apoi=apoi(1:poi);
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
xyz=coors(apoi,:);

function [message,xyz]=xyz_paradigm_residue(indices)
% returns xyz coordinates for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
na=length(atoms);
apoi=zeros(1,na);
for k=1:na,
    locs=atoms{k};
    apoi(k)=locs(1,1);
end;
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
xyz=coors(apoi,:);

function [message,xyz]=xyz_heavy_residue(indices)
% returns xyz coordinates for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
elements = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).elements;
na=length(atoms);
apoi=zeros(1,na);
for k=1:na,
    locs=atoms{k};
    apoi(k)=locs(1,1);
end;
coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
xyz=coors(apoi,:);
xyz = xyz(elements ~= 1,:);

function [message,xyz]=xyz_base(indices)
% returns xyz coordinates of heavy atoms of a nucleobase in the same
% sequence as in the PDB template connections file
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

xyz = [];

type = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).type;

if type ~= 2 % not a nucleotide
    message.error = 1;
    message.text='Residue is not a nucleotide';
    return
end

name = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).name;

switch strtrim(name)
    case {'A','DA'}
        atoms = {'N9' 'C8' 'N7' 'C5' 'C6' 'N6' 'N1' 'C2' 'N3' 'C4'};
    case {'C','DC'} 
        atoms = {'N1' 'C2' 'O2' 'N3' 'C4' 'N4' 'C5' 'C6'};
    case {'G', 'DG'}
        atoms = {'N9' 'C8' 'N7' 'C5' 'C6' 'O6' 'N1' 'C2' 'N2' 'N3' 'C4'};
    case 'U'
        atoms = {'N1' 'C2' 'O2' 'N3' 'C4' 'O4' 'C5' 'C6'};
    case 'DT'
        atoms = {'N1' 'C2' 'O2' 'N3' 'C4' 'O4' 'C5' 'C7' 'C6'};
end

na = length(atoms);
xyz = zeros(na,3);

radr = mk_address(indices);

for k = 1:na
    adr = [radr '.' atoms{k}];
    [msg,axyz] = get_object(adr,'coor');
    if msg.error
        xyz = [];
        return
    end
    xyz(k,:) = axyz;
end

function [message,xyz]=xyz_backbone(indices)
% returns xyz coordinates of heavy atoms of a nucleobase in the same
% sequence as in the PDB template connections file
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

xyz = [];

type = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).type;

if type ~= 2 % not a nucleotide
    atoms = {'N' 'CA' 'C' 'O'};
else % nucleotide
    atoms = {'P' 'O5''' 'C5''' 'C4''' 'C3'''};
end

na = length(atoms);
xyz = zeros(na,3);

radr = mk_address(indices);

for k = 1:na
    adr = [radr '.' atoms{k}];
    [msg,axyz] = get_object(adr,'coor');
    if msg.error
        xyz = [];
        return
    end
    xyz(k,:) = axyz;
end


function [message,Bfactors]=Bfactor_residue(indices)
% returns xyz coordinates for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
na=length(atoms);
apoi=zeros(1,3*na);
poi=0;
for k=1:na,
    locs=atoms{k};
    [m,n]=size(locs);
    apoi(poi+1:poi+m)=locs(:,1)';
    poi=poi+m;
end;
apoi=apoi(1:poi);
allBF=model.structures{indices(1)}(indices(2)).Bfactor{indices(3)};
Bfactors=allBF(apoi);


function [message,elements]=elements_residue(indices)
% returns atomic numbers of elements for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
na=length(atoms);
apoi=zeros(1,3*na);
poi=0;
for k=1:na,
    locs=atoms{k};
    apoi(k)=locs(1,1);
    poi=poi+m;
end;
apoi=apoi(1:poi);
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
elements=isotopes(apoi,1);

function [message,elements]=elements_heavy_residue(indices)
% returns atomic numbers of elements for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
na=length(atoms);
apoi=zeros(1,na);
for k=1:na,
    locs=atoms{k};
    apoi(k)=locs(1,1);
end;
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
elements=isotopes(apoi,1);
elements = elements(elements ~= 1);

function [message,elements]=elements_paradigm_residue(indices)
% returns atomic numbers of elements for a residue
%
% indices  index vector, length at least 4

global model

message.error=0;
message.text='';

atoms=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).atom_numbers;
na=length(atoms);
apoi=zeros(1,na);
for k=1:na,
    locs=atoms{k};
    apoi(k)=locs(1,1);
end;
isotopes=model.structures{indices(1)}(indices(2)).isotopes;
elements=isotopes(apoi,1);

function [message,dssp_info]=dssp_info(indices)
% returns the dssp information for this residue if it exists in the model,
% otherwise an empty record

global model

message.error=0;
message.txt='No error.';

dssp_info=[];

if ~isfield(model.info{indices(1)},'dssp') || isempty(model.info{indices(1)}.dssp),
    message.error=1;
    message.txt='No DSSP information for this structure.';
    return;
end;

if indices(3)~=1,
    message.error=2;
    message.txt=sprintf('DSSP information for coordinate set 1, residue in coordinate set %i.',indices(3));
end;
snum=indices(1);
dssp=model.info{indices(1)}.dssp;

for k=1:length(dssp),
    cnum=tag2id(dssp(k).chain,model.chain_tags{snum});
    if cnum==indices(2),
        rnum=tag2id(dssp(k).tag,model.structures{snum}(cnum).residues{1}.residue_tags);
        if rnum==indices(4),
            dssp_info=dssp(k);
            dssp_info.my_num=k;
            break;
        end;
    end;
end;

function [msg,labels] = label_residue(indices,label_name)

global model
global label_defs

poi=0;
for k0=1:length(model.sites),
    for k1=1:length(model.sites{k0}),
        for k=1:length(model.sites{k0}(k1).residue),
            if abs(sum(model.sites{k0}(k1).residue(k).indices - indices)) == 0, % label of the current residue
                name1 = model.sites{k0}(k1).residue(k).label;
                id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                name2 = label_defs.residues(id).short_name;
                if isempty(label_name) || strcmpi(label_name,name1) || strcmpi(label_name,name2),
                    poi=poi+1;
                    id=tag2id(model.sites{k0}(k1).residue(k).label,label_defs.restags);
                    labels(poi).name=label_defs.residues(id).short_name;
                    labels(poi).T=model.sites{k0}(k1).residue(k).T;
                    labels(poi).NOpos=model.sites{k0}(k1).residue(k).NOpos;
                    x=sum(labels(poi).NOpos(:,1).*labels(poi).NOpos(:,4));
                    y=sum(labels(poi).NOpos(:,2).*labels(poi).NOpos(:,4));
                    z=sum(labels(poi).NOpos(:,3).*labels(poi).NOpos(:,4));
                    labels(poi).xyz=[x y z];
                    labels(poi).rmsd=NOpos_rmsd(labels(poi).NOpos);
                end;
            end;
        end;
    end;
end;

if ~exist('labels','var'),
    labels = [];
    msg.error = 1;
    msg.txt = 'Residue is not labelled.';
else
    msg.error = 0;
    msg.txt = 'No error.';
end;

function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm
