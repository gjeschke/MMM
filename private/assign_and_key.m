function correspondence = assign_and_key(network,xyzname,snum,options,km)
% function correspondence = assign_and_key(network,xyzname,snum,options,km)
%
% reads Tinker xyz file made by pdbxyz and creates a correspondence table
% into MMM structure with number snum
% optionally, a Tinker key file for minimization is written
%
% network   array of Calpha coordinates of network points
% xyzname   xyz file name
% snum      number of corresponding MMM structure
% options   key file options, all fields optional
%           .solvation  Tinker solvation keyword
%           .restrain   restraint mode, currently implemented:
%                       'all_CA'  all C alpha atoms are restrained
% km        optional number of coordinate set (model) in structure snum,
%           defaults to 1
%
% correspondence    correspondence table for xyz file atom numbers to MMM
%                   internal variables
%                   .network  Tinker indices of C alpha network points
%                   .chains   array of Tinker index vectors for all chains
%                             in structure snum, .chains(kc).pointers
%                   .model    model number in structure snum, is either km
%                             or 1 (default)
%
% G. Jeschke, 2012

global model

if nargin<5,
    km=1;
end;

correspondence=[];
xyz_files=dir([xyzname '.xyz*']);
[mypath,basname] = fileparts(xyzname);

poi=0;
maxnum=0;
for k=1:length(xyz_files);
    cnum=xyz_files(k).name(length([basname '.xyz_'])+1:end);
    if isempty(cnum) && maxnum<1,
        maxnum=1;
        poi=k;
    elseif str2double(cnum)>maxnum,
        maxnum=str2double(cnum);
        poi=k;
    end;
end;

if poi>0,
    xyzname=fullfile(mypath,xyz_files(poi).name);
end;

% read Tinker xyz coordinates

fid=fopen(xyzname,'r');
if fid==-1,
    add_msg_board('ERROR: Tinker .xyz file does not exist');
    return;
end;

poi=0;
nl=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) break, end
    if ~isempty(tline),
        nl=nl+1;
        myline = textscan(tline,'%s');
        args=myline{1};
        if nl==1,
            num_atoms=str2double(char(args(1)));
            add_msg_board(sprintf('%i atoms in Tinker .xyz file %s',num_atoms,xyzname));
        else
            atnum=str2double(char(args(1)));
            id=char(args(2));
            x=str2double(char(args(3)));
            y=str2double(char(args(4)));
            z=str2double(char(args(5)));
            poi=poi+1;
            atnums(poi)=atnum;
            coor(poi,1)=x;
            coor(poi,2)=y;
            coor(poi,3)=z;
            at_ids{poi}=id;
        end;
    end;
end;

coor=coor(1:poi,:);
atnums=atnums(1:poi);
at_ids=at_ids(1:poi);

fclose(fid);

% Atom matching and creating of constraint table

[mn,~]=size(network);
correspondence.network=zeros(mn,1);

% indexing of the network
for k=1:mn,
    compare=repmat(network(k,:),poi,1);
    diff=abs(compare-coor);
    index1=find(sum(diff,2)<0.01);
    if length(index1)~=1,
        fprintf(1,'No match for network node %i.\n',k);
        continue
    else
        correspondence.network(k)=atnums(index1);
    end;
end;

% indexing of the structure
structure=model.structures{snum};
chains=length(structure);
correspondence.model=km;

for kc=1:chains,
    snum_xyz=structure(kc).xyz{km}; % use first model
    [mn,~]=size(snum_xyz);
    correspondence.chains(kc).pointers=zeros(1,mn);
    for k=1:mn,
        compare=repmat(snum_xyz(k,:),poi,1);
        diff=abs(compare-coor);
        index1=find(sum(diff,2)<0.01);
        if length(index1)~=1,
            fprintf(1,'No match for atom %i in chain %i.\n',k,kc);
            continue
        else
            correspondence.chains(kc).pointers(k)=atnums(index1);
        end;
    end;
end;

if nargin>2 && ~isempty(options),
    keyfile=fullfile(mypath,sprintf('%s.key',basname));
    % determine indices of restrained atoms in .xyz file
    switch options.restrain
        case 'all_CA'
            restrained=correspondence.network;
        otherwise
            restrained=[];
    end;
    % write key file
    fid=fopen(keyfile,'wt');

    if isfield(options,'solvation') && ~isempty(options.solvation),
        fprintf(fid,'SOLVATE %s\n',options.solvation);
    end;
    if ~isfield(options,'restraint_force'),
        options.restraint_force=100; % this is the Tinker default
    end;
    [mr,nr]=size(restrained);
    for k=1:mr,
        switch nr
            case 1
                fprintf(fid,'RESTRAIN-POSITION %i %12.6f%12.6f%12.6f%10.2f\n',restrained(k),network(k,1),network(k,2),network(k,3),options.restraint_force);
            case 2
                fprintf(fid,'RESTRAIN-DISTANCE %i %i\n',restrained(k,:));
            case 5
                fprintf(fid,'RESTRAIN-DISTANCE %i %i %5.1f %6.2f %6.2f%10.2f\n',restrained(k,:),options.restraint_force);
        end;
    end;

    fclose(fid);

    add_msg_board(sprintf('Key file written: %s',keyfile));
end;
