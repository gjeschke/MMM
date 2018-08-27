function network=update_from_Tinker(xyzname,snum,correspondence,xyznum)
% function update_from_Tinker(xyzname,snum,correspondence)
%
% reads Tinker xyz file made by minimize and updates corresponding
% coordinates in MMM structure snum, if optional argument network is given,
% the new network coordinates are also assigned, if correspondence.network
% exists
%
% xyzname           xyz file name, the file with largest number i
%                   (extension .xyz_i) is read, as by Tinker convention
%                   this is the current structure  (unless argument xyznum
%                   is provided)
% snum              number of corresponding MMM structure
% correspondence    correspondence table for xyz file atom numbers to MMM
%                   internal variables
%                   .network  Tinker indices of C alpha network points
%                   .chains   array of Tinker index vectors for all chains
%                             in structure snum, .chains(kc).pointers
%                   .model    model number in structure snum to which the
%                             Tinker minimization refers
% xyznum            optional argument specifying the number of the xyz file
%                   corresponding to file extension '.xyz_%i',xyznum, if
%                   empty, the file with largest number is used
%
% network           updated network coordinates
%
% G. Jeschke, 2012

global model

network=[];

xyz_files=dir([xyzname '.xyz*']);
[mypath,basname] = fileparts(xyzname);

if nargin<4 || isempty(xyznum),
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
else
    poi=xyznum;
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

% Update network coordinates, if requested

if nargout>0 && isfield(correspondence,'network') && ~isempty(correspondence.network),
    mn=length(correspondence.network);
    network=zeros(mn,3);
    for k=1:mn,
        network(k,:)=coor(correspondence.network(k),:);
    end;
end;

if isfield(correspondence,'chains') && isfield(correspondence,'model'),
    % update of the structure
    structure=model.structures{snum};
    chains=length(structure);
    if chains~=length(correspondence.chains),
        add_msg_board('ERROR: Mismatch in number of chains. No coordinate update.');
        return;
    end;
    km=correspondence.model;
    for kc=1:chains,
        snum_xyz=structure(kc).xyz{km}; % use first model
        [mn,~]=size(snum_xyz);
        for k=1:mn,
            atnum=correspondence.chains(kc).pointers(k);
            if atnum~=0,
                snum_xyz(k,:)=coor(atnum,:);
            else
                fprintf(1,'Correspondence information missing for atom %i in chain %i\n',k,kc);
            end;
        end;
        structure(kc).xyz{km}=snum_xyz;
    end;
end;

