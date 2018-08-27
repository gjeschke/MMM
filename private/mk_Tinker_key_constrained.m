function mk_Tinker_key_constrained(xyzname,network,DEER,direct,options)

direct_force_constant=100; %[kcal/(mol Å)]
DEER_force_constant=100; %[kcal/(mol Å)]
tolerance=1;

if nargin>4 && isfield(options,'direct_force_constant'),
    direct_force_constant=options.direct_force_constant;
end;
if nargin>4 && isfield(options,'DEER_force_constant'),
    DEER_force_constant=options.DEER_force_constant;
end;

xyz_files=dir([xyzname '.xyz*']);
[mypath,basname] = fileparts(xyzname);
keyfile=fullfile(mypath,sprintf('%s.key',basname));

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
 

% read Tinker .xyz file for matching coordinates and Tinker atom numbers

atnums=zeros(1,20000);
coor=zeros(20000,3);
at_ids=cell(1,2000);
poi=0;
nl=0;
mode=0;

fid=fopen(xyzname,'r');
if fid==-1,
    add_msg_board('ERROR: Tinker .xyz file does not exist');
    cd(my_path);
    return;
end;

while 1
    tline = fgetl(fid);
    if ~ischar(tline) || mode<0, break, end
    if ~isempty(tline),
        nl=nl+1;
        myline = textscan(tline,'%s');
        args=myline{1};
        if nl==1,
            num_atoms=str2double(char(args(1)));
            % add_msg_board(sprintf('%i atoms in Tinker .xyz file %s',num_atoms,fname_xyz));
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
    if poi>=num_atoms,
        mode=-1;
    end;
end;

coor=coor(1:poi,:);
atnums=atnums(1:poi);
at_ids=at_ids(1:poi);

fclose(fid);

% Atom matching and creating of constraint table

[num_direct,~]=size(direct);
num_DEER=length(DEER);

Tinker_constraints=zeros(num_direct+num_DEER,5);

% direct constraints (C_alpha)
for k=1:num_direct,
    compare=repmat(network(direct(k,1),:),poi,1);
    diff=abs(compare-coor);
    index1=find(sum(diff,2)<0.01);
    if length(index1)~=1,
        fprintf(1,'No match for first site in direct constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,1)=atnums(index1);
    end;
    compare=repmat(network(direct(k,2),:),poi,1);
    diff=abs(compare-coor);
    index2=find(sum(diff,2)<0.01);
    if length(index2)~=1,
        fprintf(1,'No match for second site in direct constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,2)=atnums(index2);
    end;
    Tinker_constraints(k,3)=direct_force_constant;
    Tinker_constraints(k,4)=10*(direct(k,4)-tolerance*direct(k,5)); % nm -> Å
    Tinker_constraints(k,5)=10*(direct(k,4)+tolerance*direct(k,5)); % nm -> Å
end;

% DEER constraints (N-O midpoint translated to Calpha)
for k=1:num_DEER,
    compare=repmat(network(DEER(k).res1,:),poi,1);
    diff=abs(compare-coor);
    index1=find(sum(diff,2)<0.01);
    if length(index1)~=1,
        fprintf(1,'No match for first site in DEER constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,1)=atnums(index1);
    end;
    compare=repmat(network(DEER(k).res2,:),poi,1);
    diff=abs(compare-coor);
    index2=find(sum(diff,2)<0.01);
    if length(index2)~=1,
        fprintf(1,'No match for second site in DEER constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,2)=atnums(index2);
    end;
    Tinker_constraints(k,3)=DEER_force_constant;
    % translate N-O midpoint distance to Calpha-Calpha distance constraint
    r_NO=norm(DEER(k).xyz2-DEER(k).xyz1);
    r_DEER=10*DEER(k).r;
    q=r_NO^2-r_DEER^2;
    r_CA=norm(network(DEER(k).res2,:)-network(DEER(k).res1,:));
    p=2*r_CA;
    dr=-p/2+sqrt(p^2/4-q);
    translated_r=r_CA+dr;
    Tinker_constraints(k,4)=translated_r-10*tolerance*DEER(k).sigr; % nm -> Å
    Tinker_constraints(k,5)=translated_r+10*tolerance*DEER(k).sigr; % nm -> Å
end;


% Write .key file

fid=fopen(keyfile,'wt');

if nargin>4 && isfield(options,'solvation') && ~isempty(options.solvation),
    fprintf(fid,'SOLVATE %s\n',options.solvation);
end;

for k=1:length(Tinker_constraints),
    fprintf(fid,'RESTRAIN-DISTANCE %i %i %5.1f %6.2f %6.2f\n',Tinker_constraints(k,:));
    % fprintf(fid,'RESTRAIN-DISTANCE %i %i\n',Tinker_constraints(k,1:2));
end;

fclose(fid);

add_msg_board(sprintf('Key file written: %s',keyfile));

