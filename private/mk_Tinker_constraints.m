function mk_Tinker_constraints

% prototype version
% G. Jeschke, 2012

global general

solvation=true;

diagnostics=false; % set true for test purposes

tolerance=1;
force_constant=100; % Tinker force constant in kcal/mol/Å for Calpha constraints
DEER_force_constant=100; % Tinker force constant in kcal/mol/Å for DEER constraints

add_msg_board('### Generating Tinker key file with constraints ###');

num_direct=0;
num_DEER=0;

% Get constraint file name and read constraint file
my_path=pwd;
cd(general.restraint_files);

[fname,pname]=uigetfile('*.dat','Load constraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Constraint loading cancelled by user');
    cd(my_path);
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    constraints=rd_restraints(fullfile(pname,fname));
end;

cd(general.pdb_files);
[FileName,PathName] = uigetfile({'*.pdb;*.pdb1;*.ent'},'Select PDB file');
if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Loading of PDB file canceled by user.');
    message.error=5;
    message.text='Cancelled.';
    cd(my_path);
    return;
end;
reset_user_paths(PathName);
general.pdb_files=PathName;
initial=fullfile(PathName,FileName);
initial_PDB=FileName(1:4);

[message,initial_num]=add_pdb(initial);

if message.error,
    add_msg_board(message.text);
    cd(my_path);
    return
end;

[pathstr, name] = fileparts(initial);

fname_xyz=fullfile(pathstr,sprintf('%s.xyz',name));
fname_key=fullfile(pathstr,sprintf('%s.key',name));

% get coordinates of all constraint sites

if isfield(constraints,'direct'),
    num_direct=length(constraints.direct);
    for k=1:length(constraints.direct),
        cindices1=resolve_address(sprintf('%s.CA',constraints.direct(k).adr1));
        if isempty(cindices1) || length(cindices1)~=5,
            add_msg_board(sprintf('Warning: Direct constraint site %s does not exist in template. Constraint %i will be ignored.',constraints.direct(k).adr1,k));
            continue
        else
           [message,xyz]=get_atom(cindices1,'coor');
           if message.error,
                add_msg_board(sprintf('Warning: Failed to find Calpha coordinates for direct constraint site %s (Errror: %s). Constraint %i will be ignored.',constraints.direct(k).adr1,message.text,k));
           else
               constraints.direct(k).xyz1=xyz;
           end;
        end;
        cindices2=resolve_address(sprintf('%s.CA',constraints.direct(k).adr2));
        if isempty(cindices1) || length(cindices2)~=5,
            add_msg_board(sprintf('Warning: Direct constraint site %s does not exist in template. Constraint %i will be ignored.',constraints.direct(k).adr2,k));
            continue
        else
           [message,xyz]=get_atom(cindices2,'coor');
           if message.error,
                add_msg_board(sprintf('Warning: Failed to find Calpha coordinates for direct constraint site %s (Errror: %s). Constraint %i will be ignored.',constraints.direct(k).adr2,message.text,k));
           else
               constraints.direct(k).xyz2=xyz;
           end;
        end;
    end;
elseif isfield(constraints,'DEER') && ~isempty(constraints.DEER),
    add_msg_board(sprintf('Processing %i DEER constraints.',length(constraints.DEER)));
    [DEER,cancelled]=process_DEER_restraints(constraints);
    if cancelled,
        add_msg_board('ERROR: DEER constraint processing cancelled by user. Aborting.')
        cd(my_path);
        return
    end;
    num_DEER=length(DEER);
    for k=1:length(DEER),
        cindices1=resolve_address(sprintf('%s.CA',DEER(k).adr1));
        if isempty(cindices1) || length(cindices1)~=5,
            add_msg_board(sprintf('Warning: DEER constraint site %s does not exist in template. Constraint %i will be ignored.',constraints.direct(k).adr1,k));
            continue
        else
           [message,xyz]=get_atom(cindices1,'coor');
           if message.error,
                add_msg_board(sprintf('Warning: Failed to find Calpha coordinates for DEER constraint site %s (Errror: %s). Constraint %i will be ignored.',constraints.direct(k).adr1,message.text,k));
           else
               DEER(k).Ca_xyz1=xyz;
           end;
        end;
        cindices2=resolve_address(sprintf('%s.CA',DEER(k).adr2));
        if isempty(cindices1) || length(cindices2)~=5,
            add_msg_board(sprintf('Warning: DEER constraint site %s does not exist in template. Constraint %i will be ignored.',constraints.direct(k).adr2,k));
            continue
        else
           [message,xyz]=get_atom(cindices2,'coor');
           if message.error,
                add_msg_board(sprintf('Warning: Failed to find Calpha coordinates for DEER constraint site %s (Errror: %s). Constraint %i will be ignored.',constraints.direct(k).adr2,message.text,k));
           else
               DEER(k).Ca_xyz2=xyz;
           end;
        end;
    end;    
else
    add_msg_board('Warning: No valid constraints in constraint file. Aborting.')
    cd(my_path);
    return
end;

% read Tinker .xyz file for matching coordinates and Tinker atom numbers

atnums=zeros(1,20000);
coor=zeros(20000,3);
at_ids=cell(1,2000);
poi=0;
nl=0;
mode=0;

fid=fopen(fname_xyz,'r');
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
            add_msg_board(sprintf('%i atoms in Tinker .xyz file %s',num_atoms,fname_xyz));
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

Tinker_constraints=zeros(num_direct+num_DEER,5);

% direct constraints (C_alpha)
for k=1:num_direct,
    compare=repmat(constraints.direct(k).xyz1,poi,1);
    diff=abs(compare-coor);
    index1=find(sum(diff,2)<0.01);
    if length(index1)~=1,
        fprintf(1,'No match for first site in direct constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,1)=atnums(index1);
    end;
    compare=repmat(constraints.direct(k).xyz2,poi,1);
    diff=abs(compare-coor);
    index2=find(sum(diff,2)<0.01);
    if length(index2)~=1,
        fprintf(1,'No match for second site in direct constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,2)=atnums(index2);
    end;
    Tinker_constraints(k,3)=force_constant;
    Tinker_constraints(k,4)=10*(constraints.direct(k).r-tolerance*constraints.direct(k).sigr); % nm -> Å
    Tinker_constraints(k,5)=10*(constraints.direct(k).r+tolerance*constraints.direct(k).sigr); % nm -> Å
    if diagnostics,
        fprintf(1,'Atom identifiers of constraint %i sites are: %s and %s.\n',k,char(at_ids{index1}),char(at_ids{index2}));
    end;
end;

% DEER constraints (N-O midpoint translated to Calpha)
for k=1:num_DEER,
    compare=repmat(DEER(k).Ca_xyz1,poi,1);
    diff=abs(compare-coor);
    index1=find(sum(diff,2)<0.01);
    if length(index1)~=1,
        fprintf(1,'No match for first site in DEER constraint %i.\n',k);
        continue
    else
        Tinker_constraints(k,1)=atnums(index1);
    end;
    compare=repmat(DEER(k).Ca_xyz2,poi,1);
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
    r_CA=norm(DEER(k).Ca_xyz2-DEER(k).Ca_xyz1);
    p=2*r_CA;
    dr=-p/2+sqrt(p^2/4-q);
    translated_r=r_CA+dr;
    Tinker_constraints(k,4)=translated_r-10*tolerance*DEER(k).sigr; % nm -> Å
    Tinker_constraints(k,5)=translated_r+10*tolerance*DEER(k).sigr; % nm -> Å
    if diagnostics,
        fprintf(1,'Atom identifiers of constraint %i sites are: %s and %s.\n',k,char(at_ids{index1}),char(at_ids{index2}));
    end;
end;


% Write .key file

fid=fopen(fname_key,'wt');

if solvation,
    fprintf(fid,'SOLVATE STILL\n');
end;
for k=1:length(Tinker_constraints),
    fprintf(fid,'RESTRAIN-DISTANCE %i %i %5.1f %6.2f %6.2f\n',Tinker_constraints(k,:));
    % fprintf(fid,'RESTRAIN-DISTANCE %i %i\n',Tinker_constraints(k,1:2));
end;

fclose(fid);

add_msg_board(sprintf('Key file written: %s',fname_key));

cd(my_path);
drawnow;