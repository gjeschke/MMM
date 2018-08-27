function [energy,components] = get_Tinker_energy(full_bas,full_force_field)
% function [energy,components] = get_Tinker_energy(full_bas,full_force_field)
%
% get Tinker energy and its components

energy=0;
components.unrestrained=0;

cmd=sprintf('analyze %s %s E',full_bas,full_force_field);

dospath=which('analyze.exe');
if isempty(dospath),
    message.error=2;
    message.text='Tinker software not found on Matlab path.';
    add_msg_board('This feature requires Tinker''s analyze from the Ponder lab');
    add_msg_board('ERROR: Tinker could not be found on the Matlab path.');
    add_msg_board('Please check whether Tinker is installed and the path set.');
    % add_msg_board('(see also help browser)');
    % webcall(entry,'-helpbrowser');
    return
end;
[modpath, modcmd] = fileparts(dospath);
my_dir=pwd;
cd(modpath);
tic;
[s, w] = dos(cmd);
runtime=toc;
% add_msg_board(sprintf('Tinker analyze was running %5.1f s',runtime));
if s~=0,
    rem=w;
    while ~isempty(rem),
        [token,rem]=strtok(rem,char(10));
        if ~isempty(token),
            add_msg_board(token);
        end;
    end;
    message.error=2;
    message.text='Tinker error.';
    add_msg_board('ERROR: Tinker did not run successfully.');
    % set(hwin,'Pointer','arrow');
    cd(my_dir);
    return
end;

comments=textscan(w,'%s','Delimiter','\n');
lines=comments{1};
for k=1:length(lines),
    tline=char(lines(k));
    if ~isempty(tline),
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmpi(char(args(1)),'Total') && strcmpi(char(args(2)),'Potential') && strcmpi(char(args(3)),'Energy'),
            energy=str2double(char(args(5)));
        end;
        if strcmpi(char(args(1)),'Bond') && strcmpi(char(args(2)),'Stretching'),
            components.bond_stretching=str2double(char(args(3)));
        end;
        if strcmpi(char(args(1)),'Angle') && strcmpi(char(args(2)),'Bending'),
            components.angle_bending=str2double(char(args(3)));
        end;
        if strcmpi(char(args(1)),'Improper') && strcmpi(char(args(2)),'Torsion'),
            components.improper_torsion=str2double(char(args(3)));
        end;
        if strcmpi(char(args(1)),'Van') && strcmpi(char(args(2)),'der') && strcmpi(char(args(3)),'Waals'),
            components.van_der_Waals=str2double(char(args(4)));
        end;
        if strcmpi(char(args(1)),'Charge-Charge'),
            components.charge_charge=str2double(char(args(2)));
        end;
        if strcmpi(char(args(1)),'Implicit') && strcmpi(char(args(2)),'Solvation'),
            components.implicit_solvation=str2double(char(args(3)));
        end;
        if strcmpi(char(args(1)),'Geometric') && strcmpi(char(args(2)),'Restraints'),
            components.geometric_restraints=str2double(char(args(3)));
        end;
    end;
end;

if isfield(components,'geometric_restraints'),
    components.unrestrained=energy-components.geometric_restraints;
else
    components.unrestrained=energy;
end;
cd(my_dir);