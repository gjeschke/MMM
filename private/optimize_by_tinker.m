function [msg,snum,energy] = optimize_by_tinker(basname,molecule,selected,options)
% [msg,snum,energy] = optimize_by_tinker(basname,molecule,selected,options)
%
% Optimizes geometry of a selected part of a structure by calling Tinker
%
% basname   name of the input molecular system file, if empty, a default
%           name tinker_mol is used
% molecule  cell array of index arrays of objects that make up the whole 
%           input molecule, defaults to current structure
% selected  cell array of index arrays of objects that are active or 
%           inactive during the Tinker computation (depending on options), 
%           defaults to the current selection
% options   Tinker computation options, contains also restraints, see
%           wr_tinker_key source code for possibilities and defaults
%
% msg       message on success, msg.error is error code (0 for success) and
%           msg.text the error message
% snum      structure index of the optimized structure
% energy    energy in J/mol of the optimized structure, NaN if the
%           optimization was not successful
%
% G. Jeschke, 12.7.-31.8.2017

global model

if isfield(options,'initial') && options.initial
    initial_structure = true; 
else
    initial_structure = false; 
end
force_newton = false; % only for debugging, set fals for production

kcal2J = 4186.8;

snum = [];
energy = NaN;

if ~isfield(options,'tolerance') || isempty(options.tolerance)
    options.tolerance = 0.1;
end

set(gcf,'Pointer','watch');

[PDBind,Tink2PDB,snum0,msg] = prepare_tinker_input(basname,molecule,selected,options);

if msg.error > 0
    add_msg_board(sprintf('ERROR (prepare_tinker_input): %s',msg.text));
    set(gcf,'Pointer','arrow');
    return
elseif msg.error < 0
    add_msg_board(sprintf('Warning (prepare_tinker_input): %s',msg.text));
end
    
tinker_path = get_tinker_path;

if isempty(tinker_path)
    msg.error = 1;
    msg.text = 'Tinker distribution (file amber99.prm) not found on Matlab path.';
    set(gcf,'Pointer','arrow');
    return
end

optimizer = 'newton';
if ~force_newton && options.tolerance >= 0.1
    optimizer = 'minimize';
end

cmd=sprintf('%s %s ..%sparams%s%s.prm %6.3f',optimizer,fullfile(tinker_path,basname),filesep,filesep,options.forcefield,options.tolerance);

if strcmpi(optimizer,'newton')
    cmd = sprintf('%s A 0.01',cmd);
end

dospath=which(sprintf('%s.exe',optimizer));
if isempty(dospath)
    msg .error = 1;
    msg.text = sprintf('Tinker module %s not found on Matlab path.',optimizer);
    set(gcf,'Pointer','arrow');
    return
end
[modpath, ~] = fileparts(dospath);
my_dir=pwd;
if ~initial_structure
    cd(modpath);
    add_msg_board(sprintf('Now calling Tinker %s with tolerance %6.3f kcal/mol',optimizer,options.tolerance));
    tic;
    [s, w] = dos(cmd);
    runtime=toc;
    add_msg_board(sprintf('Tinker module %s was running %5.1f s\n',optimizer,runtime));
    if s~=0
        rem=w;
        while ~isempty(rem)
            [token,rem]=strtok(rem,char(10));
            if ~isempty(token)
                add_msg_board(token);
            end
        end
        msg.error = 2;
        msg.text = sprintf('ERROR: Tinker module %s did not run successfully.',optimizer);
        cd(my_dir);
        set(gcf,'Pointer','arrow');
        return
    else % get the energy
        poi = strfind(w,'Final Function Value :');
        if ~isempty(poi)
            args = textscan(w(poi+length('Final Function Value :'):end),'%s');
            energy = kcal2J*str2double(args{1}(1));
        end
    end
    add_msg_board(sprintf('Tinker energy is %8.2f kJ/mol\n',energy/1000));

    cd(my_dir);
    if isnan(energy)
        return
    end
end

invalid_tag = true;
ipoi = 0;
while invalid_tag && ipoi < 999
    ipoi = ipoi + 1;
    inum = sprintf('%i',ipoi);
    while length(inum)<3
        inum= strcat('0',inum);
    end
    tag = strcat('+',inum);
    id = tag2id(tag,model.structure_tags);
    if isempty(id)
        invalid_tag = false;
    end
end

if invalid_tag
    add_msg_board('ERROR: Have run out of structure tags.');
    set(gcf,'Pointer','arrow');
    cd(my_dir);
    return
end
    

[snum,message]=copy_structure(snum0,tag);

if message.error
    add_msg_board(sprintf('ERROR (copy_structure): %s',message.text));
    set(gcf,'Pointer','arrow');
    cd(my_dir);
    return
end

if ~initial_structure
    % write updated coordinates to copied structure
    ecoor = read_tinker_xyz(basname,2);
    [nt,~] = size(ecoor);
    for k = 1:nt
        if Tink2PDB(k) > 0
            lind = PDBind(Tink2PDB(k),:);
            lind(1) = snum;
            ncoor{1} = ecoor(k,2:4);
            set_location(lind,'xyz',ncoor);
        end
    end
end

% dospath=which('xyzpdb.exe');
% if isempty(dospath)
%     add_msg_board('Tinker module xyzpdb not found on Matlab path.');
%     set(gcf,'Pointer','arrow');
%     return
% end
% [modpath, ~] = fileparts(dospath);
% add_msg_board('Now calling Tinker xyzpdb');
% 
% cmd = sprintf('xyzpdb ..%s%s ..%sparams%s%s.prm',filesep,basname,filesep,filesep,options.forcefield);
% 
% my_dir = pwd;
% cd(modpath);
% tic,
% [s, w] = dos(cmd);
% runtime=toc;
% add_msg_board(sprintf('Tinker xyzpdb was running %5.1f s\n',runtime));
% if s~=0
%     rem=w;
%     while ~isempty(rem)
%         [token,rem]=strtok(rem,char(10));
%         if ~isempty(token)
%             add_msg_board(token);
%         end
%     end
%     add_msg_board('ERROR: Tinker xyzpdb did not run successfully.');
%     cd(my_dir);
%     set(gcf,'Pointer','arrow');
%     return
% end

cd(my_dir);
set(gcf,'Pointer','arrow');
