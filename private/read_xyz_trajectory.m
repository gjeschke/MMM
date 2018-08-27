function trajectory = read_xyz_trajectory(fname)
% function read_xyz_trajectory(fname,first,last)
%
% Reads an ASCII archive file of a molecular dynamics simulation performed
% in the form of concatenated .xyz coordinate files, as output by OpenBabel,
% or by the Perl script xtd2xyz.pl for Materials Studio
% if arguments first and last are provide, only frames first to last are read
% output is empty, if an error occurs
% set maxframes to the maximum number of trajectory frames for good memory
% management
%
% fname         file name of the trajectory file with extension
%
% trajectory    structure with fields
%               .atoms      atoms per frame
%               .elements   atomic numbers for all atoms
%               .data       data array [frames*atoms,3] of Cartesian
%                           coordinates
%
% G. Jeschke, 2008-2013

global chemistry

maxframes = 100000;

trajectory = [];

fid=fopen(fname);
if fid==-1,
    fprintf(2,'Warning! Trajectory file does not exist. Empty output.\n');
    return;
end;

trajectory.atoms = 0;
record = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if trajectory.atoms == 0, % initialization, first record
        trajectory.atoms = str2double(tline);
        tline = fgetl(fid);
        if ~ischar(tline), trajectory = []; break, end
        record = record + 1;
        trajectory.data = zeros(maxframes*trajectory.atoms,3);
        trajectory.elements = zeros(1,trajectory.atoms);
        for k = 1:trajectory.atoms,
            tline = fgetl(fid);
            if ~ischar(tline), trajectory = []; break, end
            [element, ~, ~, ni] = sscanf(tline,'%s',1);
            tline = tline(ni:end);
            atnum = str2double(element);
            if isnum(atnum),
                trajectory.elements(k) = atnum;
            else
                trajectory.elements(k) = tag2id(upper(element),upper(chemistry.element_tags));
            end;
            coor = sscanf(tline,'%f',3);
            trajectory.data(k,:) = coor;
        end;
    else
        tline = fgetl(fid);
        if ~ischar(tline), trajectory = []; break, end
        bas = record*trajectory.atoms;
        record = record + 1;
        for k = 1:trajectory.atoms,
            tline = fgetl(fid);
            if ~ischar(tline), trajectory = []; break, end
            [~, ~, ~, ni] = sscanf(tline,'%s',1);
            tline = tline(ni:end);
            coor = sscanf(tline,'%f',3);
            trajectory.data(bas+k,:) = coor;
        end;
    end;
end;

