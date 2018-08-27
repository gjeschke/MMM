function [coor,residues,frames]=rd_pdb_trajectory_coarse(fname)
% function [coor,residues,frames]=rd_pdb_trajectory_coarse(fname)
%
% Reads Calpha coordinates from a PDB trajectory file into a Matlab 
% coordinate array
% each frame should be defined as a MODEL in the PDB file
% fname         filename of the trajectory file
%
% coor          coordinate array [F*R,3] of Calpah atom Cartesian
%               coordinates for F frames with R residues each, frame by
%               frame
% residues      number R of residues per frame
% frames        number F of frames
%
% G. Jeschke, 2013


coor=[];
residues=0;
frames=0;
cres=0;

if isempty(fname), 
    return
end;

fid=fopen(fname);
if fid==-1,
    return;
end;

coor=zeros(5000000,3);

maxpoi=0;

nl=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    while length(tline)<6, % catches too short end line of NMsim files
        tline=[tline ' '];
    end;
    record=tline(1:6);
    switch record
        case 'MODEL '
            frames=frames+1;
            if length(tline)>=11,
                curr_model=sscanf(tline(11:end),'%i');
            else
                curr_model=frames;
            end;                
        case 'ENDMDL'
            residues=cres;
            cres=0;
        case {'ATOM  '}
            atom_tag=tline(13:16);
            atom_tag=strtrim(atom_tag);
            if strcmpi(atom_tag,'CA'),
                cres=cres+1;
                % store the atom coordinates and B-factor in the appropriate array
                x=sscanf(tline(31:38),'%f');
                y=sscanf(tline(39:46),'%f');
                z=sscanf(tline(47:54),'%f');
                poi=residues*(curr_model-1)+cres;
                coor(poi,:)=[x,y,z];
                if poi>maxpoi,
                    maxpoi=poi;
                end;
            end;
        case 'END   '
            break;
    end;
end
coor=coor(1:maxpoi,:);

fclose(fid);
