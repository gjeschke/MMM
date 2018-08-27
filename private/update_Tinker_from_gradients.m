function update_Tinker_from_gradients(xyzname,gradients,smallstep)
% function update_Tinker_from_gradients(xyzname,gradients,smallstep)
%
% updates the current Tinker xyz file from gradients
%
% xyzname           xyz file name, the file with largest number i
%                   (extension .xyz_i) is read, as by Tinker convention
%                   this is the current structure, a file .xyz_i+1 is
%                   written
% gradients         coordinate gradients as obtained with testgrad
% smallstep         maximum coordinate step, defaults to 0.05 Å
%
% G. Jeschke, 2012

if nargin<3,
    smallstep=0.05;
end;

steps=sqrt(sum(gradients.^2,2));
sc=smallstep/max(steps);
dxyz=-sc*gradients;

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

outfile=fullfile(mypath,sprintf('%s.xyz_%i',basname,maxnum+1));

% read and rewrite Tinker xyz coordinates

fid=fopen(xyzname,'r');
if fid==-1,
    add_msg_board('ERROR: Tinker .xyz file does not exist');
    return;
end;
fid_out=fopen(outfile,'wt');
if fid_out==-1,
    add_msg_board('ERROR: Tinker .xyz file cannot be written');
    return;
end;
poi=0;
nl=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break, 
    end
    if ~isempty(tline),
        nl=nl+1;
        myline = textscan(tline,'%s');
        args=myline{1};
        if nl==1,
            num_atoms=str2double(char(args(1)));
            add_msg_board(sprintf('%i atoms in Tinker .xyz file %s',num_atoms,xyzname));
            fprintf(fid_out,'%s\n',tline);
        else
            atnum=str2double(char(args(1)));
            id=char(args(2));
            x=str2double(char(args(3)));
            y=str2double(char(args(4)));
            z=str2double(char(args(5)));
            poi=poi+1;
            x=x+dxyz(atnum,1);
            y=y+dxyz(atnum,2);
            z=z+dxyz(atnum,3);
            fprintf(fid_out,'%6i  %s',atnum,id);
            for ks=1:3-length(id),
                fprintf(fid_out,' ');
            end;
            fprintf(fid_out,'%12.6f%12.6f%12.6f',x,y,z);
            fprintf(fid_out,'%s\n',tline(48:end));
        end;
    end;
end;

fclose(fid);
fclose(fid_out);

% fprintf(1,'Out of %i network nodes, %i were not found in Tinker file.\n',length(correspondence.network),sum(not_found));