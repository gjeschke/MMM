function update_Tinker_from_MMM(xyzname,network,snum,correspondence,options)
% function update_Tinker_from_MMM(xyzname,network,snum,correspondence,options)
%
% updates the current Tinker xyz file from internal MMM coordinates
%
% xyzname           xyz file name, the file with largest number i
%                   (extension .xyz_i) is read, as by Tinker convention
%                   this is the current structure, a file .xyz_i+1 is
%                   written
% network           C alpha network coordinates
% snum              number of corresponding MMM structure
% correspondence    correspondence table for xyz file atom numbers to MMM
%                   internal variables
%                   .network  Tinker indices of C alpha network points
%                   .chains   array of Tinker index vectors for all chains
%                             in structure snum, .chains(kc).pointers
%                   .model    model number in structure snum to which the
%                             Tinker minimization refers
% options           optional argument specifying the extent of coordinate
%                   update
%                   .full  flag indicating full coordinate update, if true
%                   if options is missing, or .full is missing, only
%                   coordinates of C alpha network points are updated
%
% G. Jeschke, 2012

global model

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

% fprintf(1,'Updating %s to %s\n',xyzname,outfile);

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
not_found=ones(size(correspondence.network));
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
            if nargin>4 && isfield(option,'full') && options.full,
                structure=model.structures{snum};
                chains=length(structure);
                for kc=1:chains,
                    snum_xyz=structure(kc).xyz{km}; % use first model
                    pointers=correspondence.chains(kc).pointers;
                    poic=find(pointers-atnum);
                    if length(poic)==1,
                        x=snum_xyz(poic,1);
                        y=snum_xyz(poic,2);
                        z=snum_xyz(poic,3);
                        break
                    end;
                end;
            else
                poic=find(correspondence.network-atnum==0);
                if length(poic)==1,
                    x=network(poic,1);
                    y=network(poic,2);
                    z=network(poic,3);
                    not_found(poic)=0;
                end;
            end;
            fprintf(fid_out,'%7i  %s',atnum,id);
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