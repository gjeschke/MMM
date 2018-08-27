function [tri,x,y,z,info,tnorm]=get_SAS(coor,elements,radius,MSMS_path)
% function get_SAS(coor,elements)
%
% returns the solvent-accessible surface using MSMS
% if the routine fails, empty outputs are returned
%
% Input:
% coor      Cartesian coordinates of atoms
% elements  corresponding atom numbers in the PSE
% radius    optional, radius of the probe sphere, defaults to 1.5 Å
% MSMS_path optional, path to MSMS (for speed), determined if not given
% 
%
% Output:
% tri       triangle vertex numbers in coordinate vectors
% x         x coordinates of vertices
% y         y coordinates of vertices
% z         z coordinates of vertices
% info      additional information on SAS computation
%
% G. Jeschke, 2010

global chemistry
global general

tri=[];
x=[];
y=[];
z=[];
info=[];

if nargin<3,
    radius=1.5;
end;

if nargin<4,
    MSMS_path=which('msms.exe');
    if isempty(MSMS_path),
        return
    end;
end;

vdW=zeros(size(elements));
for k=1:length(elements),
    if elements(k)>0 && elements(k)<length(chemistry.pse),
        vdW(k)=chemistry.pse(elements(k)).vdW;
    end;
end;

outfile=[general.tmp_files 'tmp.xyzr'];
fid=fopen(outfile,'w');
if fid==-1,
    add_msg_board('ERROR: Coordinate file could not be opened for writing.');
    return
end;
for k=1:length(vdW),
    if vdW(k)>0,
        fprintf(fid,'%12.3f%12.3f%12.3f%6.2f',coor(k,1),coor(k,2),coor(k,3),vdW(k));
        if k<length(vdW),
            fprintf(fid,'\n');
        end;
    end;
end;
fclose(fid);

dstr=' -density 1';
radstr=sprintf(' -probe_radius %4.1f',radius);

[pathstr, name, ext, versn] = fileparts(outfile);
outfile=fullfile(pathstr,[name ext versn]);

msmsfile=[general.tmp_files 'tmp'];
[pathstr, name, ext, versn] = fileparts(msmsfile);
msmsfile=fullfile(pathstr,[name ext versn]);
cmd=[MSMS_path ' -if ' outfile ' -of ' msmsfile dstr radstr];
[s, w] = dos(cmd);
if s~=0,
    return
end;

[tri,x,y,z,info]=rd_msms(msmsfile);

vertices=[x',y',z'];

e1=vertices(tri(:,2),:)-vertices(tri(:,1),:);
e2=vertices(tri(:,3),:)-vertices(tri(:,1),:);

tnorm=cross(e1,e2); % vertex normals
[mt,nt]=size(tnorm);
% normalize vertex normals and catch degenerate triangles
for k=1:mt,
    nv=norm(tnorm(k,:));
    if nv<eps,
        tnorm(k,:)=[0,0,0];
    else
        tnorm(k,:)=tnorm(k,:)/nv;
    end;
end;
