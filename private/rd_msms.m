function [tri,x,y,z,info]=rd_msms(filename)
% function [tri,x,y,z]=rd_msms(filename)
% 
% Reads the triangle definitions tri and vertices (x,y,z) that define the
% solvent excluded surface as computed with the MSMS program 
% by Michel F. Sanner 
% (see: http://www.scripps.edu/~sanner/html/msms_home.html)
%
% the proper citation is:
% Sanner, M. F., Spehner, J.-C., Olson, A. J. (1996) Reduced surface: an
% efficient way to compute molecular surfaces. Biopolymers, 38(3), 305-320.
%
% the program MSMS is not (and must not be) bundled with MMM, users have to
% obtain it for themselves or for their institution
%
% files filname.face and filename.vert must exist, otherwise the
% corresponding variables are empty, info provides parameters
%
% G. Jeschke, 2009

tri=[];
x=[];
y=[];
z=[];

info.header='';
info.faces=[];
info.spheres=[];
info.density=[];
info.probe_r=[];

fid=fopen([filename '.face'],'r');

if fid~=-1,
    tline = fgetl(fid);
    info.header=strtrim(tline(2:end));
    tline = fgetl(fid);
    info.faces = fscanf(fid,'%i',1);
    info.spheres = fscanf(fid,'%i',1);
    info.density = fscanf(fid,'%f',1);
    info.probe_r = fscanf(fid,'%f',1);
    tline = fgetl(fid);
    tri=zeros(info.faces,3);
    for k=1:info.faces,
        tline = fgetl(fid);
        tri(k,:)=sscanf(tline,'%i',3);
    end;
end;

fclose(fid);

fid=fopen([filename '.vert'],'r');

if fid~=-1,
    tline = fgetl(fid);
    tline = fgetl(fid);
    info.vertices = fscanf(fid,'%i',1);
    coor=zeros(info.vertices,3);
    tline = fgetl(fid);
    for k=1:info.vertices,
        tline = fgetl(fid);
        coor(k,:)=sscanf(tline,'%f',3);
    end;
end;

fclose(fid);

x=coor(:,1)';
y=coor(:,2)';
z=coor(:,3)';
