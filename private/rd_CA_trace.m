function [coor,N,n] = rd_CA_trace(fname)
% [coor,N,n] = rd_CA_trace(fname)
%
% Reads Calpha coordinates from a PDB file into an [N*n,3] array coor,
% where N is the number of models and n the number of CA atoms per model
%
% fname     file name, if no extension is given ' .pdb'  is appended
%
% if not all models have the same number of CA atoms, a warning is output
% and output parameter n is empty
%
% G. Jeschke, 24.12.2019

N = 0;
n = [];
coor = zeros(1000000,3);
nvec = zeros(1,10000);

if ~contains(fname,'.')
    fname = strcat(fname,'.pdb');
end

ifid=fopen(fname);
if ifid==-1
    warning('rd_CA_trace:Input file could not be opened');
    return;
end

cpoi = 0;
nc = 0;
while 1
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if length(tline) >= 54 
        if strcmpi(tline(1:4),'ATOM') && strcmpi(tline(14:15),'CA')
            cpoi = cpoi + 1;
            nc = nc + 1;
            coor(cpoi,1) = str2double(tline(31:38)); 
            coor(cpoi,2) = str2double(tline(39:46)); 
            coor(cpoi,3) = str2double(tline(47:54)); 
        end
    end
    if length(tline) >= 6 && strcmpi(tline(1:6),'ENDMDL')
        N = N + 1;
        nvec(N) = nc;
        nc = 0;
    end
end
fclose(ifid);

coor = coor(1:cpoi,:);

if N == 0
    nvec(1) = nc;
    N = 1;
end

nvec = nvec(1:N);
if max(abs(nvec - mean(nvec))) ~= 0
    warning('rd_CA_trace: MIsmatch in the number of CA atoms per model');
    return
end

n = nvec(1);