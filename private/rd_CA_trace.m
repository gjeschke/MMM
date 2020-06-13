function [coor,N,n,resnum,sequence] = rd_CA_trace(fname,chain_ID)
% [coor,N,n] = rd_CA_trace(fname)
%
% Reads Calpha coordinates from a PDB file into an [N*n,3] array coor,
% where N is the number of models and n the number of CA atoms per model
% 
% fname     file name, if no extension is given ' .pdb'  is appended
% chain_ID  optional chain selector, if given, only the chain with this
%           identifier is read
%
% resnum    residue number axis
% sequence  sequence in single-letter code
%
% if not all models have the same number of CA atoms, a warning is output
% and output parameter n is empty
%
% G. Jeschke, 24.12.2019

global residue_defs

if ~exist('chain_ID','var') || isempty(chain_ID)
    chain_ID = '*';
end

N = 0;
n = [];
coor = zeros(1000000,3);
nvec = zeros(1,10000);
resnum = zeros(1,10000);
sequence = char(32*ones(1,10000));

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
        cid = tline(22);
        if chain_ID == '*' || strcmpi(cid,chain_ID)
            read_chain = true;
        else
            read_chain = false;
        end
        if strcmpi(tline(1:4),'ATOM') && strcmpi(tline(14:15),'CA') && read_chain
            cpoi = cpoi + 1;
            nc = nc + 1;
            coor(cpoi,1) = str2double(tline(31:38)); 
            coor(cpoi,2) = str2double(tline(39:46)); 
            coor(cpoi,3) = str2double(tline(47:54));
            if N == 0
                resnum(cpoi) = str2double(tline(23:26));
                amino_id=tag2id(tline(18:20),upper(residue_defs.restags),residue_defs.single_letter_code); % test whether this is an amino acid
                if isempty(amino_id)
                    amino_id = '?';
                end
                sequence(cpoi) = amino_id;
            end
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
resnum = resnum(1:n);
sequence = sequence(1:n);
