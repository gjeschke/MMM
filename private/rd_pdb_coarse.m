function [protein,chain_tags,chain_ids]=rd_pdb_coarse(pdbid)
% function protein=rd_pdb_coarse(fname)
%
% Reads partial structure information from a PDB file into a Matlab structure
% only the 
%   residue numbers, residue types, C_alpha and C_beta coordinates and 
%   database references of all chains are stored
% pdbid         PDB identifier
% protein       protein structure information, array of structures for all
%               chains (index k)
%               protein(k).dbref    data base reference for sequence
%               protein(k).resnum   residue numbers
%               protein(k).restype  residue types, single letter code
%               protein(k).Ca       Calpha coordinates
%               protein(k).Cb       Cbeta coordinates
% chain_tags    list of chain tags in MMM format
% chain_ids     list of chain identifiers
%
% G. Jeschke, 2012


global residue_defs

fname = get_pdb_file(pdbid);

if isempty(fname), 
    protein=[];
    chain_tags=':';
    chain_ids=[];
    return
end;

db_tags=':';
db_access={};
db_cutout=zeros(50,2);
chain_tags=':';
chain_id=0;
chain_ids=[];

pointers=zeros(100,2);


fid=fopen(fname);
if fid==-1,
    return;
end;

nl=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    if length(tline)>=6, % catches too short end line of MolProbity files
        record=tline(1:6);
    else
        record='END   ';
    end;
    switch record
        case 'DBREF '
            if length(tline)>=67,
                tag=tline(13);
                db_id=tag2id(tag,db_tags);                
                if isempty(db_id),
                    db_id=length(db_access)+1;
                    db_tags=[db_tags tag ':'];
                    db_access{db_id}=[strtrim(tline(27:32)) ':' strtrim(tline(34:41))];
                    sstart=sscanf(tline(56:60),'%i');
                    send=sscanf(tline(63:67),'%i');
                end;
                if ~isnan(sstart) && ~isnan(send),
                    db_cutout(db_id,1)=sstart;
                    db_cutout(db_id,2)=send;
                end;
            end;
        case 'MODEL '
            if length(tline)>=11,
                curr_model=sscanf(tline(11:end),'%i');
            else
                curr_model=1;
            end;                
            if curr_model>1, 
                break; 
            end;
        case {'ATOM  '}
            atom_tag=tline(13:16);
            atom_tag=strtrim(atom_tag);
            chain_tag=tline(22);
            if chain_tag==' '; 
                chain_tag=get_chain_tag(terminated_chains); 
            end;
            c_id=tag2id(chain_tag,chain_tags);
            if isempty(c_id), % catch cases without SEQRES records
                chain_id=chain_id+1;
                chain_tags=[chain_tags chain_tag ':'];
                chain_ids=[chain_ids chain_id];
                c_id=chain_id;
                protein(c_id).Ca=zeros(2000,3);
                protein(c_id).Cb=zeros(2000,3);
                protein(c_id).resnum=nan(1,2000);
                protein(c_id).restype=char(32*ones(1,2000));
            end;            
            if strcmpi(atom_tag,'CA'),
                curr_resname=tline(18:20);
                amino_id=tag2id(curr_resname,upper(residue_defs.restags),residue_defs.single_letter_code); % test whether this is an amino acid
                % store the atom coordinates and B-factor in the appropriate array
                x=sscanf(tline(31:38),'%f');
                y=sscanf(tline(39:46),'%f');
                z=sscanf(tline(47:54),'%f');
                coor=[x,y,z];
                curr_res=sscanf(tline(23:26),'%i');
                if isnan(curr_res), continue; end;
                ind=find(protein(c_id).resnum==curr_res);
                if isempty(ind),
                    pointers(c_id,1)=pointers(c_id,1)+1;
                    ind=pointers(c_id,1);
                    protein(c_id).resnum(ind)=curr_res;
                end;
                protein(c_id).Ca(ind,:)=coor;
                if ~isempty(amino_id),
                    protein(c_id).restype(ind)=amino_id;
                else
                    protein(c_id).restype(ind)='?';
                end;
            end;
            if strcmpi(atom_tag,'CB'),
                % store the atom coordinates and B-factor in the appropriate array
                x=sscanf(tline(31:38),'%f');
                y=sscanf(tline(39:46),'%f');
                z=sscanf(tline(47:54),'%f');
                coor=[x,y,z];
                curr_res=sscanf(tline(23:26),'%i');
                if isnan(curr_res), continue; end;
                ind=find(protein(c_id).resnum==curr_res);
                if isempty(ind),
                    pointers(c_id,1)=pointers(c_id,1)+1;
                    ind=pointers(c_id,1);
                    protein(c_id).resnum(ind)=curr_res;
                end;
                protein(c_id).Cb(ind,:)=coor;
            end;
        case 'END   '
            break;
    end;
end

fclose(fid);

% process data base reference for sequences

for k=1:length(chain_ids),
    cid=chain_ids(k);
    ctag=id2tag(cid,chain_tags);
    id=tag2id(ctag,db_tags);
    access=db_access{id};
    protein(k).dbref=access;
end;

% terminate arrays
for k=1:length(chain_ids),
    cid=chain_ids(k);
    poi=pointers(cid,1);
    protein(k).Ca=protein(k).Ca(1:poi,:);
    protein(k).Cb=protein(k).Cb(1:poi,:);
    protein(k).resnum=protein(k).resnum(1:poi);
    protein(k).restype=protein(k).restype(1:poi);
end;

delete(fname);