function [coor,col]=get_nucleotide_coor(indices)
% function [coor,col]=get_nucleotide_coor(indices)
%   returns coordinate array coor for a nucleotide residue in RNA or DNA
%   together with the specified default color 
%
% indices   index vector with at least 4 indices
%
% coor      nx3 coordinate array with following assignements
%           coor(1:2,:)     two backbone atoms (C4', C3')
%           coor(3:n,:)     edge polygon of base, starting at connection to
%                           ribose
% col       RGB vector in range 0...1 for standard color of nucleotide
%
% empty array is returned if indices do not address a single residue, 
% residue is not a nucleotide, or definitions are missing
%
% G. Jeschke, 2010

global model
global residue_defs
global graph_settings

coor=[];
col = [];

[m,n]=size(indices);
if n<4 || m~=1, % not a single residue addressed
    return;
end;

info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4));
if info.type~=2, % not a nucleotide
    add_msg_board('Residue is not a nucleotide');
    return;
end;

resname=info.name;
id=tag2id(resname,residue_defs.nucleotide_tags);
if isempty(id), % definitions missing
    add_msg_board(sprintf('Definitions missing for residue :%s.',resname));
    return; 
end;
nucleotide=residue_defs.nucleotides(id);

coor=zeros(20,3);

nonsense=textscan(nucleotide.backbone,'%s','Delimiter',':');
atag_list=nonsense{1};

if length(atag_list)~=3,
    disp('Backbone definition mismatch!');
    return
end;

poi=0;
for k=2:length(atag_list),
    curr_id=tag2id(char(atag_list{k}),info.atom_tags);
    curr_atnum=info.atom_numbers{curr_id};
    poi=poi+1;
    [mm,nn]=size(curr_atnum);
    if nn==1, curr_atnum=[curr_atnum 1]; end;
    ccoor=zeros(1,3);
    pop=0;
    for kk=1:mm,
        ccoor=ccoor+curr_atnum(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(curr_atnum(kk,1),:);
        pop=pop+curr_atnum(kk,2);
    end;
    coor(poi,:)=ccoor/pop;
end;

nonsense=textscan(nucleotide.polygon,'%s','Delimiter',':');
atag_list=nonsense{1};

if length(atag_list)<4,
    disp('Polygon definition mismatch!');
    return
end;

for k=2:length(atag_list),
    curr_id=tag2id(char(atag_list{k}),info.atom_tags);
    curr_atnum=info.atom_numbers{curr_id};
    poi=poi+1;
    [mm,nn]=size(curr_atnum);
    if nn==1, curr_atnum=[curr_atnum 1]; end;
    ccoor=zeros(1,3);
    pop=0;
    for kk=1:mm,
        ccoor=ccoor+curr_atnum(kk,2)*model.structures{indices(1)}(indices(2)).xyz{indices(3)}(curr_atnum(kk,1),:);
        pop=pop+curr_atnum(kk,2);
    end;
    coor(poi,:)=ccoor/pop;
end;

coor=coor(1:poi,:);

colind=tag2id(nucleotide.color,graph_settings.color_tags);
if ~isempty(colind),
    col=graph_settings.colors(colind,:);
else
    col=[192,192,192]; % silver as default
end;
