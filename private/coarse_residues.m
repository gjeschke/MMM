function [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(address,modnum)
% function [Ca_coor,masses,rindices,Bfactors,restypes]=coarse_residues(address,modelnum)
%
% extracts Calpha coordinates and residue masses for an elastic network
% model from an MMM atomistic model
% the ENM can encompass a whole structure, a chain, a chain model 
% (in structures with several models), or a range of residues as specified 
% by address
% instead of the address, an MMM index vector can be provided
%
% if a chain is addressed, for which several models (coordinate sets) exist, 
% the first model is used, unless argument modnum is given
%
% note that modnum is also used, if an MMM index vector is provided with
% length < 3
%
% if argument modnum is empty, model 1 is addressed
%
% empty arrays are returned if the address is invalid
% if a Ca atom has alternate locations, an average coordinate is returned
% any residues other than amino acids are ignored
% if the structure does not contain hydrogen atoms, only heavy-atom masses
% are considered
%
% indices of the residues and Bfactors are returned for later reference
% residue types (internal number codes for the amino acids) are returned in
% restypes, unknown amino acids have type 0
%
%
% G. Jeschke, 2010-2012

global model
global residue_defs

if nargin<2
    modnum=1;
end

if isempty(modnum)
    modnum=1;
end

Ca_coor=[];
rindices=[];
masses=[];
Bfactors=[];
restypes=[];

% check hierarchy level that is addressed
if isa(address,'char')
    indices=resolve_address(address); 
else
    indices=address;
    address=mk_address(indices);
end
if isempty(indices), return; end % invalid address
[~,n]=size(indices); % n is the hierarchy level

% extend address appropriately
switch n
    case 1 % structure
        address=sprintf('%s(:){%i}.CA',address,modnum);
    case 2 % chain
        address=sprintf('%s{%i}.CA',address,modnum);
    case {3,4} % chain model or range of residues
        address=[address '.CA'];
end
% fprintf(1,'Addressing: %s.\n',address);

% create index array of all residues
indices=resolve_address(address);
if isempty(indices), return; end % may happen for DNA, RNA, generally non-proteins
[m,~]=size(indices); % m is the number of residues

% extract coordinates and masses
Ca_coor=zeros(m,3);
rindices=zeros(m,4);
masses=zeros(1,m);
Bfactors=zeros(1,m);
restypes=zeros(1,m);
poi=0;
for k=1:m
    info=model.structures{indices(k,1)}(indices(k,2)).residues{indices(k,3)}.info(indices(k,4));
    if info.type==1 % only amino acids are considered
        poi=poi+1;
        [~,coor]=get_atom(indices(k,:),'coor');
        Ca_coor(poi,:)=coor;
        [~,mass]=get_residue(indices(k,1:4),'mass');
        masses(poi,:)=mass;
        [~,Bfactor]=get_atom(indices(k,:),'Bfactor');
        Bfactors(poi,:)=mean(Bfactor);
        rindices(poi,:)=indices(k,1:4);
        id=tag2id(info.name,upper(residue_defs.restags));
        if ~isempty(id)
            restypes(poi)=id;
        end
    end
end
if poi==0 % no amino acid residues found
    Ca_coor=[];
    masses=[];
    rindices=[];
    Bfactors=[];
    restypes=[];
    return;
end

% restrict arrays to actual number of amino acid residues
Ca_coor = Ca_coor(1:poi,:);
rindices = rindices(1:poi,:);
masses = masses(1:poi);
Bfactors=Bfactors(1:poi);
restypes=restypes(1:poi);



