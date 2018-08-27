function shortlist=residue_pair_score_thermal(maxn,mask,target)
% function [shortlist,score]=residue_pair_score(maxn,,mask,target)
%
% determines up to maxn residue pairs that are predicted to be most useful
% for characterizing state transitions via low-frequency modes of an
% elastic network model, assuming fit mode with equipartitioning of thermal
% energy to the normal modes
% chain identifiers of template and target structure
% must correspond if a target structure is specified
%
% maxn      maximum number of residue pairs requested
% mask      (optional) mask of allowed residues, 0 disallowed, 1 allowed
%           indices into must correspond to indices into
%           model.coarse.Ca_coor,
%           can be empty
% target    (optional) MMM address of target structure, if present, mask is
%           automatically generated from template/target correspondence and
%           shortlist has a fourth column with target distances
%
% shortlist a short list of non-redundant residue pairs sorted by
%           descending score, 1st column: score, 2nd/3rd column residue
%           numbers, 4th column: Calpha-Calpha distance in target, if
%           target is given, otherwise distance change for small
%           deformation
%           structure, if input argument target was present
% if target structure is specified, but not found, the function returns
% empty output
%
%
% G. Jeschke, 2012


global model

tnum=[];
if nargin>2,
    [tnum,msg]=resolve_address(target);
    if isempty(tnum) || msg.error,
        shortlist=[];
        add_msg_board('ERROR: Target structure was specified but not found.');
        add_msg_board('No constraints generated.');
        return
    end;
end;

[m,n]=size(model.ANM(model.current_structure).u);
network0=model.coarse(model.current_structure).Ca_coor;
if nargin<2 || isempty(mask),
    mask=ones(1,m/3);
    mask(1:2)=0;
    mask(end-1:end)=0;
end;

if ~isempty(tnum),
    add_msg_board('Determining template/target correspondence.');
    rindices=model.coarse(model.current_structure).indices;
    [mr,nr]=size(rindices);
    target_Ca_xyz=zeros(mr,3);
    for k=1:mr,
        [stag,ctag,modelnum,resnum]=mk_address_parts(rindices(k,:));
        adr=sprintf('%s(%s)%i.CA',target,ctag,resnum);
        [cindices,msg]=resolve_address(adr);
        if isempty(cindices) || msg.error,
            mask(k)=0;
        else
            [msg,xyz]=get_residue(cindices,'xyz');
            if msg.error || isempty(xyz),
                mask(k)=0;
            else
                [mx,nx]=size(xyz);
                if mx>1,
                    xyz=mean(xyz);
                end;
                target_Ca_xyz(k,:)=xyz;
            end;
        end;
    end;
    figure(3); clf;
    plot(mask);
    title('Correspondence mask');
end;
mmask=kron(mask,mask');

add_msg_board('Determining best site pairs');
shortlist=zeros(maxn,4);
poi=0;
dmat0=coor2dmat(network0); % distance matrix for initial structure
for k=7:6+maxn,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    network=network0+mode'; % network change by move along this mode
    dmat=coor2dmat(network); % corresponding distance matrix 
    dmat1=dmat-dmat0; % distance change matrix
    dmat=abs(dmat1).*mmask; % mask all pairs that should not be selected at all
    for k0=1:poi, % mask all pairs that were selected before and very close pairs
        k1=shortlist(k0,2);
        k2=shortlist(k0,3);
        dmat(k1-2:k1+2,k2-2:k2+2)=0;
        dmat(k2-2:k2+2,k1-2:k1+2)=0;
    end;
    [ma,k0]=max(dmat); % determine site pair that corresponds to 
    [m2,k2]=max(ma);   % maximum distance change for this mode
    k1=k0(k2);
    poi=poi+1;
    shortlist(poi,1)=dmat1(k1,k2);
    shortlist(poi,2)=k1;
    shortlist(poi,3)=k2;
    shortlist(poi,4)=dmat0(k1,k2);
end;


if ~isempty(tnum),
    add_msg_board('Extracting target distances');
    for k=1:maxn,
        i1=shortlist(k,2);
        i2=shortlist(k,3);
        xyz1=target_Ca_xyz(i1,:);
        xyz2=target_Ca_xyz(i2,:);
        r12=norm(xyz1-xyz2);
        shortlist(k,4)=r12;
    end;
end;