function shortlist=residue_pair_score_vectors(maxn,mask,target)
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
%
%           score for the first pair is distance change in first mode
%           score for other pairs is a figure of merit that combines
%           sum of sines of angles between this pair vector and all other
%           pair vectors and sum of minimum distances between this pair
%           vector and all other pair vectors with minimum distances
%           normalized by radius of gyration of the protein
%
% if target structure is specified, but not found, the function returns
% empty output
%
%
% G. Jeschke, 2012


global model


min_dist=1.7; %minimum and maximum distance [nm]
max_dist=5;
min_dist=min_dist*10; % working in Angstroem
max_dist=max_dist*10;

svec=zeros(maxn,3);
svecn=svec;

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
% mask the first two and last two residues
if nargin<2 || isempty(mask),
    mask=ones(1,m/3);
    mask(1:2)=0;
    mask(end-1:end)=0;
end;

% determine radius of gyration
[mc,nc]=size(network0);
meanc=mean(network0,1);
centered=network0-repmat(meanc,mc,1);
dr2=sum(sum(centered.^2));
rg=sqrt(dr2/mc);
fprintf(1,'Radius of gyration: %4.2f A\n',rg);



if ~isempty(tnum),
    add_msg_board('Determining template/target correspondence.');
    rindices=model.coarse(model.current_structure).indices;
    [mr,nr]=size(rindices);
    target_Ca_xyz=zeros(mr,3);
    for k=1:mr,
        [stag,ctag,modelnum,resnum]=mk_address_parts(rindices(k,:));
        ctag='A';
        %resnum=resnum-200;
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
evec=model.ANM(model.current_structure).u(:,7);
mode=reshape(evec,3,m/3);
network=network0+mode'; % network change by move along this mode
dmat=coor2dmat(network); % corresponding distance matrix 
dmat1=dmat-dmat0; % distance change matrix
dmat=abs(dmat1).*mmask; % mask all pairs that should not be selected at all
dmat=dmat.*(dmat0>=min_dist); % mask all distances that are too short
dmat=dmat.*(dmat0<=max_dist); % mask all distances that are too long
% First pair is the one with largest distance change in the first
% significant normal mode
[ma,k0]=max(dmat); % determine site pair that corresponds to 
[m2,k2]=max(ma);   % maximum distance change for this mode
k1=k0(k2);
poi=poi+1;
shortlist(poi,1)=m2;
shortlist(poi,2)=k1;
shortlist(poi,3)=k2;
shortlist(poi,4)=dmat0(k1,k2);
svec(1,:)=network0(k2,:)-network0(k1,:);
svecn(1,:)=svec(1,:)/norm(svec(1,:));

% make pair list
ind=0; % pair pointer
pvec=zeros((m/3)*(m/3-1)/2,3); % pair distance vectors
pa=zeros((m/3)*(m/3-1)/2,3); % pair starting points
pe=zeros((m/3)*(m/3-1)/2,3); % pair end points
pairs=zeros((m/3)*(m/3-1)/2,2); % residue index pairs
rpairs=zeros((m/3)*(m/3-1)/2,1); % pair distances
pindices=zeros(1,maxn);
add_msg_board('Setting up distance change vectors');
for i=3:m/3-3,
    for j=i+1:m/3-2,
        if mmask(i,j),
            rvec=network0(i,:)-network0(j,:);
            r=norm(rvec);
            if r>=min_dist && r<=max_dist,
                ind=ind+1;
                pairs(ind,1)=i;
                pairs(ind,2)=j;
                pa(ind,:)=network0(i,:);
                pe(ind,:)=network0(j,:);
                pvec(ind,:)=rvec;
                rpairs(ind)=r;
                if i==k1 && j==k2,
                    pindices(1)=ind;
                elseif i==k2 && j==k1,
                    pindices(1)=ind;
                end;
            end;
        end;
    end;
end;
pairs=pairs(1:ind,:);
rpairs=rpairs(1:ind,:);
pvec=pvec(1:ind,:);
pa=pa(1:ind,:);
pe=pe(1:ind,:);
dvec=sqrt(sum(pvec.^2,2)); % pair distances
pvecn=pvec./repmat(dvec,1,3); % unit vectors along pairs
merit=zeros(ind,maxn);  % matrix of figure of merit of all pairs with 
                        % respect to the pairs already selected
merit2=zeros(ind,1);

% make additional pairs that are linearly independent
for k=1:maxn-1, % index relates to already selected pair for which the
                % figure of merit is computed
    % cosine of angle between all pair vectors and vector of selected pair
    dircos=sum(pvecn.*repmat(svecn(k,:),ind,1),2);
    dirsin=sqrt(1-dircos.^2);
    for kk=1:ind,
        merit2(kk)=vec_dist(pa(pairs(kk,1),:),pe(pairs(kk,1),:),...
            pa(pairs(kk,2),:),pe(pairs(kk,2),:))/rg;
    end;
    merit(:,k)=dirsin+merit2; % full figure of merit with respect to pair k
    fom=prod(merit(:,1:k),2); % total figure of merit 
    pindsel=pindices(pindices>0); % indices of selected pairs in pair list
    fom(pindsel)=0; % make sure that a pair is not selected for a second
                    % time by setting figure of merit of selected pairs to
                    % zero
    [mfom,sel]=max(fom);
    shortlist(k+1,1)=mfom;
    shortlist(k+1,2)=pairs(sel,1);
    shortlist(k+1,3)=pairs(sel,2);
    shortlist(k+1,4)=dvec(sel);
    svec(k+1,:)=network0(pairs(sel,2),:)-network0(pairs(sel,1),:);
    svecn(k+1,:)=svec(k+1,:)/norm(svec(k+1,:));
    pindices(k+1)=sel;
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
