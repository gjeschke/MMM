function shortlist=residue_pair_score_basis(maxn,mask,target,threshold)
% function [shortlist,score]=residue_pair_score_basis(maxn,mask,target,threshold)
%
% determines maxn residue pairs that are predicted to be most useful
% for characterizing state transitions via low-frequency modes of an
% elastic network model, assuming fit mode with equipartitioning of thermal
% energy to the normal modes
% chain identifiers of template and target structure
% must correspond if a target structure is specified
%
% maxn      maximum number of residue pairs requested, at the same time
%           number of normal modes included in the basis 
% mask      (optional) mask of allowed residues, 0 disallowed, 1 allowed
%           indices into must correspond to indices into
%           model.coarse.Ca_coor,
%           can be empty
% target    (optional) MMM address of target structure, if present, mask is
%           automatically generated from template/target correspondence and
%           shortlist has a fourth column with target distances
% threshold threshold (minimum distance change) for including a site pair
%           into the constraint list, defined via the fraction of the maximum
%           distance change of any site pair for a given mode, to be
%           included, a pair must reach or exceed the threshold for at
%           least one mode
%           optional, defaults to 0.75
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

vec_norm=true; % minimum distance change threshold corresponds to norm of distance change vector
               % if true, to maximum change for at least one mode, if false

min_dist=1.7; % minimum and maximum distance [nm]
max_dist=5;
min_dist=min_dist*10; % working in Angstroem
max_dist=max_dist*10;

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

if nargin<4,
    if vec_norm,
        threshold=0.5;
    else
        threshold=0.5;
    end;
end;

[m,n]=size(model.ANM(model.current_structure).u);
network0=model.coarse(model.current_structure).Ca_coor;
frq=sqrt(model.ANM(model.current_structure).lambda);
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
        % ctag='A';
        adr=sprintf('%s(%s)%i.CA',target,ctag,resnum-200);
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

add_msg_board(sprintf('%i site pairs are considered.',(m/3)*(m/3-1)/2));
% make pair list
ind=0; % pair pointer
pairs=zeros((m/3)*(m/3-1)/2,2); % residue index pairs
rpairs=zeros((m/3)*(m/3-1)/2,1); % pair distances
pindices=zeros(size(mmask));
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
                rpairs(ind)=r;
                pindices(i,j)=ind;
            end;
        end;
    end;
end;
pairs=pairs(1:ind,:);
rpairs=rpairs(1:ind,:);
add_msg_board(sprintf('%i site pairs are within distance range.',ind));

drvecs=zeros(ind,maxn); % matrix of distance changes

add_msg_board('Determining best site pairs');
% set up the matrix of distance changes for all site pairs with respect to
% all modes in the basis
dmat0=coor2dmat(network0); % distance matrix for initial structure
for k=7:6+maxn,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    network=network0+mode'; % network change by move along this mode
    dmat=coor2dmat(network); % corresponding distance matrix 
    dmat1=dmat-dmat0; % distance change matrix
    for i=3:m/3-3,
        for j=i+1:m/3-2,
            if pindices(i,j),
                drvecs(pindices(i,j),k-6)=dmat1(i,j);
            end;
        end;
    end;
end;

if vec_norm,
    % Determine site pairs that are above the threshold for norm of distance
    % change vector
    above_thresh=zeros(1,ind);
    for k=1:ind,
        above_thresh(k)=sqrt(sum(drvecs(k,:).^2)); % norm of distance change vector
    end;
    above_thresh=above_thresh/max(above_thresh); % normalization
    selected=find(above_thresh>=threshold);
    thresh_ext=threshold;
    while length(selected)>36000,
        thresh_ext=thresh_ext+0.05;
        selected=find(above_thresh>=thresh_ext);
    end;
    while length(selected)<maxn && thresh_ext>0.05,
        thresh_ext=thresh_ext-0.05;
        selected=find(above_thresh>=thresh_ext);
    end;
else
    % Determine site pairs that are above the threshold for minimum distance
    % change with respect to at least one mode
    above_thresh=zeros(ind,maxn);
    for k=1:maxn,
        thr=threshold*max(abs(drvecs(:,k))); % minimum required distance change
        above_thresh(:,k)=abs(drvecs(:,k))>=thr;
    end;
    selected=find(sum(above_thresh,2)>0);
    thresh_ext=threshold;
    while length(selected)>37000,
        thresh_ext=thresh_ext+0.05;
        for k=1:maxn,
            thr=thresh_ext*max(abs(drvecs(:,k))); % minimum required distance change
            above_thresh(:,k)=abs(drvecs(:,k))>=thr;
        end;
        selected=find(sum(above_thresh,2)>0);
    end;
    while length(selected)<maxn && thresh_ext>0.05,
        thresh_ext=thresh_ext-0.05;
        for k=1:maxn,
            thr=thresh_ext*max(abs(drvecs(:,k))); % minimum required distance change
            above_thresh(:,k)=abs(drvecs(:,k))>=thr;
        end;
        selected=find(sum(above_thresh,2)>0);
    end;
end;
add_msg_board(sprintf('%i site pairs are above distance change threshold of %5.2f.',length(selected),thresh_ext));

drvecsn=zeros(length(selected),maxn);
% Normalize distance change vectors
for ks=1:length(selected),
    drvecsn(ks,:)=drvecs(selected(ks),:)/norm(drvecs(selected(ks),:));
end;

% Determine linear dependence of site pairs
mindep=1e6;
sellist=zeros(1,maxn);
selfom=sellist;
lindep=eye(length(selected));
for k1=1:length(selected)-1,
    for k2=k1+1:length(selected),
        lindep(k1,k2)=abs(sum(drvecsn(k1,:).*drvecsn(k2,:)));
        lindep(k2,k1)=lindep(k1,k2);
        % The first two site pairs are the one that are least dependent with
        % respect to the basis of modes
        if lindep(k1,k2)<mindep,
            mindep=lindep(k1,k2);
            sellist(1)=k1;
            sellist(2)=k2;
        end;
    end;
end;
fprintf(1,'Minimum linear dependence of first two site pairs: %5.3f\n',mindep);
selfom(2)=mindep;
poi=2;
% add additional site pairs that are least linearly dependent with respect
% to existing pairs
while poi<maxn,
    selmat=lindep(:,sellist(1:poi));
    selvec=max(selmat,[],2);
    for k=1:poi,
        selvec(sellist(k))=1;
    end;
    [mindep,sel]=min(selvec);
    poi=poi+1;
    sellist(poi)=sel;
    selfom(poi)=mindep;
    fprintf(1,'Maximum linear dependence of new site pair %i: %5.3f\n',poi,mindep);
end;

% Create short list from list of selections
shortlist=zeros(maxn,4);
for k=1:maxn,
    pair=selected(sellist(k));
    shortlist(k,1)=selfom(k);
    shortlist(k,2)=pairs(pair,1);
    shortlist(k,3)=pairs(pair,2);
    shortlist(k,4)=rpairs(pair);
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