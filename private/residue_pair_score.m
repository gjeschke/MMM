function [shortlist,score,redundant,pairs,cutoff]=residue_pair_score(maxn,corr_cut,slowmodes,mask,target)
% function [shortlist,score,redundant]=residue_pair_score(maxn,corr_cutoff,slowmodes,mask,target)
%
% determines up to maxn residue pairs that are predicted to be most useful
% for characterizing state transitions via low-frequency modes of an
% elastic network model, chain identifiers of template and target structure
% must correspond if a target structure is specified
%
% maxn      maximum number of residue pairs requested
% corr_cut  optional cutoff value for correlation between radical pairs
%           that are considered redundant
% slowmodes number of low frequency modes to consider in fitting
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
% score     relative expected mean square distance changes for propagation
%           along low-frequency modes, sorted vector (descending) for
%           residue pairs
% redundant vector that shows whether a residue pair is (largely) redundant
%           with a previous one, if so the number of the previous (better)
%           pair is given, if not a zero is given
% pairs     indices into model.coarse(model.current_structure).Ca_coor for
%           the site pairs
% cutoff    correlation cutoff at which the requested (or maximum) number
%           of pairs was generated
%
% the total number of returned pairs nn can be larger than the number maxn
% of requested non-redundant pairs
%
% if target structure is specified, but not found, the function returns
% empty output
%
% for algorithm, see:
% W. Zheng, B. R. Brooks, Biophys. J. 2006, 90: 4327-4336, p.4328
%
% G. Jeschke, 2010-2012


global model

maxlow=10; % maximum number of low-frequency modes
corr_cutoff=0.3; % correlation cutoff for considering residue pair as redundant 

if nargin>1,
    corr_cutoff=corr_cut;
end;

if nargin>2,
    maxlow=slowmodes;
end;

tnum=[];
if nargin>4,
    [tnum,msg]=resolve_address(target);
    if isempty(tnum) || msg.error,
        shortlist=[];
        score=[];
        redundant=[];
        pairs=[];
        add_msg_board('ERROR: Target structure was specidied but not found.');
        add_msg_board('No constraints generated.');
        return
    end;
end;

[m,n]=size(model.ANM(model.current_structure).u);
network0=model.coarse(model.current_structure).Ca_coor;
network=network0;
scal=sqrt(model.ANM(model.current_structure).lambda(7));
for k=7:6+maxlow,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    network=network+scal*mode'/sqrt(model.ANM(model.current_structure).lambda(k));
end;
diff=zeros((m/3)*(m/3-1)/2,1);
pairs=zeros((m/3)*(m/3-1)/2,2);
rpairs=zeros((m/3)*(m/3-1)/2,1);
ind=0;
if nargin<4 || isempty(mask),
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

add_msg_board('Setting up distance changes');
for i=3:m/3-3,
    for j=i+1:m/3-2,
        if mask(i) && mask(j),
            ind=ind+1;
            pairs(ind,1)=i;
            pairs(ind,2)=j;
            r0=norm(network0(i,:)-network0(j,:));
            r=norm(network(i,:)-network(j,:));
            diff(ind)=(r-r0);
            rpairs(ind)=r0;
        end;
    end;
end;
diff=diff(1:ind);
pairs=pairs(1:ind,:);
add_msg_board('Sorting by distance change');
[diff1,poi]=sort(diff.^2,'descend');
if maxn>ind,
    maxn=ind;
end;
pointers=zeros(1,maxn);
rijm=zeros(maxn,maxlow);
rijm0=zeros(1,maxlow);
slowmodes=zeros(maxlow,m/3,3);
for k=7:6+maxlow,
    evec=model.ANM(model.current_structure).u(:,k);
    mode=reshape(evec,3,m/3);
    slowmodes(k-6,:,:)=mode';
end;

add_msg_board('Finding pairs...');
n_old=0;
n=0;
% scal=model.ANM(model.current_structure).lambda(7);
while n<maxn && corr_cutoff<=0.76,
    n=0;
    nn=0;
    score=zeros(1,(m/3)*(m/3-1)/2);
    redundant=zeros(1,(m/3)*(m/3-1)/2);
    shortlist=zeros(maxn,4);
    while n<maxn && nn <ind,
        nn=nn+1;
        score(nn)=diff1(nn);
        for k=7:6+maxlow,
            mode=reshape(slowmodes(k-6,:,:),m/3,3);
            i=pairs(poi(nn),1);
            j=pairs(poi(nn),2);
            r0=norm(network0(i,:)-network0(j,:));
            r=norm(network(i,:)+mode(i,:)-network(j,:)-mode(j,:));
            rijm0(k-6)=(r-r0)/sqrt(model.ANM(model.current_structure).lambda(k));
        end;
       rijm0=rijm0/norm(rijm0);
        if nn>1 && n>=1,
            maxcorr=0;
            for j=1:n,
                num=sum(rijm0.*rijm(j,:));
                denom=sum(rijm0.^2);
                corr=abs(num/denom);
                if corr>corr_cutoff && corr>maxcorr,
                    redundant(nn)=j;
                    maxcorr=corr;
                end;
            end;
        end;
        if redundant(nn)==0,
            n=n+1;
            add_msg_board(sprintf('Cutoff %5.2f: %i pair(s) found out of %i tested.\n',corr_cutoff,n,nn));
            shortlist(n,1)=score(nn);
            shortlist(n,2:3)=pairs(poi(nn),:);
            shortlist(n,4)=rpairs(poi(nn));
            pointers(n)=poi(nn);
            rijm(n,:)=rijm0;
        end;
    end;
    if n>n_old,
        shortlist_old=shortlist;
        score_old=score;
        redundant_old=redundant;
        pairs_old=pairs;
        n_old=n;
        cutoff_old=corr_cutoff;
    else
        shortlist=shortlist_old;
        score=score_old;
        redundant=redundant_old;
        pairs=pairs_old;
        n=n_old;
        add_msg_board('Number of pairs did not increase.');
        add_msg_board(sprintf('Reverting to solution for correlation cutoff %5.2f\n',cutoff_old));
    end;
    corr_cutoff=corr_cutoff+0.05;
end;
shortlist=shortlist(1:n,:);
score=score(1:nn);
redundant=redundant(1:nn);
pairs=pairs(1:nn,:);
cutoff=cutoff_old;

if ~isempty(tnum),
    add_msg_board('Extracting target distances');
    for k=1:n,
        i1=shortlist(k,2);
        i2=shortlist(k,3);
        xyz1=target_Ca_xyz(i1,:);
        xyz2=target_Ca_xyz(i2,:);
        r12=norm(xyz1-xyz2);
        shortlist(k,4)=r12;
    end;
end;