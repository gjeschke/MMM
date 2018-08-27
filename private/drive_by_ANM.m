function [rmsd,coor1] = drive_by_ANM(coor0,target,basis0,update_mode,thermal)
% update_mode   0 no mode update
%               1 reorientation
%               2 recomputation

max_step=0.2; % maxium coordinate change (Å)

if nargin<4,
    update_mode=2;
end;

if nargin<5,
    thermal=false;
end;

switch update_mode
    case 0
        update_string='original';
    case 1
        update_string='reoriented';
    case 2
        update_string='recomputed';
end;

iter=0;
[rmsd,coor]=rmsd_superimpose(target,coor0);

basis=basis0+6;

figure(1); clf;
set(gca,'FontSize',14);
title('Deviation (C^\alpha r.m.s.d.) from target structure');
h=line(iter,rmsd,'Color','k');

rmsd_diff=rmsd;
nit=0;

[u,lambda]=get_ANM(coor);
ur=u(:,1:basis);
if thermal,
    ut=u(:,7:end);
    lambdat=lambda(7:end);
    [mt,nt]=size(ut);
    phases=zeros(nt,1);
end;

nu=sqrt(lambda(7:end));
inv_nu=ones(size(nu))./nu;
inv_nu=inv_nu/inv_nu(1);
[mi,maxbas]=min(abs(inv_nu-0.1));

if maxbas<basis0,
    minbas=maxbas;
else
    minbas=basis0;
end;

basis_size=maxbas*ones(1,400);
basis_size(1:100)=linspace(minbas,maxbas,100);
basis_size=round(basis_size);

scsum=0;
while (nit<400 && rmsd_diff > 5e-6) || nit<10,
    nit=nit+1;
    diff=target-coor;
    [m,m1]=size(u);
    [n,n1]=size(diff);
    change=reshape(diff',1,3*n);
    if ~thermal,
        coeff=ur\change';
    else
        for k=1:nt,
            uk=ut(:,k);
            uk=uk/norm(uk);
            nchange=change'/norm(change);
            phases(k)=sum(uk.*nchange);
        end;
    end;
    if nit==1 && ~thermal,
        fprintf(1,'Basis is reduced to %i modes, corresponding to %4.1f%%.\n',basis0,100*basis0/(m-6));
        diff_test=zeros(n,3);
        for k=1:basis,
            evec=u(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test=diff_test+mode1'*coeff(k);
        end;
        coor_b=coor+diff_test;
        rmsd_f=rmsd_superimpose(target,coor0);
        rmsd_e=rmsd_superimpose(target,coor_b);
        fprintf(1,'%i modes account for %4.2f%% of change. Deviation is %4.2f Å\n',basis0,100*(rmsd_f-rmsd_e)/rmsd_f,rmsd_e);
    end;
%     step=zeros(n,3);
%     for k=1:basis,
%         evec=u(:,k);
%         mode1=reshape(evec,3,m/3);
%         step=step+mode1'*coeff(k);
%     end;
    if ~thermal,
        step=ur*coeff;
    else
        coeff=sqrt(ones(size(lambdat))./lambdat);
        step=ut(:,1:basis_size(nit))*(phases(1:basis_size(nit)).*coeff(1:basis_size(nit)));
    end;
    step=reshape(step,3,m/3);
    step=step';
    per_residue=sqrt(sum(step.^2,2));
    sc=max_step/max(per_residue);
    scsum=scsum+sc;
    if sc>1, sc=1; end;
    coorm=coor;
    coor=coor+sc*step;
    crmsd=rmsd_superimpose(target,coor);
    iter=[iter nit];
    rmsd=[rmsd crmsd];
    if nit>10,
        rmsd_diff=(rmsd(end-10)-rmsd(end))/10;
    end;
    set(h,'XData',iter,'YData',rmsd);
    drawnow;
    if update_mode==2,
        [u,lambda]=get_ANM(coor);
        if ~thermal,
            ur=u(:,1:basis);
        else
            ut=u(:,7:end);
            lambdat=lambda(7:end);
        end;
    end;
    if update_mode==1,
        if ~thermal,
            ur=reorientate_eigenvectors(coor,coorm,ur);
        else
            ut=reorientate_eigenvectors(coor,coorm,ut);
        end;
    end;
end;
fprintf(1,'Sum of scaling coefficients over path is: %6.4f\n',scsum);
axis([0,max(iter),0,1.1*max(rmsd)]);
coor1=coor;
rmsd_f=rmsd_superimpose(target,coor0);
rmsd_e=rmsd_superimpose(target,coor1);
fprintf(1,'%i %s modes account for fraction %5.3f of change. Deviation is %4.2f Å\n',basis0,update_string,(rmsd_f-rmsd_e)/rmsd_f,rmsd_e);
fprintf(1,'%i iterations were required.\n',nit);


function [u,lambda]=get_ANM(coor)

Hessian=setup_ANM_bonded(coor);
[u,d]=eig(Hessian);
clear Hessian
lambda=diag(d);

function Hessian=setup_ANM_bonded(Ca_coor)
% function Hessian=setup_ANM_bonded(Ca_coor)
%
% returns the Hessian for a fully connected anisotropic network model (ANM)
% when given the Calpha coordinates
% the force constants have inverse 6th power dependence on interresidue
% distance, except for next neighbors (i,i+1) and second neighbors (i,i+2)
% along the polypetide chain
%
% ### This version assumes a single-chain peptide without gaps in the sequence ###
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% and:
% L. Yang, G. Song, R L. Jernigan, Proc. Natl. Acad. Sci. 106 (2009)
% 12347–12352.
%
% the force constant at a distance of 1 Å p_ANM_gamma is defined in
% initialize_MMM.m  and stored in global variable ENM_param
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param

bond_force_1=10000*3.8^(-6); % 1200
bond_force_2=10000*6^(-6); % 80
exponent=6;

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Hessian=zeros(3*m);
r1=zeros(1,m);
p1=0;
r2=zeros(1,m);
p2=0;
for i=1:m-1,
    basi=3*(i-1);
    for j=i+1:m,
        basj=3*(j-1);
        rvec=Ca_coor(i,:)-Ca_coor(j,:); % distance vector of Calpha atoms
        r=norm(rvec);
        diff=i-j;
        switch diff
            case 1
                cgamma=-bond_force_1*ENM_param.p_ANM_gamma;
                p1=p1+1;
                r1(p1)=r;
            case 2
                cgamma=-bond_force_2*ENM_param.p_ANM_gamma;
                p2=p2+1;
                r2(p2)=r;
            otherwise
                cgamma=-ENM_param.p_ANM_gamma*r^(-exponent);
        end;
        if cgamma~=0,
            submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
            Hessian(basi+1:basi+3,basj+1:basj+3)=cgamma*submat;
            Hessian(basj+1:basj+3,basi+1:basi+3)=cgamma*submat;
        end;
    end;
end;

conn=0;
for i=1:m,
    basi=3*(i-1);
    submat=zeros(3);
    for j=1:m,
        basj=3*(j-1);
        submat=submat+Hessian(basi+1:basi+3,basj+1:basj+3);
    end;
    conn=conn+trace(submat);
    Hessian(basi+1:basi+3,basi+1:basi+3)=-submat;
end;

if conn==0, % ANM is unconnected
    Hessian=[];
end;

function basis=reorientate_eigenvectors(network,network0,basis0)

[mb,nb]=size(basis0);
basis=basis0;
dmat=coor2dmat(network);
[m,n]=size(dmat);
dr=network-network0;
dr=sum(dr.^2,2);
for k=1:m, % loop over residues
    basn=3*(k-1);
    distances=dmat(k,:); % distances to other residues
    distances(k)=1e6; % distances to 
    candidates=find(distances<=10);
    change=dr(candidates);
    [schange,poi]=sort(change);
    loc=4;
    if length(candidates)<loc,
        loc=length(candidates);
    end;
    rigid=[candidates(poi(1:loc)) k];
    coor0=network0(rigid,:);
    coor=network(rigid,:);
    [rms,coor0b,transmat]=rmsd_superimpose(coor,coor0);
    transmat(1:3,4)=[0;0;0];
    for kk=1:nb,
        xyz=[basis0(basn+1:basn+3,kk);1];
        xyz=transmat*xyz;
        basis(basn+1:basn+3,kk)=xyz(1:3);
    end;
end;
