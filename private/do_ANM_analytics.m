function do_ANM_analytics(template,target)
% G. Jeschke, 2012

check_size=20;

clc;

[u,lambda]=get_ANM(template);
[m,m1]=size(u);
diff=target-template;
[n,n1]=size(diff);
[rmsd0,coor]=rmsd_superimpose(target,template);
fprintf(1,'Structure pair has %i aligned Calpha atoms\n',n);
fprintf(1,'The r.m.s.d. between the two structures is %5.2f Å\n',rmsd0);
change=reshape(diff',1,3*n);
change=change';
coeff=zeros(m,1);
for k=1:m,
    coeff(k)=sum(u(:,k).*change);
end;

xax=1:m-6;
acoeff=abs(coeff(7:end));
acoeff=acoeff/sum(acoeff);
figure(3); clf;
plot(xax,acoeff,'k');
hold on;
plot([20,20],[0,max(acoeff)],'k:');
axis([0,100,0,0.13]);
set(gca,'FontSize',14);

kappa = mode_collectivity(u);
mean_kappa=sum(kappa./lambda')/sum(ones(size(lambda))./lambda);
fprintf(1,'Normal modes have mean collectivity %5.3f\n',mean_kappa);
[kappas,poi_coll]=sort(kappa,'descend');
[coeffs,poi_ideal]=sort(abs(coeff),'descend');
kappa_DR=mode_collectivity(change);
fprintf(1,'Structural change has collectivity %5.3f\n',kappa_DR);
nmode_ax=7:m;
rmsd_ideal=zeros(size(nmode_ax));
rmsd_lowfrq=rmsd_ideal;
rmsd_collectivity=rmsd_ideal;
collectivity=rmsd_ideal;
for basis=nmode_ax,
    diff_ideal=u(:,poi_ideal(1:basis))*coeff(poi_ideal(1:basis));
    diff_ideal=reshape(diff_ideal,3,m/3);
    coor_ideal=template+diff_ideal';
    rmsd_ideal(basis-6)=rmsd_superimpose(target,coor_ideal);
    diff_coll=u(:,poi_coll(1:basis))*coeff(poi_coll(1:basis));
    diff_coll=reshape(diff_coll,3,m/3);
    coor_coll=template+diff_coll';
    rmsd_collectivity(basis-6)=rmsd_superimpose(target,coor_coll);
    diff_real=u(:,1:basis)*coeff(1:basis);
    collectivity(basis-6)=mode_collectivity(diff_real);
    diff_real=reshape(diff_real,3,m/3);
    coor_real=template+diff_real';
    rmsd_lowfrq(basis-6)=rmsd_superimpose(target,coor_real);
end;
nmode_ax=[6 nmode_ax];
rmsd_ideal=[rmsd0 rmsd_ideal];
rmsd_lowfrq=[rmsd0 rmsd_lowfrq];
rmsd_collectivity=[rmsd0 rmsd_collectivity];
collectivity=[0 collectivity];
rdecay=real(rmsd_ideal)/rmsd0;
bdecay=real(rmsd_lowfrq)/rmsd0;
figure(7);
plot(rdecay);
see=rdecay-1/2;
[~,poi]=min(abs(see));
% v0=[poi,1];
% v=fminsearch(@rms_stretched_exp,v0,[],rdecay);
fprintf(1,'Half change is covered by %i best normal modes.\n',poi-1);
% fprintf(1,'Stretch exponent is %4.2f\n',v(2));
% arg=0:length(rdecay)-1;
% arg=(arg/v(1)).^v(2);
% fct=rmsd0*exp(-arg);
see=bdecay-2/3;
[~,poi]=min(abs(see));
fprintf(1,'1/3 of the structural change is fitted by %i lowest normal modes.\n',round(poi)-1);
see=bdecay-1/2;
[~,poi]=min(abs(see));
fprintf(1,'1/2 of the structural change is fitted by %i lowest normal modes.\n',round(poi)-1);
see=rmsd_lowfrq-2;
[~,poi]=min(abs(see));
fprintf(1,'2 Å r.m.s.d. is achieved with %i lowest normal modes.\n',round(poi)-1);
fprintf(1,'%i lowest normal modes fit a fraction of %5.3f of the structural change.\n',check_size,1-bdecay(check_size+1));
fprintf(1,'These %i modes fit the final structure with an r.m.s.d. of %4.2f Å.\n',check_size,rmsd_lowfrq(check_size+1));
fprintf(1,'10 lowest normal modes fit a fraction of %5.3f of the structural change.\n',1-bdecay(11));


figure(1); clf;
plot(nmode_ax,real(rmsd_ideal),'k');
hold on;
% plot(nmode_ax,fct,'g:');
plot(nmode_ax,real(rmsd_lowfrq),'b');
plot(nmode_ax,real(rmsd_collectivity),'r');
title('Convergence behavior');
figure(2); clf;
plot(nmode_ax,real(collectivity),'k');
hold on;
plot([nmode_ax(1),nmode_ax(end)],[kappa_DR,kappa_DR],'r:');
title('Collectivity convergence');



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

function rms = rms_stretched_exp(v,decay)
%rms_stretched_exp	Root mean square error of function exp(-(n/N)^(ksi)).
% function rms = rms_stretched_exp(v,decay)
%	
%  Parameter: v(1) N decay constant (as point index)
%             v(2) ksi stretch exponent
%  	      decay decay function to be fitted
%             is assumed to be normalized at firts point decay(1)=1
%             and to decay to zero for infinite point index

%	G. Jeschke, 2012

if v(1)<0, rms=1.0e6; return; end;
arg=0:length(decay)-1;
arg=(arg/v(1)).^v(2);
fct=exp(-arg);

diff=fct-decay;
rms=sqrt(sum(diff.^2)/length(diff));
