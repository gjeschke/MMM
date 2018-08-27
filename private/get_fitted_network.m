function network=get_fitted_network(template,target,basis,direct,correspondence)
% G. Jeschke, 2012

if nargin<3,
    basis=20;
end;

basis=basis+6;

[rmsd0,coor]=rmsd_superimpose(target(correspondence(2,:),:),template(correspondence(1,:),:));
add_msg_board(sprintf('The r.m.s.d. between the two structures is %5.2f Å\n',rmsd0));
u=get_ANM(coor);
[m,m1]=size(u);
diff=target(correspondence(2,:),:)-coor;
[n,n1]=size(diff);
add_msg_board(sprintf('Structure pair has %i aligned Calpha atoms\n',n));
change=reshape(diff',1,3*n);
change=change';
coeff=zeros(m,1);
for k=1:m,
    coeff(k)=sum(u(:,k).*change);
end;
diff_real=u(:,1:basis)*coeff(1:basis);
diff_real=reshape(diff_real,3,m/3);
network=coor+diff_real';
rmsd=rmsd_superimpose(target(correspondence(2,:),:),network);
add_msg_board(sprintf('%i modes allow for minimum r.m.s.d. of %5.2f Å.\n',basis-6,rmsd));

network1=template;
network1(correspondence(1,:),:)=network;
if nargin>3,
    [md,nd]=size(direct); % m number of residues
    if md>0,
        rmsd_dir=0;
        for k=1:md, % set up unit vector between restraint residue pairs
            vec=network1(direct(k,2),:)-network1(direct(k,1),:);
            rp=norm(vec);
            rmsd_dir=rmsd_dir+(10*direct(k,4)-rp)^2;    
        end;
        rmsd_dir=sqrt(rmsd_dir/md);
        add_msg_board(sprintf('Best possible direct distance constraint r.m.s.d. %6.2f',rmsd_dir));
    end;
end;


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
