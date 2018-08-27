function [fom,rmsd,network,dxmat] = fit_by_ANM_1(network0,restr)
% function [fom,rmsd,network,dxmat] = fit_by_ANM_1(network0,restr)
%
% fits a network model to distance restraints based on the algorithm from
% W. Zheng, B. R. Brooks,Biophys. J. 90 (2006) 4327-4336.
%
% the number of slow modes to be used for the transition is defined in
% global variable ENM_param.fit_basis
% specifics of local structure restraints to fix C_alpha-C_alpha distances
% or to stabilize secondary structure are defined in global variable
% ENM_param.fix_local
% see initialize_MMM.m for initialization of these variables
%
% Input:
% network0  matrix [n,3] of coordinates of network points (knots),
%           restraint variable DEER must refer labelled sites to network
%           knots (C_alpha coordinates of residues), these references are
%           given for the k-th restraint by DEER(k).res1 and DEER(k).res2
%           indexing of residues in network is supposed to conform to
%           model.coarse(model.current_structure).indices
% restr     restraint variable, n-by-4 or n-by-5 array
%           first two columns: indices of network points 
%           third column     : initial distance
%           forth column     : final (target) distance
%           fifth column     : optional standrad deviation of distances
%
% Output:
% fom       figure of merit, quantifies deviation of the fitted model
%           from restraints, maximum likelihood estimated based on mean
%           distance and standard deviation of mean distance for all
%           restraints
% rmsd      root mean square error of fitted structuer w.r.t. restraints
% network   the fitted network
% dxmat     matrix with coordinate displacement vectors for all steps
%           
% G. Jeschke, 2010

global model
global target

global ENM_param

f=ENM_param.fix_force; % (relative) force constant for restraining local structure
f2=ENM_param.sec_force; % (relative) force constant for restraining local structure
                        % for second neighbor
tol_dist=ENM_param.tol_dist; % minimum tolerance in distance restraint fitting to avoid overfitting

distortion=0.75; % threshold of coordinate change (Å) in relaxation at which fit is assumed to become unstable 

modes=7:6+ENM_param.fit_basis;

maxtime=9600; % maximum time in seconds

[m,n]=size(network0); % m number of residues

dxmat=zeros(m,300);

% set up table of internal restraints for stabilizing local structure
internal=zeros(2*m,4); % initialize table of internal restraints
poi=0;
CaCa_flag=ENM_param.fix_local>=1;
sec_flag=abs(ENM_param.fix_local)>=2;
for k=1:m-2,
    cind0=model.coarse(model.current_structure).indices(k,:);
    cind1=model.coarse(model.current_structure).indices(k+1,:);
    cind2=model.coarse(model.current_structure).indices(k+2,:);
    if CaCa_flag && sum(abs(cind1(1:3)-cind0(1:3)))==0 && cind1(4)-cind0(4)==1, % consecutive residues
        poi=poi+1;
        internal(poi,1)=k;
        internal(poi,2)=k+1;
        internal(poi,3)=norm(network0(k,:)-network0(k+1,:));
        internal(poi,4)=f;
    end;
    if sec_flag && sum(abs(cind2(1:3)-cind0(1:3)))==0 && cind2(4)-cind0(4)==2, % consecutive residues
        poi=poi+1;
        internal(poi,1)=k;
        internal(poi,2)=k+2;
        internal(poi,3)=norm(network0(k,:)-network0(k+2,:));
        internal(poi,4)=f2;
    end;
end;
if poi>0,
    internal=internal(1:poi,:);
else
    internal=[];
end;
mi=poi;
disp(sprintf('There are %i local restraints',poi));

[mr,nr]=size(restr); % m number of residues

if nr>4,
    weights=restr(:,5);
else
    weights=ones(mr,1);
end;

fom0=0;
for k=1:mr, 
    vec=network0(restr(k,2),:)-network0(restr(k,1),:);
    rp=norm(vec);
    fom0=fom0+weights(k)*(restr(k,4)-rp)^2/mr;    
end;
for k=1:mi, 
    vec=network0(internal(k,2),:)-network0(internal(k,1),:);
    rp=norm(vec);
    fom0=fom0+internal(k,4)*(internal(k,3)-rp)^2/mi;    
end;
fom0=sqrt(fom0);
yf=fom0;
rmsd=0;
for k=1:mr, % set up unit vector between restraint residue pairs
    vec=network0(restr(k,2),:)-network0(restr(k,1),:);
    rp=norm(vec);
    rmsd=rmsd+(restr(k,4)-rp)^2;    
end;
rmsd=sqrt(rmsd/mr);

x=0;
y=rmsd;

figure(1); clf;
title('Distance restraint r.m.s.d./ Figure of merit');

h=line(x,y,'Color','k');
hf=line(x,yf,'Color','r');
axis([0,100,0,1.1*y]);

rmsd_struct=rmsd_superimpose(target,network0);

ys=rmsd_struct;

figure(2); clf;
title('Structure r.m.s.d.');

hs=line(x,ys,'Color','b');
axis([0,100,0,1.1*ys]);

numit=0;
runtime=0;
numsat=5;
satflag=false;
dfmax=0;
rmsd=1e6;
best_fit=1e6;
best_it=0;
l2d=0;
tic;
while ~satflag && runtime<maxtime && l2d<0.5 && (rmsd>tol_dist || l2d<0.01),
    [fom,dx,l2d,network]=linear_regression(modes,network0,restr,internal);
    bas=3*numit;
    dxmat(:,bas+1:bas+3)=dx;
    yf=[yf fom];
    dfom=fom0-fom;
    numit=numit+1;
    [rmsd,maxstep,network]=relax_structure(network,internal,numit);
    disp(sprintf('Local contribution to E: %6.3f%%',100*l2d));
    if maxstep>distortion || dfom<0,
       network=network0;
       numit=numit-1;
       break;
    end;
    if dfom<0.05*dfmax, 
        satflag=true;
        disp(sprintf('Saturation detected at step %i',numit));
    end;
    if dfom>dfmax, dfmax=dfom; end;
    fom0=fom;
    network0=network;
    [mr,nr]=size(restr);
    rmsd=0;
    for k=1:mr, % set up unit vector between restraint residue pairs
        vec=network(restr(k,2),:)-network(restr(k,1),:);
        rp=norm(vec);
        rmsd=rmsd+(restr(k,4)-rp)^2;    
    end;
    rmsd=sqrt(rmsd/mr);
    x=[x numit];
    y=[y rmsd];
    set(h,'XData',x,'YData',y);
    rmsd_struct=rmsd_superimpose(target,network0);
    if rmsd_struct<best_fit,
        best_fit=rmsd_struct;
        best_it=numit;
        best_dist=rmsd;
        best_net=network0;
        best_l2d=l2d;
    end;
    ys=[ys rmsd_struct];
    set(hs,'XData',x,'YData',ys);
    set(hf,'XData',x,'YData',yf);
    drawnow;
    runtime=toc;
end;
dxmat=dxmat(:,1:3*numit);
[m,n]=size(dx);
rmsd=0;
[mr,nr]=size(restr);
for k=1:mr, % set up unit vector between restraint residue pairs
    vec=network(restr(k,2),:)-network(restr(k,1),:);
    rp=norm(vec);
    rmsd=rmsd+(restr(k,4)-rp)^2;    
end;
rmsd=sqrt(rmsd/mr);

ys=ys(1:numit);
y=y(1:numit);
figure(4); clf;
plot(y,ys,'k');
title('Structure r.m.s.d. vs. distance r.m.s.d.');
set(gca,'FontSize',14);
axis([0,1.05*max(y),0,1.05*max(ys)]);

figure(1);
axis([0,numit,0,1.05*max(yf)]);
set(gca,'FontSize',14);

figure(2);
axis([0,numit,0,1.05*max(ys)]);
set(gca,'FontSize',14);

disp(sprintf('Fitting took %i iterations. Best fit at iteration %i with local contribution of %5.2f%%',numit,best_it,100*best_l2d));
disp(sprintf('Final fit %4.2f Å with d.r.m.s.d. %4.2f Å.',rmsd_struct,rmsd));
disp(sprintf('Local contribution %5.2f%%',100*l2d));
disp(sprintf('Fitting required %5.1f s',runtime));

% network=best_net;

function [fom,dx,l2d,network]=linear_regression(modes,network,restr,internal)
% sets up and solves linear regression problem for displacement vector
% see Eq. (5), (7), (9), (10) of Zheng/Brooks

global u
global lambda

Hessian=setup_ANM_bonded(network);
[u,d]=eig(Hessian);
clear Hessian
[m,n]=size(d);
lambda=zeros(1,m);
for k=1:m,
    lambda(k)=d(k,k);
end;
clear d

smallstep=0.01; % defines an incremental step within the linear regime

network0=network;
[m,n]=size(network);

if isempty(internal),
    mi=0;
else
    [mi,ni]=size(internal);
end;

[mr,nr]=size(restr);
if nr>4,
    weights=restr(:,5);
else
    weights=ones(mr,1);
end;

nij=zeros(mr+mi,3);
dR=zeros(mr+mi,1);
for k=1:mr, % set up unit vector between restraint residue pairs
    vec=network(restr(k,2),:)-network(restr(k,1),:);
    rp=norm(vec);
    nij(k,:)=vec/rp;
    dR(k)=sqrt(weights(k)/mr)*(restr(k,4)-rp);
end;
if mi>0,
    for k=1:mi, % set up unit vector between local (internal) restraint residue pairs
        vec=network(internal(k,2),:)-network(internal(k,1),:);
        rp=norm(vec);
        nij(k+mr,:)=vec/rp;
        dR(k+mr)=sqrt(internal(k,4)/mi)*(internal(k,3)-rp);
    end;
end;

nm=length(modes);
F=zeros(mr,nm);
for j=1:nm,
    evec=u(:,modes(j));
    mode=reshape(evec,3,m);
    wm=lambda(modes(j));
    for k=1:mr,
        nu=mode(:,restr(k,1))-mode(:,restr(k,2));
        F(k,j)=sqrt(weights(k)/mr)*sum(-nij(k,:).*nu'/wm);
    end;
    if mi>0,
        for k=1:mi,
            nu=mode(:,internal(k,1))-mode(:,internal(k,2));
            F(k+mr,j)=sqrt(internal(k,4)/mi)*sum(-nij(k,:).*nu'/wm);
        end;
    end;
end;

A=F'*F;
y=F'*dR;
% dA=linsolve(A,y);
dA=A'\y;
dx=zeros(m,3);
for j=1:nm,
    evec=u(:,modes(j));
    mode=reshape(evec,3,m);
    mode=mode';
    wm=lambda(modes(j));
    dx=dx+dA(j)*mode/wm;
end;

E=0;
E0=0;
lE=0;
dE=0;
network=network0+smallstep*dx;
for k=1:mr, % set up unit vector between restraint residue pairs
    vec0=network0(restr(k,2),:)-network0(restr(k,1),:);
    rp0=norm(vec0);
    E0=E0+weights(k)*(restr(k,4)-rp0)^2/mr;
    vec=network(restr(k,2),:)-network(restr(k,1),:);
    rp=norm(vec);
    dE=dE+weights(k)*(restr(k,4)-rp)^2/mr;
    E=E+weights(k)*(restr(k,4)-rp)^2/mr;    
end;
for k=1:mi, % set up unit vector between restraint residue pairs
    vec0=network0(internal(k,2),:)-network0(internal(k,1),:);
    rp0=norm(vec0);
    E0=E0+internal(k,4)*(internal(k,3)-rp0)^2/mi;
    vec=network(internal(k,2),:)-network(internal(k,1),:);
    rp=norm(vec);
    lE=lE+internal(k,4)*(internal(k,3)-rp)^2/mi;
    E=E+internal(k,4)*(internal(k,3)-rp)^2/mi;    
end;

l2d=lE/(lE+dE);
fE=(E0-E)/E0;
sc=0.1*smallstep/fE; % error function should decrease by about 10%
dx=sc*dx;
steps=sqrt(sum(dx.^2,2));
sc2=0.5/max(steps);
if sc2<1,
    dx=sc2*dx;
    disp(sprintf('Reduced step size by factor %5.3f',sc2));
end;

% disp(sprintf('Original m.s.e. : %6.4f',E0));
network=network0+dx;
% compute actual error function
fom=0;
for k=1:mr, 
    vec=network(restr(k,2),:)-network(restr(k,1),:);
    rp=norm(vec);
    fom=fom+weights(k)*(restr(k,4)-rp)^2/mr;    
end;
for k=1:mi, 
    vec=network(internal(k,2),:)-network(internal(k,1),:);
    rp=norm(vec);
    fom=fom+internal(k,4)*(internal(k,3)-rp)^2/mi;    
end;
% disp(sprintf('New      m.s.e. : %6.4f',fom));
fom=sqrt(fom);

function [rmsd,maxstep,network]=relax_structure(network,internal,iteration)
% sets up and solves linear regression problem for local structure
% relaxation
% see Eq. (5), (7), (9), (10) of Zheng/Brooks and remarks at the end of p.
% 4329

global u
global lambda
global ENM_param


tol_local=iteration*0.5/50; % targeted maximum r.m.s.d. (Å) after 50 iterations
stepmax=0.2; % no coordinate can change by more than 0.5 Å in relaxation

[m,n]=size(network);
modes=7+ENM_param.fit_basis:m; % all modes not used for fitting distance restraints

if isempty(internal),
    rmsd=0;
    dx=zeros(m,3);
    sc=0;
    return
else
    [mi,ni]=size(internal);
end;

sumdx=zeros(m,3);

E0=0;
for k=1:mi, % set up unit vector between restraint residue pairs
    vec0=network(internal(k,2),:)-network(internal(k,1),:);
    rp0=norm(vec0);
    E0=E0+internal(k,4)*(internal(k,3)-rp0)^2/mi;
end;
rmsd0=sqrt(E0);
rmsd=rmsd0;

if rmsd<=tol_local,
    dx=zeros(m,3);
end;

maxit=25;
it=0;
maxstep=0;
ms0=0;
while rmsd>tol_local && it<maxit,
    it=it+1;
    nij=zeros(mi,3);
    dR=zeros(mi,1);
    for k=1:mi, % set up unit vector between local (internal) restraint residue pairs
        vec=network(internal(k,2),:)-network(internal(k,1),:);
        rp=norm(vec);
        nij(k,:)=vec/rp;
        dR(k)=sqrt(internal(k,4)/mi)*(internal(k,3)-rp);
    end;

    nm=length(modes);
    F=zeros(mi,nm);
    for j=1:nm,
        evec=u(:,modes(j));
        mode=reshape(evec,3,m);
        wm=lambda(modes(j));
        for k=1:mi,
            nu=mode(:,internal(k,1))-mode(:,internal(k,2));
            F(k,j)=sqrt(internal(k,4)/mi)*sum(-nij(k,:).*nu'/wm);
        end;
    end;

    A=F'*F;
    y=F'*dR;
    dA=(A')\y;
    dx=zeros(m,3);
    for j=1:nm,
        evec=u(:,modes(j));
        mode=reshape(evec,3,m);
        mode=mode';
        wm=lambda(modes(j));
        dx=dx+dA(j)*mode/wm;
    end;

    E=0;
    E0=0;
    network0=network;
    network=network+0.05*dx;
    sumdx=sumdx+0.05*dx;
    maxstep=max(sqrt(sum(sumdx.^2,2)));
    for k=1:mi, % set up unit vector between restraint residue pairs
        vec0=network0(internal(k,2),:)-network0(internal(k,1),:);
        rp0=norm(vec0);
        E0=E0+internal(k,4)*(internal(k,3)-rp0)^2/mi;
        vec=network(internal(k,2),:)-network(internal(k,1),:);
        rp=norm(vec);
        E=E+internal(k,4)*(internal(k,3)-rp)^2/mi;    
    end;
    rmsd=sqrt(E);
    if E>=E0 || maxstep>stepmax,
        network=network0;
        rmsd=sqrt(E0);
        break;
    end;
    ms0=maxstep;
end;
maxstep=ms0;
disp(sprintf('Original rmsd : %5.3f at iteration %i',rmsd0,iteration));
disp(sprintf('Targeted rmsd : %5.3f and maximum step: %4.2f',tol_local,maxstep));
disp(sprintf('Achieved rmsd : %5.3f with %i iterations',rmsd,it));
