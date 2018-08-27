function [fom,rmsd,network,dxmat,best] = fit_displacements_by_ANM(network0,displacements,test_mode,correspondence,target)
% function [fom,rmsd,network,dxmat,best] = fit_by_ANM_transition(network0,displacements,test_mode,correspondence,target,maxbas)
%
% fits a network model to DEER and optional direct Calpha-Calpha distance 
% restraints based on equipartitioning of energy over the modes, mode signs
% are derived from direction of the improvement in distance constraint
% r.m.s.d.
%
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
% displacements     Calpha displacement constraint variable, n-by-5 array
%                   1st column     : indices of network points 
%                   2nd-4th column : displacement vectors (Å) 
%                   5th column     : weights of network points
% test_mode optional flag, if true, a target Calpha coordinate array must 
%           be provided as argument target and a correspondence table to 
%           network0, information is then displayed on fit
%           progress and final fit quality w.r.t. this target structure,
%           defaults to false (no target structure known)
% correspondence    correspondebce table between template and target residue
%                   coordinates (only in test mode)
% target            target Calpha coordinates (only in test mode)         
%
% Output:
% fom       figure of merit, quantifies deviation of the fitted model
%           from restraints, maximum likelihood estimated based on mean
%           distance and standard deviation of mean distance for all
%           restraints
% rmsd      root mean square error of fitted structuer w.r.t. restraints
% network   the fitted network
% dxmat     matrix with coordinate displacement vectors for all steps
% best      information on best approach to target structure (in test_mode
%           only)
%           best.it     iteration cycle
%           best.rmsd   structure rmsd
%           best.drmsd  distance rmsd
%
% NOTE that further parameters are defined in global variable ENM_param,
% which is initialized in initialize_MMM.m
%           
% G. Jeschke, 2012

global model

global ENM_param

if nargin<3,
    test_mode=false;
end;
maximum_likelihood=true;

f=ENM_param.fix_force; % (relative) force constant for restraining local structure
f2=ENM_param.sec_force; % (relative) force constant for restraining local structure
                        % for second neighbor

maxtime=12*3600; % maximum time in seconds

[m,n]=size(network0); % m number of residues

dxmat=zeros(m,150);

[md,nd]=size(displacements);
basis=md*nd;
if basis>m,
    basis=m;
end;

indices=displacements(:,1);
weights=displacements(:,5);

active0=network0(indices,:);
vecs=displacements(:,2:4);
active_target=active0+vecs;

vdiff=sum(vecs.^2,2);
fom0=sqrt(sum(vdiff.*weights)/sum(weights));
yf=fom0;
x=0;
best.initial_DEER=0;

figure(1); clf;
title('Constraint r.m.s.d./ Figure of merit');

hf=line(x,yf,'Color','r');
axis([0,ENM_param.cycles,0,1.1*yf]);


if test_mode,
    rmsd_struct=rmsd_superimpose(target(correspondence(2,:),:),network0(correspondence(1,:),:));

    ys=rmsd_struct;
    fprintf(1,'Initial structure r.m.s.d. is %5.2f Å.\n',ys);
    figure(2); clf;
    title('Structure r.m.s.d.');

    hs=line(x,ys,'Color','b');
    axis([0,100,0,1.1*ys]);
end;

numit=0;
runtime=0;
dfmax=0;
ddmax=0;
best_fit=1e6;
best_it=0;
l2d=0;
struct_step=zeros(1,ENM_param.cycles);
change=zeros(1,ENM_param.cycles);
scaling=zeros(1,ENM_param.cycles);
dx_step=zeros(1,ENM_param.cycles);
fom_vec=zeros(1,ENM_param.cycles);
best_network=network0;

network00=network0;
no_stop=10;
ldfom_vec=zeros(1,200);
no_sat=200;

[ma,na]=size(network0);

fit_limit=0.2;

limit_found=false;
tic;
while (no_stop || no_sat) && runtime<maxtime, %(~satflag && runtime<maxtime && rmsd>fit_limit && l2d<0.5 && (rmsd>tol_dist || l2d<0.01)),
    if no_stop>0,
        no_stop=no_stop-1;
    else
        no_stop=0;
    end;
    if no_sat>0,
        no_sat=no_sat-1;
    else
        no_sat=0;
    end;
    [fom,dx,network,sc]=transition_step(network0,indices,active_target,weights,basis);
    bas=3*numit;
    dxmat(:,bas+1:bas+3)=dx;
    yf=[yf fom];
    dfom=fom0-fom;
    numit=numit+1;
    ldfom_vec(numit)=log(fom0)-log(fom);
    fom_vec(numit)=fom;
    scaling(numit)=sc;
    networks{numit}=network;
    struct_step(numit)=rmsd_superimpose(network,network0);
    change(numit)=rmsd_superimpose(network,network00);
    if dfom>dfmax, 
        dfmax=dfom; 
    end;
    if dfom<dfmax/100 && no_sat>10,
        no_sat=10;
    end;
    fom0=fom;
    network0=network;
    x=[x numit];
    set(hf,'XData',x,'YData',yf);
    if length(x)>100, nnx=length(x); else nnx=100; end;
    axis([0,nnx,0,1.1*max(yf)]);
    if test_mode,
        rmsd_struct=rmsd_superimpose(target(correspondence(2,:),:),network0(correspondence(1,:),:));
        if rmsd_struct<best_fit,
            best_fit=rmsd_struct;
            best_it=numit;
            best.rmsd=rmsd_struct;
            best.it=numit;
        end;
        ys=[ys rmsd_struct];
        set(hs,'XData',x,'YData',ys);
        axis([0,nnx,0,1.1*max(ys)]);
    end;
    drawnow;
    runtime=toc;
end;

best_guess=numit;

if runtime>=maxtime,
    add_msg_board('Fit stopped as maximum time was exceeded.');
end;
dxmat=dxmat(:,1:3*numit);
[m,n]=size(dx);


der2=diff(yf,2);
der2(yf(2:end-1)>1.2*min(yf))=0;
der2s=der2;
der2s(1)=(2*der2(1)+der2(2))/3;
der2s(end)=(2*der2(end)+der2(end-1))/3;
for k=2:length(der2)-1,
    der2s(k)=(der2(k-1)+2*der2(k)+der2(k+1))/4;
end;
[ma,pois]=max(der2s(1:end-1));
converges=pois+1;


figure(3); clf;
plot(der2,'k');
hold on;
plot(der2s(1:end-1),'g');

network=networks{converges};

active=network(indices,:);
vdiff=active-active_target;
vdiff=sum(vdiff.^2,2);
fom=sqrt(sum(vdiff.*weights)/sum(weights));
rmsd=sqrt(sum(vdiff)/length(vdiff));
add_msg_board(sprintf('Final figure of merit         : %4.2f Å',fom));
add_msg_board(sprintf('Final r.m.s.d. of active atoms: %4.2f Å',rmsd));

if test_mode,
    ys=ys(1:numit);
    y=y(1:numit);
    axis([0,1.05*max(y),0,1.05*max(ys)]);    
    add_msg_board(sprintf('Final structure r.m.s.d. %5.2f\n',rmsd_struct));
end;


figure(1);
axis([0,numit,0,1.05*max(yf)]);
set(gca,'FontSize',14);
hold on;
plot([best_guess best_guess],[0,max(yf)],'r:');
plot(converges-1,yf(converges),'ko');

if test_mode,
    figure(2);
    axis([0,numit,0,1.05*max(ys)]);
    set(gca,'FontSize',14);
    hold on;
    plot([best_guess best_guess],[0,max(ys)],'r:');
end;

if test_mode,
    network=networks{converges};
    rmsd_struct=rmsd_superimpose(target(correspondence(2,:),:),network(correspondence(1,:),:));
    fprintf(1,'Best r.m.s.d. during whole fit was %5.2f at iteration %i.\n',best_fit,best_it);
    fprintf(1,'Converged r.m.s.d. was %5.2f at iteration %i.\n',ys(converges),converges);
    fprintf(1,'Selected structure has r.m.s.d. of %5.2f at iteration %i.\n',rmsd_struct,converges);
    add_msg_board(sprintf('Final fit %4.2f Å with d.r.m.s.d. %4.2f Å.',rmsd_struct,rmsd));
else
    add_msg_board(sprintf('Fitting took %i iterations.',numit));
end;
add_msg_board(sprintf('Fitting required %5.1f s',runtime));


function [fom,dx,network,sc]=transition_step(network,indices,active_target,weights,basis)
% computes one transition step

global model

if isfield(model.ANM(model.current_structure),'u'),
    model.ANM(model.current_structure).u=[];
end;
Hessian=setup_ANM_bonded(network);
clear u
[u,d]=eig(Hessian);
clear Hessian
lambda=diag(d);
clear d
ut=u(:,7:basis+6);
lambdat=lambda(7:basis+6);
nt=length(lambdat);

smallstep=0.2; % defines an incremental step within the linear regime

network0=network;
[m,n]=size(network);

active=network(indices,:);
dR=active_target-active;
[md,nd]=size(dR);

phases=zeros(basis,1);

eweights=repmat(weights',3,1);
eweights=reshape(sqrt(eweights),3*length(weights),1);
dR=reshape(dR,nd*md,1);
dR2=dR.*eweights;
dR2=dR2/norm(dR2);

for j=1:nt,
    evec=ut(:,j);
    mode=reshape(evec,3,m);
    mode=mode.';
    smode=mode(indices,:);
    [ms,ns]=size(smode);
    smode=reshape(smode,ns*ms,1);
    smode=smode.*eweights;
    smode=smode/norm(smode);
    phases(j)=sum(smode.*dR2);
end;

coeff=sqrt(ones(size(lambdat))./lambdat);
step=ut(:,1:basis)*(phases(1:basis).*coeff(1:basis));
step=reshape(step,3,m);
step=step';
per_residue=sqrt(sum(step.^2,2));
sc=smallstep/max(per_residue);
dx=sc*step;
network=network0+dx;

active=network(indices,:);
vdiff=active_target-active;
vdiff=sum(vdiff.^2,2);
fom=sqrt(sum(vdiff.*weights)/sum(weights));

