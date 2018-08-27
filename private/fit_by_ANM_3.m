function [fom,rmsd,network,dxmat,DEER] = fit_by_ANM_3(network0,DEER,direct,test_mode)
% function [fom,rmsd,network,dxmat,DEER] = fit_by_ANM_3(network0,DEER,direct,test_mode)
%
% fits a network model to DEER and optional direct Calpha-Calpha distance 
% restraints based on the algorithm from
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
% DEER      DEER restraint variable 1xn structure array for n constraints
%           all in  units of nm
%           .r          mean distances
%           .sigr       uncertainties of mean distances
%           .indices    2x4 arrays of MMM indices for residue pair
%           .xyz1       Cartesian coordinates of mean position first label
%           .rmsd1      position r.m.s.d. of first label
%           .xyz2       Cartesian coordinates of mean position second label
%           .rmsd2      position r.m.s.d. of second label
%           .l1         handle to graphics of first label in 3D model (in
%                       fit_from_template)
%           .rl1        handle to graphics of first label tether in 3D 
%                       model (in fit_from_template)
%           .l2         handle to graphics of second label in 3D model (in
%                       fit_from_template)
%           .rl2        handle to graphics of second label tether in 3D 
%                       model (in fit_from_template)
%           .ll         handle to graphics of distance vector in 3D model
%                       (in fit_from_template)
%           .res1       number of first residue in Calpha network model
%           .res2       number of second residue in Calpha network model
% direct    optional Calpha-Calpha restraint variable, n-by-5 array
%           all in units of nm
%           first two columns: indices of network points 
%           third column     : initial distance (nm)
%           forth column     : final (target) distance
%           fifth column     : standard deviation of distance, strength of
%                              restraint
% test_mode optional flag, if true, a target structure must be defined in
%           global variable 'target', information is then displayed on fit
%           progress and final fit quality w.r.t. this target structure,
%           defaults to false (no target structure known)
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
% diagnostic information output:
% matrices with size [number of iterations, number of modes in fit basis]
% mode_spectrum     distribution of structural changes over the first slow
%                   modes in the template structure
% mode_similarity   similarity of modes in later iteration steps to modes
%                   of the template structure (dot product)
%           
% G. Jeschke, 2010

global model
global target

global ENM_param

global u

if nargin<4,
    test_mode=false;
end;
maximum_likelihood=true;
saturation_threshold=0.01;

f=ENM_param.fix_force; % (relative) force constant for restraining local structure
f2=ENM_param.sec_force; % (relative) force constant for restraining local structure
                        % for second neighbor
tol_dist=ENM_param.tol_dist; % minimum tolerance in distance restraint fitting to avoid overfitting

distortion=0.75; % threshold of coordinate change (Å) in relaxation at which fit is assumed to become unstable 

modes=7:6+ENM_param.fit_basis;

maxtime=12*3600; % maximum time in seconds

[m,n]=size(network0); % m number of residues

dxmat=zeros(m,150);


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
add_msg_board(sprintf('There are %i local restraints',poi));

mr=length(DEER); % m number of residues

if maximum_likelihood && mr>0,
    weights=zeros(mr,1);
    for k=1:mr,
        weights(k)=1/DEER(k).sigr;
    end;
else
    weights=ones(mr,1);
end;
if ~isempty(DEER),
    weights=weights.^2;
    sc=mr/sum(weights);
    weights=sc*weights; % renormalization to conform to Eq. (4) Zheng/Brooks
else
    weights=0;
    sc=1;
end;

if nargin>2 && ~isempty(direct)
    [md,nd]=size(direct); % m number of residues
    if maximum_likelihood,
        dir_weights=zeros(md,1);
        for k=1:md,
            dir_weights(k)=1/direct(k,5);
        end;
    else
        dir_weights=ones(md,1);
    end;
    dir_weights=dir_weights.^2;
    weights=weights/sc;
    sc2=(md+mr)/(sum(dir_weights)+sum(weights));
    dir_weights=sc2*dir_weights; % renormalization to conform to Eq. (4) Zheng/Brooks
    weights=sc2*weights; % renormalization to conform to Eq. (4) Zheng/Brooks
else
    md=0;
    dir_weights=[];
end;

fom0=0;
if ~isempty(DEER),
    for k=1:mr, 
        vec=DEER(k).xyz2-DEER(k).xyz1;
        rp=norm(vec);
        fom0=fom0+weights(k)*(10*DEER(k).r-rp)^2/(mr+md);    
    end;
end;
if md>0, % direct constraints exist
    for k=1:md,
        vec=network0(direct(k,2),:)-network0(direct(k,1),:);
        rp=norm(vec);
        fom0=fom0+dir_weights(k)*(10*direct(k,4)-rp)^2/(mr+md);    
    end;
end;
for k=1:mi, 
    vec=network0(internal(k,2),:)-network0(internal(k,1),:);
    rp=norm(vec);
    fom0=fom0+internal(k,4)*(internal(k,3)-rp)^2/mi;    
end;
fom0=sqrt(fom0);
yf=fom0;
rmsd=0;
fit_limit=0;
x=0;
if ~isempty(DEER),
    for k=1:mr, % set up unit vector between restraint residue pairs
        vec=DEER(k).xyz2-DEER(k).xyz1;
        rp=norm(vec);
        rmsd=rmsd+(10*DEER(k).r-rp)^2;    
        fit_limit=fit_limit+(10*DEER(k).sigr)^2;
    end;
    fit_limit=sqrt(fit_limit/mr^2);
    rmsd=sqrt(rmsd/mr);
    rmsd0=rmsd;
    add_msg_board(sprintf('Initial DEER distance constraint r.m.s.d. %6.2f',rmsd));
else
    rmsd=0.001;
end;
y=rmsd;
if fit_limit>2,
    fit_limit=2;
end;

figure(1); clf;
title('DEER constraint r.m.s.d./ Figure of merit');

if md>0,
    rmsd_dir=0;
    for k=1:md, % set up unit vector between restraint residue pairs
        vec=network0(direct(k,2),:)-network0(direct(k,1),:);
        rp=norm(vec);
        rmsd_dir=rmsd_dir+(10*direct(k,4)-rp)^2;    
    end;
    rmsd_dir=sqrt(rmsd_dir/md);
    add_msg_board(sprintf('Initial direct distance constraint r.m.s.d. %6.2f',rmsd_dir));
    ydir=rmsd_dir;
    rmsd_dir0=rmsd_dir;
    hdir=line(x,ydir,'Color','b');
else
    ydir=-1e6;
end;

if ~isempty(DEER),
    h=line(x,y,'Color','k');
end;
hf=line(x,yf,'Color','r');
axis([0,100,0,1.1*max([y,yf,ydir])]);


if test_mode,
    rmsd_struct=rmsd_superimpose(target,network0);

    ys=rmsd_struct;
    fprintf(1,'Initial structure r.m.s.d. is %5.2f Å.\n',ys);
    figure(2); clf;
    title('Structure r.m.s.d.');

    hs=line(x,ys,'Color','b');
    axis([0,100,0,1.1*ys]);
end;

numit=0;
runtime=0;
satflag=false;
dfmax=0;
ddmax=0;
rmsd=1e6;
best_fit=1e6;
best_it=0;
l2d=0;
struct_step=zeros(1,100);
scaling=zeros(1,100);
dx_step=zeros(1,100);
mode_spectrum=zeros(100,length(modes));
mode_similarity=zeros(100,length(modes));
best_network=network0;
stop_rmsd=1e6;
stop_it=0;
maxsc=0;
tic;
while ~satflag && runtime<maxtime && rmsd>fit_limit && l2d<0.5 && (rmsd>tol_dist || l2d<0.01),
    [fom,dx,l2d,network,DEER2,maxstep,mode_weights,sc]=linear_regression(modes,network0,DEER,internal,weights,direct,dir_weights);
%     if numit==0,
%         disp(size(u));
%         disp(size(modes));
%         disp(modes);
%         mode_bas=u(:,modes);
%         mode_spectrum(1,:)=mode_weights;
%         mode_similarity(1,:)=ones(size(mode_weights));
%     else
%         % mode following
%         simmat=zeros(length(modes));
%         for k=1:length(modes),
%             for kk=1:length(modes),
%                 simmat(k,kk)=sum(mode_bas(:,k).*u(:,modes(kk)));
%             end;
%         end;
%         for k=1:length(modes),
%             % determine mode most similar to mode in template structure
%             [cp,poi1]=max(simmat(k,:));
%             cp2=max(simmat(:,poi1));
%             if cp>=cp2,
%                 mode_spectrum(numit+1,k)=mode_weights(poi1);
%                 mode_similarity(numit+1,k)=cp;
%             end;
%         end;
%     end;
    bas=3*numit;
    dxmat(:,bas+1:bas+3)=dx;
    yf=[yf fom];
    dfom=fom0-fom;
    numit=numit+1;
    scaling(numit)=sc;
    dx_step(numit)=maxstep;
    % [rmsd,maxstep,network]=relax_structure(network,internal,numit);
    struct_step(numit)=rmsd_superimpose(network,network0);
    if sc>maxsc,
        maxsc=sc;
        best_network=network;
        if test_mode,
            stop_rmsd=rmsd_superimpose(target,network);
            stop_it=numit;
        end;
    end;
%    disp(sprintf('Local contribution to E: %6.3f%%',100*l2d));
%     if maxstep>distortion % || dfom<0,
%        network=network0;
%        numit=numit-1;
%        break;
%     end;
    DEER=DEER2;
    if dfom<saturation_threshold*dfmax,
        satflag=true;
        add_msg_board(sprintf('Total energy function saturation detected at step %i',numit));
    end;
    if dfom>dfmax, 
        dfmax=dfom; 
    end;
    fom0=fom;
    network0=network;
    x=[x numit];
    if ~isempty(DEER),
        mr=length(DEER);
        rmsd=0;
        for k=1:mr, % set up unit vector between restraint residue pairs
            vec=DEER(k).xyz2-DEER(k).xyz1;
            rp=norm(vec);
            rmsd=rmsd+(10*DEER(k).r-rp)^2;    
        end;
        rmsd=sqrt(rmsd/mr);
        ddis=rmsd0-rmsd;
        if ddis>ddmax, ddmax=ddis; end;
        rmsd0=rmsd;
        if ddis<0.05*ddmax, 
            satflag=true;
            add_msg_board(sprintf('DEER constraint saturation detected at step %i',numit));
        end;
        y=[y rmsd];
        set(h,'XData',x,'YData',y);
    end;
    set(hf,'XData',x,'YData',yf);
    if md>0,
        rmsd_dir=0;
        for k=1:md, % set up unit vector between restraint residue pairs
            vec=network(direct(k,2),:)-network(direct(k,1),:);
            rp=norm(vec);
            rmsd_dir=rmsd_dir+(10*direct(k,4)-rp)^2;    
        end;
        rmsd_dir=sqrt(rmsd_dir/md);
        ddisdir=rmsd_dir0-rmsd_dir;
        rmsd_dir0=rmsd_dir;
        ydir=[ydir rmsd_dir];
        set(hdir,'XData',x,'YData',ydir);
    end;
    if length(x)>100, nnx=length(x); else nnx=100; end;
    axis([0,nnx,0,1.1*max([max(y),max(yf),max(ydir)])]);
    if test_mode,
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
        axis([0,nnx,0,1.1*max(ys)]);
    end;
    drawnow;
    runtime=toc;
end;
if satflag,
    add_msg_board('Fit converged according to energy or distance constraint criterion.');
end;
if runtime>=maxtime,
    add_msg_board('Fit stopped as maximum time was exceeded.');
end;
if rmsd<=fit_limit,
    add_msg_board('Fit stopped since distance r.m.s.d. fell below limit determined by specified distance ranges');
end;
if l2d>=0.5,
    add_msg_board('Fit stopped as local structure distortions were too large.');
end;
if (rmsd<=tol_dist && l2d<0.01),
    add_msg_board('Fit stopped as distance r.m.s.d. is reasonably good and local structure distortion becomes significant.');
end;
% mode_spectrum=mode_spectrum(1:numit,:);
% mode_similarity=mode_similarity(1:numit,:);
dxmat=dxmat(:,1:3*numit);
scaling=scaling(1:numit);
[m,n]=size(dx);
rmsd=10000;
if ~isempty(DEER),
    rmsd=0;
    mr=length(DEER);
    for k=1:mr, % set up unit vector between restraint residue pairs
        vec=DEER(k).xyz2-DEER(k).xyz1;
        rp=norm(vec);
        rmsd=rmsd+(10*DEER(k).r-rp)^2;    
    end;
    rmsd=sqrt(rmsd/mr);
end;

if md>0,
    rmsd_dir=0;
    for k=1:md, % set up unit vector between restraint residue pairs
        vec=network0(direct(k,2),:)-network0(direct(k,1),:);
        rp=norm(vec);
        rmsd_dir=rmsd_dir+(10*direct(k,4)-rp)^2;    
    end;
    rmsd_dir=sqrt(rmsd_dir/md);
    add_msg_board(sprintf('Final direct distance constraint r.m.s.d. %6.2f',rmsd_dir));
end;

if test_mode && ~isempty(DEER),
    ys=ys(1:numit);
    y=y(1:numit);
    figure(4); clf;
    plot(y,ys,'k');
    title('Structure r.m.s.d. vs. distance r.m.s.d.');
    set(gca,'FontSize',14);

    axis([0,1.05*max(y),0,1.05*max(ys)]);
    dx_step=dx_step(1:numit);
    ys=ys(1:numit);
    xx=1:numit;
    figure(6); clf;
    plot(xx,dx_step,'k');
    hold on;
    plot(xx,ys,'r');
end;

struct_step=struct_step(1:numit);
figure(5);
clf;
plot(struct_step,'k');
title('R.m.s.d. between consecutive structures');

figure(6);
clf;
plot(scaling,'k');
title('Step scaling factor');
hold on
[masc,best]=max(scaling);
plot([best,best],[0,masc],'g:');
fprintf(1,'Extracting structure at iteration %i.\n',best);
%rmsd=stop_rmsd;
network=best_network;

figure(1);
axis([0,numit,0,1.05*max(yf)]);
set(gca,'FontSize',14);

if test_mode,
    figure(2);
    axis([0,numit,0,1.05*max(ys)]);
    set(gca,'FontSize',14);
end;

if test_mode,
    fprintf(1,'Best r.m.s.d. during whole fit was %5.2f at iteration %i.\n',best_fit,best_it);
    fprintf(1,'Last r.m.s.d. was %5.2f at iteration %i.\n',ys(numit),numit);
    fprintf(1,'Selected structure has r.m.s.d. of %5.2f at iteration %i.\n',stop_rmsd,stop_it);
    add_msg_board(sprintf('Fitting took %i iterations. Best fit at iteration %i with local contribution of %5.2f%%',numit,best_it,100*best_l2d));
    add_msg_board(sprintf('Final fit %4.2f Å with d.r.m.s.d. %4.2f Å.',rmsd_struct,rmsd));
    add_msg_board(sprintf('Local contribution %5.2f%%',100*l2d));
else
    add_msg_board(sprintf('Fitting took %i iterations.',numit));
end;
add_msg_board(sprintf('Fitting required %5.1f s',runtime));
% network=best_net;

function [fom,dx,l2d,network,DEER2,maxstep,mode_weights,sc2]=linear_regression(modes,network,DEER,internal,weights,direct,dir_weights)
% sets up and solves linear regression problem for displacement vector
% see Eq. (5), (7), (9), (10) of Zheng/Brooks

global u
global lambda
global model

transrot1=false;

md=0;
if nargin>5,
    if ~isempty(direct),
        [md,nd]=size(direct);
    end;
end;
if isfield(model.ANM(model.current_structure),'u'),
    model.ANM=rmfield(model.ANM(model.current_structure),'u');
end;
Hessian=setup_ANM_bonded(network);
clear u
[u,d]=eig(Hessian);
clear Hessian
[m,n]=size(d);
lambda=zeros(1,m);
for k=1:m,
    lambda(k)=d(k,k);
end;
clear d

% make modes for rotation/translation of first chain
chains=model.coarse(model.current_structure).chains;
% disp(sprintf('First chain: %i from residue %i to residue %i',chains(1,1),chains(1,2),chains(1,3)));
% disp(sprintf('Second chain: %i from residue %i to residue %i',chains(2,1),chains(2,2),chains(2,3)));
smallnet=network(chains(1,2):chains(1,3),:);
Hessian=setup_ANM_bonded(smallnet);
[us,ds]=eig(Hessian);
clear Hessian
clear ds
uchain1=us(:,1:6);

smallstep=0.01; % defines an incremental step within the linear regime

network0=network;
[m,n]=size(network);
[m1,n1]=size(smallnet);

if isempty(internal),
    mi=0;
else
    [mi,ni]=size(internal);
end;

mr=length(DEER);

nij=zeros(mr+md+mi,3);
dR=zeros(mr+md+mi,1);
if mr>0,
    for k=1:mr, % set up unit vector between restraint residue pairs
        vec=DEER(k).xyz2-DEER(k).xyz1;
        rp=norm(vec);
        nij(k,:)=vec/rp;
        dR(k)=sqrt(weights(k)/(mr+md))*(10*DEER(k).r-rp);
    end;
end;
if md>0,
    for k=1:md, % set up unit vector between restraint residue pairs
        vec=network0(direct(k,2),:)-network0(direct(k,1),:);
        rp=norm(vec);
        nij(k+mr,:)=vec/rp;
        dR(k+mr)=sqrt(dir_weights(k)/(mr+md))*(10*direct(k,4)-rp);
    end;
end;
if mi>0,
    for k=1:mi, % set up unit vector between local (internal) restraint residue pairs
        vec=network(internal(k,2),:)-network(internal(k,1),:);
        rp=norm(vec);
        nij(k+md+mr,:)=vec/rp;
        dR(k+md+mr)=sqrt(internal(k,4)/mi)*(internal(k,3)-rp);
    end;
end;

nm=length(modes);
nc1=0;
if transrot1, nc1=6; end;
F=zeros(mr+md+mi,nm+nc1);
for j=1:nm,
    evec=u(:,modes(j));
    mode=reshape(evec,3,m);
    wm=lambda(modes(j));
    if mr>0,
        for k=1:mr,
            nu=mode(:,DEER(k).res1)-mode(:,DEER(k).res2);
            F(k,j)=sqrt(weights(k)/(mr+md))*sum(-nij(k,:).*nu'/wm);
        end;
    end;
    if md>0,
        for k=1:md,
            nu=mode(:,direct(k,1))-mode(:,direct(k,2));
            F(k+mr,j)=sqrt(dir_weights(k)/(mr+md))*sum(-nij(k+mr,:).*nu'/wm);
        end;
    end;
    if mi>0,
        for k=1:mi,
            nu=mode(:,internal(k,1))-mode(:,internal(k,2));
            F(k+mr+md,j)=1e5*sqrt(internal(k,4)/mi)*sum(-nij(k+mr+md,:).*nu'/wm);
        end;
    end;
end;
if transrot1,
    for j=1:nc1,
        evec=uchain1(:,j);
        mode0=reshape(evec,3,m1);
        mode=zeros(3,m);
        mode(:,1:m1)=mode0;
        wm=lambda(modes(1))/10;
        if mr>0,
            for k=1:mr,
                nu=mode(:,DEER(k).res1)-mode(:,DEER(k).res2);
                F(k,nm+j)=sqrt(weights(k)/(mr+md))*sum(-nij(k,:).*nu'/wm);
            end;
        end;
        if md>0,
            for k=1:md,
                nu=mode(:,direct(k,1))-mode(:,direct(k,2));
                F(k+mr,nm+j)=sqrt(dir_weights(k)/(mr+md))*sum(-nij(k+mr,:).*nu'/wm);
            end;
        end;
        if mi>0,
            for k=1:mi,
                nu=mode(:,internal(k,1))-mode(:,internal(k,2));
                F(k+mr+md,nm+j)=1e5*sqrt(internal(k,4)/mi)*sum(-nij(k+mr+md,:).*nu'/wm);
            end;
        end;
    end;
end;

A=F'*F;
y=F'*dR;
% dA=linsolve(A,y);
dA=A'\y;
dx=zeros(m,3);
mode_weights=dA(1:nm);
for j=1:nm,
    evec=u(:,modes(j));
    mode=reshape(evec,3,m);
    mode=mode';
    wm=lambda(modes(j));
    dx=dx+dA(j)*mode/wm;
end;
if transrot1,
    for j=1:nc1,
        evec=uchain1(:,j);
        mode0=reshape(evec,3,m1);
        mode=zeros(3,m);
        mode(:,1:m1)=mode0;
        mode=mode';
        wm=lambda(modes(1))/5;
        dx=dx+dA(nm+j)*mode/wm;
    end;
end;

E=0;
E0=0;
lE=0;
dE=0;
network=network0+smallstep*dx;
if isempty(DEER),
    DEER2=DEER;
else
    DEER2=update_labels(DEER,network0,network);
    for k=1:mr, % set up unit vector between restraint residue pairs
        vec0=DEER(k).xyz2-DEER(k).xyz1;
        rp0=norm(vec0);
        E0=E0+weights(k)*(10*DEER(k).r-rp0)^2/(mr+md);
        vec=DEER2(k).xyz2-DEER2(k).xyz1;
        rp=norm(vec);
        dE=dE+weights(k)*(10*DEER(k).r-rp)^2/(mr+md);
        E=E+weights(k)*(10*DEER(k).r-rp)^2/(mr+md);    
    end;
end;
if md>0,
    for k=1:md, % set up unit vector between restraint residue pairs
        vec0=network0(direct(k,2),:)-network0(direct(k,1),:);
        rp0=norm(vec0);
        E0=E0+dir_weights(k)*(10*direct(k,4)-rp0)^2/(mr+md);
        vec=network(direct(k,2),:)-network(direct(k,1),:);
        rp=norm(vec);
        dE=dE+dir_weights(k)*(10*direct(k,4)-rp)^2/(mr+md);
        E=E+dir_weights(k)*(10*direct(k,4)-rp)^2/(mr+md);    
    end;
end;
if mi>0,
    for k=1:mi, % set up unit vector between restraint residue pairs
        vec0=network0(internal(k,2),:)-network0(internal(k,1),:);
        rp0=norm(vec0);
        E0=E0+internal(k,4)*(internal(k,3)-rp0)^2/mi;
        vec=network(internal(k,2),:)-network(internal(k,1),:);
        rp=norm(vec);
        lE=lE+internal(k,4)*(internal(k,3)-rp)^2/mi;
        E=E+internal(k,4)*(internal(k,3)-rp)^2/mi;    
    end;
end;

l2d=lE/(lE+dE);
fE=(E0-E)/E0;
sc=0.1*smallstep/fE; % error function should decrease by about 10%
dx=sc*dx;
steps=sqrt(sum(dx.^2,2));
sc2=0.5/max(steps);
maxstep=max(steps);
if sc2<1,
    dx=sc2*dx;
    add_msg_board(sprintf('Reduced step size by factor %5.3f',sc2));
end;

% disp(sprintf('Original m.s.e. : %6.4f',E0));
network=network0+dx;
if ~isempty(DEER),
    DEER2=update_labels(DEER,network0,network);
end;
% compute actual error function
fom=0;
if mr>0,
    for k=1:mr, 
        vec=DEER2(k).xyz2-DEER2(k).xyz1;
        rp=norm(vec);
        fom=fom+weights(k)*(10*DEER2(k).r-rp)^2/(mr+md);    
    end;
end;
if md>0,
    for k=1:md, 
        vec=network(direct(k,2),:)-network(direct(k,1),:);
        rp=norm(vec);
        fom=fom+dir_weights(k)*(10*direct(k,4)-rp)^2/(mr+md);    
    end;
end;
if mi>0,
    for k=1:mi, 
        vec=network(internal(k,2),:)-network(internal(k,1),:);
        rp=norm(vec);
        fom=fom+internal(k,4)*(internal(k,3)-rp)^2/mi;    
    end;
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

function DEER=update_labels(DEER,network0,network)
% updates mean spin label cocordinates after a change of C_alpha
% coordinates of the network

global model

coarse=false;

[maxnum,m]=size(network);
scarce=0;
% update of label coordinates
for k=1:length(DEER),
    xyz=DEER(k).xyz1;
    l=DEER(k).res1;
    cindices=model.coarse(model.current_structure).indices; % actual indices of residues in network
    local_template=zeros(5,3);
    local_template_0=zeros(5,3);
    % make a local template to fit rotation and translation
    poi=0;
    if ~coarse,
        for kk=-2:2,
            if l+kk>0 && l+kk<=maxnum, % is addressed residue a network point?
                diff=cindices(l+kk,4)-cindices(l,4);
                if diff==kk, % is addressed residue part of a continuous segment?
                    poi=poi+1;
                    local_template(poi,:)=network(l+kk,:);
                    local_template_0(poi,:)=network0(l+kk,:);
                end;
            end;
        end;
    end;
    if poi>=3, % found sufficient number of points to determine local rotation and translation
        [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz=[xyz 1];
        xyz=transmat*xyz';
        xyz=xyz';
        DEER(k).xyz1=xyz(1:3);
    else
        DEER(k).xyz1=xyz+network(l,:)-network0(l,:);
        scarce=scarce+1;
    end;
    xyz=DEER(k).xyz2;
    l=DEER(k).res2;
    cindices=model.coarse(model.current_structure).indices; % actual indices of residues in network
    local_template=zeros(5,3);
    local_template_0=zeros(5,3);
    % make a local template to fit rotation and translation
    poi=0;
    if ~coarse,
        for kk=-2:2,
            if l+kk>0 && l+kk<=maxnum, % is addressed residue a network point?
                diff=cindices(l+kk,4)-cindices(l,4);
                if diff==kk, % is addressed residue part of a continuous segment?
                    poi=poi+1;
                    local_template(poi,:)=network(l+kk,:);
                    local_template_0(poi,:)=network0(l+kk,:);
                end;
            end;
        end;
    end;
    if poi>=3, % found sufficient number of points to determine local rotation and translation
        [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz=[xyz 1];
        xyz=transmat*xyz';
        xyz=xyz';
        DEER(k).xyz2=xyz(1:3);
    else
        DEER(k).xyz2=xyz+network(l,:)-network0(l,:);
        scarce=scarce+1;
    end;
end;
