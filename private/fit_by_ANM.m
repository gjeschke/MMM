function [fom,rmsd,network,dxmat,DEER,best,diagnostics] = fit_by_ANM(network0,DEER,direct,test_mode,correspondence,target,model,ENM_param)
% function [fom,rmsd,network,dxmat,DEER,best] = fit_by_ANM(network0,DEER,direct,test_mode,correspondence,target)
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
% test_mode optional flag, if true, a target Calpha coordinate array must 
%           be provided as argument target and a correspondence table to 
%           network0, information is then displayed on fit
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
% best      information on best approach to target structure (in test_mode
%           only)
%           best.it     iteration cycle
%           best.rmsd   structure rmsd
%           best.drmsd  distance rmsd
%
% diagnostics   diagnostics.x   vector with iteration numbers
%               diagnostics.y   vector with constraint rmsd values for each
%                               iteration
%               diagnostics.d2  2nd derivative of constraint r.m.s.d.
%               diagnostics.d2s smoothed 2nd derivative of constraint
%                               r.m.s.d.
%               diagnostics.ys  r.m.s.d. to target structure for each
%                               iteration (only in test mode, otherwise empty)
%               diagnostics.converged   iteration number where fit is
%                                       assumed to have converged
%               .rmsd0          initial target structure r.m.s.d.
%               .rmsd           final target structure r.m.s.d.
%               .drmsd0         initial constraint r.m.s.d.
%               .drmsd          final constraint r.m.s.d.
%
%
% NOTE that further parameters are defined in global variable ENM_param,
% which is initialized in initialize_MMM.m
%           
% G. Jeschke, 2010-2011

diagnostics.x=[];
diagnostics.y=[];
diagnostics.d2=[];
diagnostics.ds=[];
diagnostics.ys=[];
diagnostics.converged=0;
diagnostics.rmsd0=[];
diagnostics.drmsd0=[];
diagnostics.rmsd=[];
diagnostics.drmsd0=[];

old_basis=[];

modes=7:6+ENM_param.fit_basis;

if isfield(model.ANM(model.current_structure),'u') && ~isempty(model.ANM(model.current_structure).u),
    basis=model.ANM(model.current_structure).u(:,1:max(modes));
    lbasis=model.ANM(model.current_structure).lambda(1:max(modes));
    u=model.ANM(model.current_structure).u;
else
    Hessian=setup_ANM_bonded_parallel(network0,model,ENM_param);
    clear u
    [u,d]=eig(Hessian);
    clear Hessian
    basis=u(:,1:max(modes));
    lbasis=zeros(1,max(modes));
    for k=1:max(modes),
        lbasis(k)=d(k,k);
    end;
    clear d
end;

if nargin<4,
    test_mode=false;
end;
maximum_likelihood=true;

f=ENM_param.fix_force; % (relative) force constant for restraining local structure
f2=ENM_param.sec_force; % (relative) force constant for restraining local structure
                        % for second neighbor

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
weights=ones(size(weights));
dir_weights=ones(size(dir_weights));
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
x=0;
if ~isempty(DEER),
    for k=1:mr, % set up unit vector between restraint residue pairs
        vec=DEER(k).xyz2-DEER(k).xyz1;
        rp=norm(vec);
        rmsd=rmsd+(10*DEER(k).r-rp)^2;    
    end;
    rmsd=sqrt(rmsd/mr);
    rmsd0=rmsd;
    add_msg_board(sprintf('Initial DEER distance constraint r.m.s.d. %6.2f',rmsd));
    best.initial_DEER=rmsd;
    diagnostics.drmsd0=rmsd;
else
    rmsd=0.001;
end;
fit_limit=rmsd/4+3*ENM_param.tol_dist/4;
if fit_limit<ENM_param.tol_dist + rmsd/8 && ~isempty(DEER),
    fit_limit=ENM_param.tol_dist + rmsd/8;
end;
fit_limit=0.2;
best.fit_limit=fit_limit;


add_msg_board(sprintf('Fit limit is %5.2f Å',fit_limit));
y=rmsd;

figure(1); clf;
title('Constraint r.m.s.d./ Figure of merit');

if md>0,
    rmsd_dir=0;
    for k=1:md, % set up unit vector between restraint residue pairs
        vec=network0(direct(k,2),:)-network0(direct(k,1),:);
        rp=norm(vec);
        rmsd_dir=rmsd_dir+(10*direct(k,4)-rp)^2;    
    end;
    rmsd_dir=sqrt(rmsd_dir/md);
    add_msg_board(sprintf('Initial direct distance constraint r.m.s.d. %6.2f',rmsd_dir));
    best.initial_direct=rmsd_dir;
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
axis([0,ENM_param.cycles,0,1.1*max([y,yf,ydir])]);


if test_mode,
    rmsd_struct=rmsd_superimpose(target(correspondence(2,:),:),network0(correspondence(1,:),:));

    ys=rmsd_struct;
    diagnostics.rmsd0=rmsd_struct;
    
%     fprintf(1,'Initial structure r.m.s.d. is %5.2f Å.\n',ys);
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
struct_step=zeros(1,ENM_param.cycles);
change=zeros(1,ENM_param.cycles);
energy_change=zeros(1,ENM_param.cycles);
scaling=zeros(1,ENM_param.cycles);
dx_step=zeros(1,ENM_param.cycles);
fom_vec=zeros(1,ENM_param.cycles);
% mode_spectrum=zeros(100,length(modes));
% mode_similarity=zeros(100,length(modes));
best_network=network0;
maxsc=0;

fom_init=fom0;
network00=network0;
no_stop=10;
ldfom_vec=zeros(1,ENM_param.cycles);
no_sat=ENM_param.cycles;
if no_sat>100,
    no_sat=100;
end;

reorientate=ENM_param.reorientate; % save initial setting
if ~reorientate 
    if ~ENM_param.diagonalize,
        add_msg_board(sprintf('Threshold for switching to normal mode reorientation: %4.1f Å',threshold));
    else
        add_msg_board('Normal modes are updated by recomputation and diagonalization of Hessian.');
    end;
else
    add_msg_board('Normal modes are updated by reorientation.');
end;

limit_found=false;
half_found=false;
tol_found=false;
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
    [fom,dx,l2d,network,DEER2,maxstep,mode_weights,sc,basis,lbasis,fE,lambda,u]=linear_regression(modes,network0,DEER,internal,weights,direct,dir_weights,basis,lbasis,old_basis,model,ENM_param);
    if ~ENM_param.diagonalize && fom<threshold,
        ENM_param.reorientate=true;
    end;
    if fE<0,
        network=network0;
        break;
    end;
        
    bas=3*numit;
    dxmat(:,bas+1:bas+3)=dx;
    yf=[yf fom];
    dfom=fom0-fom;
    numit=numit+1;
    energy_change(numit)=fE;
    ldfom_vec(numit)=log(fom0)-log(fom);
%     fprintf(1,'Saturation parameter: %6.4f\n',ldfom_vec(numit)/max(ldfom_vec(1:numit)));
%     if ~sat_log && ldfom_vec(numit)<0.01*max(ldfom_vec(1:numit)), % see Zheng/Brooks, p. 4329, Saturation detection
%         if no_sat>10,
%             no_sat=10;
%         end;
%         sat_log=true;
%         add_msg_board(sprintf('Total energy function saturation detected at step %i',numit));
%     end;
    fom_vec(numit)=fom;
    scaling(numit)=sc;
    dx_step(numit)=maxstep;
    if ENM_param.relax,
        [rmsd,maxstep,network]=relax_structure(network,internal,numit,u,lambda,ENM_param);
    end;
    networks{numit}=network;
    struct_step(numit)=rmsd_superimpose(network,network0);
    change(numit)=rmsd_superimpose(network,network00);
    if sc>maxsc,
        maxsc=sc;
        best_network=network;
        if test_mode,
            stop_rmsd=rmsd_superimpose(target(correspondence(2,:),:),network(correspondence(1,:),:));
            stop_it=numit;
        end;
    end;
    DEER=DEER2;
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
        if rmsd<fit_limit,
            no_stop=0;
            no_sat=0;
        end;
        ddis=rmsd0-rmsd;
        if ddis>ddmax, ddmax=ddis; end;
        rmsd0=rmsd;
%         if ddis<0.05*ddmax, 
%             satflag=true;
%             add_msg_board(sprintf('DEER constraint saturation detected at step %i',numit));
%         end;
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
        rmsd_struct=rmsd_superimpose(target(correspondence(2,:),:),network0(correspondence(1,:),:));
        if rmsd_struct<best_fit,
            best_fit=rmsd_struct;
            best_it=numit;
            best_dist=rmsd;
            best.rmsd=rmsd_struct;
            best.it=numit;
            best.drmsd=rmsd;
            best_net=network0;
            best_l2d=l2d;
        end;
        if rmsd<fit_limit && ~limit_found,
            limit_found=true;
            best.limit_rmsd=rmsd_struct;
            add_msg_board(sprintf('At fit limit of %5.2f Å the structure r.m.s.d. is %5.2f Å',fit_limit,rmsd_struct));
        end;
        ys=[ys rmsd_struct];
        set(hs,'XData',x,'YData',ys);
        axis([0,nnx,0,1.1*max(ys)]);
    end;
    if rmsd<fit_limit,
        limit_found=true;
    end;
    drawnow;
    runtime=toc;
end;

if ~limit_found,
    add_msg_board('Constraint limit not reached.');
    if test_mode,
        best.limit_rmsd=rmsd_struct;
    end;
end;

ENM_param.reorientate=reorientate; % restore initial setting

best_guess=numit;

if ~isempty(DEER),
    rmsd_DEER=y(1:numit);
end;
xx=1:numit;
change=change(1:numit);

if runtime>=maxtime,
    add_msg_board('Fit stopped as maximum time was exceeded.');
end;
if rmsd<=fit_limit,
    add_msg_board('Fit r.m.s.d. fell below limit determined by specified distance ranges');
end;
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
    if test_mode,
        [minys,poi]=min(ys);
        if poi<numit,
%             fprintf(1,'Overfitting started at iteration %i and r.m.s.d. %5.2f Å.\n',poi,minys);
%             fprintf(1,'Constraint r.m.s.d. was %5.2f Å\n',ydir(poi));
%             fprintf(1,'Energy change was %5.1f%%\n',100*energy_change(poi));
        end;
    end;
end;

if test_mode && ~isempty(DEER),
    ys=ys(1:numit);
    y=y(1:numit);
%     figure(4); clf;
%     plot(y,ys,'k');
%     title('Structure r.m.s.d. vs. distance r.m.s.d.');
%     set(gca,'FontSize',14);
%     hold on
%     plot([best_guess best_guess],[0,max(ys)],'r:');
    axis([0,1.05*max(y),0,1.05*max(ys)]);
    dx_step=dx_step(1:numit);
    ys=ys(1:numit);
    xx=1:numit;
    
    
%     figure(6); clf;
%     plot(xx,dx_step,'k');
%     hold on;
%     plot(xx,ys,'r');
%     plot([best_guess best_guess],[0,max(ys)],'r:');
    add_msg_board(sprintf('Final structure r.m.s.d. %5.2f\n',rmsd_struct));
end;

diagnostics.x=xx;
diagnostics.y=yf;

struct_step=struct_step(1:numit);
% figure(5);
% clf;
% plot(struct_step,'k');
% title('R.m.s.d. between consecutive structures');
% hold on;
% plot([best_guess best_guess],[0,max(struct_step)],'r:');

% figure(6);
% clf;
% plot(scaling,'k');
% title('Step scaling factor');
% hold on
% [masc,best]=max(scaling);
% plot([best_guess,best_guess],[0,masc],'r:');
% fprintf(1,'Extracting structure at iteration %i.\n',best_guess);
%rmsd=stop_rmsd;
network=best_network;

figure(1);
axis([0,numit,0,1.05*max(yf)]);
set(gca,'FontSize',14);
hold on;
plot([best_guess best_guess],[0,max(yf)],'r:');

diagnostics.converged=numit;

diagnostics.drmsd=yf(numit);
if test_mode,
    figure(2);
    axis([0,numit,0,1.05*max(ys)]);
    set(gca,'FontSize',14);
    hold on;
    plot([best_guess best_guess],[0,max(ys)],'r:');
    final_srmsd=ys(numit-1);
    diagnostics.rmsd=final_srmsd;
    diagnostics.ys=ys;
end;

if test_mode,
    network=networks{best_guess};
    rmsd_struct=rmsd_superimpose(target(correspondence(2,:),:),network(correspondence(1,:),:));
%     fprintf(1,'Best r.m.s.d. during whole fit was %5.2f at iteration %i.\n',best_fit,best_it);
    if ~isempty(DEER),
%         fprintf(1,'At this point the DEER constraint r.m.s.d. was %5.2f Å.\n',y(best_it));
%         fprintf(1,'Last r.m.s.d. was %5.2f at iteration %i.\n',ys(numit),numit);
%         fprintf(1,'Selected structure has r.m.s.d. of %5.2f at iteration %i.\n',rmsd_struct,best_guess);
    end;
    add_msg_board(sprintf('Final fit %4.2f Å with d.r.m.s.d. %4.2f Å.',rmsd_struct,rmsd));
else
    add_msg_board(sprintf('Fitting took %i iterations.',numit));
end;
add_msg_board(sprintf('Fitting required %5.1f s',runtime));
% network=best_net;

if md>0,
    rmsd_dir=0;
    for k=1:md, % set up unit vector between restraint residue pairs
        vec=network(direct(k,2),:)-network(direct(k,1),:);
        rp=norm(vec);
        rmsd_dir=rmsd_dir+(10*direct(k,4)-rp)^2;    
    end;
    rmsd_dir=sqrt(rmsd_dir/md);
    add_msg_board(sprintf('Final direct distance constraint r.m.s.d. %6.2f',rmsd_dir));
end;

if ~test_mode,
    rmsd_change=max(rmsd_DEER)-min(rmsd_DEER);
    [mi,best_it]=min(abs(change-rmsd_change));
end;

network=networks{numit};

if ~isempty(DEER),
    DEER=update_labels(DEER,network0,network,model);
    rmsd=rmsd_DEER(numit);
else
    rmsd=rmsd_dir;
end;
% compute actual error function
fom=0;
weights=ones(size(weights));
dir_weights=ones(size(dir_weights));
if mr>0,
    for k=1:mr, 
        vec=DEER(k).xyz2-DEER(k).xyz1;
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
        % fom=fom+internal(k,4)*(internal(k,3)-rp)^2/mi;    
    end;
end;
% disp(sprintf('New      m.s.e. : %6.4f',fom));
fom=sqrt(fom);



function [fom,dx,l2d,network,DEER2,maxstep,mode_weights,sc,basis,lbasis,fE,lambda,u]=linear_regression(modes,network,DEER,internal,weights,direct,dir_weights,basis,lbasis,old_basis,model,ENM_param)
% sets up and solves linear regression problem for displacement vector
% see Eq. (5), (7), (9), (10) of Zheng/Brooks

reorientate=ENM_param.reorientate;

transrot1=false;

if ENM_param.track,
    if ~isempty(old_basis),
        track=true;
    else
%         fprintf(1,'Tracking initialized.\n');
        track=false;
    end;
else
    track=false;
end;

md=0;
if nargin>5,
    if ~isempty(direct),
        [md,nd]=size(direct);
    end;
end;

if ~reorientate,
    if isfield(model.ANM(model.current_structure),'u'),
        model.ANM(model.current_structure).u=[];
    end;
    Hessian=setup_ANM_bonded_parallel(network,model,ENM_param);
    clear u
    [u,d]=eig(Hessian);
    clear Hessian
    basis=u(:,1:max(modes));
    lbasis=zeros(1,max(modes));
    for k=1:max(modes),
        lbasis(k)=d(k,k);
    end;
    [mm,nn]=size(u);
    lambda=zeros(1,mm);
    for k=1:mm,
        lambda(k)=d(k,k);
    end;
    clear d
end;

% new_modes=modes;
% rmsd_vec=modes;
% if track, % mode tracking
%     [m,n]=size(u);
%     assigned=zeros(1,m);
%     for k=1:length(modes),
%         minrmsd=1e6;
%         best=0;
%         old=reshape(old_basis(:,k),3,m/3);
%         for kk=7:6+2*length(modes),
%             if ~assigned(kk),
%                 rmsd=rmsd_superimpose(old,reshape(u(:,kk),3,m/3));
%                 if rmsd<minrmsd,
%                     best=kk;
%                     minrmsd=rmsd;
%                 end;
%             end;
%         end;
%         new_modes(k)=best;
%         rmsd_vec(k)=minrmsd;
%         assigned(best)=1;
%     end;
% end;
% if ENM_param.track,
%     modes=new_modes;
%     fprintf(1,'New mode basis: ');
%     for k=1:length(new_modes),
%         fprintf(1,'%i (%8.6f)',new_modes(k),rmsd_vec(k));
%     end;
%     fprintf(1,'\n');
%     old_basis=u(:,modes);
% end;


% [m,n]=size(d);
% lambda=zeros(1,m);
% for k=1:m,
%     lambda(k)=d(k,k);
% end;
% clear d

smallstep=0.05; % defines an incremental step within the linear regime

network0=network;
[m,n]=size(network);
% [m1,n1]=size(smallnet);

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
    evec=basis(:,modes(j));
    mode=reshape(evec,3,m);
    wm=lbasis(modes(j));
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
    evec=basis(:,modes(j));
    mode=reshape(evec,3,m);
    mode=mode';
    wm=lbasis(modes(j));
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
    DEER2=update_labels(DEER,network0,network,model);
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

% fprintf(1,'Old energy function %5.2f\n',E0)
% fprintf(1,'Predicted new energy function %5.2f\n',E)
l2d=lE/(lE+dE);
fE=(E0-E)/E0;
sc=0.1/(smallstep*fE); % error function should decrease by about 10%
if sc>1, sc=1; end;
if sc<0, sc=0.1; end;
dx=sc*dx;
steps=sqrt(sum(dx.^2,2));
% sc=1/max(steps);
% if reorientate,
%     sc2=0.2/max(steps);
% else
%     sc2=0.5/max(steps);
% end;

maxstep=max(steps);
if sc<1,
    dx=sc*dx;
    add_msg_board(sprintf('Reduced step size by factor %5.3f',sc));
end;

% disp(sprintf('Original m.s.e. : %6.4f',E0));
network=network0+dx;
if ~isempty(DEER),
    DEER2=update_labels(DEER,network0,network,model);
end;


basis=reorientate_eigenvectors(network,network0,basis);

% compute actual energy function change
E=0;
E0=0;
lE=0;
dE=0;
network=network0+smallstep*dx;
if isempty(DEER),
    DEER2=DEER;
else
    DEER2=update_labels(DEER,network0,network,model);
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
fE=(E0-E)/E0;
% fprintf(1,'Old energy function %5.2f\n',E0)
% fprintf(1,'New energy function %5.2f\n',E)
% fprintf(1,'Energy function decreased by %5.1f%%\n',100*fE)

% compute actual error function
fom=0;
weights=ones(size(weights));
dir_weights=ones(size(dir_weights));
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
        % fom=fom+internal(k,4)*(internal(k,3)-rp)^2/mi;    
    end;
end;
% disp(sprintf('New      m.s.e. : %6.4f',fom));
fom=sqrt(fom);

function [rmsd,maxstep,network]=relax_structure(network,internal,iteration,u,lambda,ENM_param)
% sets up and solves linear regression problem for local structure
% relaxation
% see Eq. (5), (7), (9), (10) of Zheng/Brooks and remarks at the end of p.
% 4329


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
% fprintf(1,'Original rmsd : %5.3f at iteration %i',rmsd0,iteration);
% fprintf(1,'Targeted rmsd : %5.3f and maximum step: %4.2f',tol_local,maxstep);
% fprintf(1,'Achieved rmsd : %5.3f with %i iterations',rmsd,it);

function DEER=update_labels(DEER,network0,network,model)
% updates mean spin label cocordinates after a change of C_alpha
% coordinates of the network

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

% function Hessian=setup_ANM_bonded(Ca_coor,model,ENM_param)
% % function Hessian=setup_ANM_bonded(Ca_coor)
% %
% % returns the Hessian for a fully connected anisotropic network model (ANM)
% % when given the Calpha coordinates
% % the force constants have inverse 6th power dependence on interresidue
% % distance, except for next neighbors (i,i+1) and second neighbors (i,i+2)
% % along the polypetide chain
% %
% % see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% % Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% % Membrane Proteins. 
% % Chem. Rev. 110 (2010) 1463-1497.
% % section 2.3.1
% % and:
% % L. Yang, G. Song, R L. Jernigan, Proc. Natl. Acad. Sci. 106 (2009)
% % 12347–12352.
% %
% % the force constant at a distance of 1 Å p_ANM_gamma is defined in
% % initialize_MMM.m  and stored in global variable ENM_param
% %
% % an empty matrix is returned if the routine fails (no connection between
% % residues)
% %
% % G. Jeschke, 2010
% 
% % Hessian=setup_ANM_Hinsen(Ca_coor);
% % return
% 
% test_mode=false;
% 
% bond_force_1=10000*3.8^(-6); % 1200
% bond_force_2=10000*6^(-6); % 80
% exponent=6;
% 
% [m,n]=size(Ca_coor); % m is number of Calpha atoms
% 
% Hessian=zeros(3*m);
% r1=zeros(1,m);
% p1=0;
% r2=zeros(1,m);
% p2=0;
% for i=1:m-1,
%     basi=3*(i-1);
%     for j=i+1:m,
%         basj=3*(j-1);
%         rvec=Ca_coor(i,:)-Ca_coor(j,:); % distance vector of Calpha atoms
%         r=norm(rvec);
%         cind0=model.coarse(model.current_structure).indices(i,:);
%         cind1=model.coarse(model.current_structure).indices(j,:);
%         if sum(abs(cind1(1:3)-cind0(1:3)))==0, % both C_alpha in the same chain
%             diff=cind1(4)-cind0(4);
%             switch diff
%                 case 1
%                     cgamma=-bond_force_1*ENM_param.p_ANM_gamma;
%                     p1=p1+1;
%                     r1(p1)=r;
%                 case 2
%                     cgamma=-bond_force_2*ENM_param.p_ANM_gamma;
%                     p2=p2+1;
%                     r2(p2)=r;
%                 otherwise
%                     cgamma=-ENM_param.p_ANM_gamma*r^(-exponent);
%             end;
%         else
%             cgamma=-ENM_param.p_ANM_gamma*r^(-exponent);
%         end;
%         if cgamma~=0,
%             submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
%             Hessian(basi+1:basi+3,basj+1:basj+3)=cgamma*submat;
%             Hessian(basj+1:basj+3,basi+1:basi+3)=cgamma*submat;
%         end;
%     end;
% end;
% 
% if test_mode,
%     r1=r1(1:p1);
%     r2=r2(1:p2);
%     x=2:0.1:8;
% %     disp(sprintf('First  neighbor mean distance is %5.3f with std. dev. %5.3f',mean(r1),std(r1))); 
% %     disp(sprintf('Second neighbor mean distance is %5.3f with std. dev. %5.3f',mean(r2),std(r2))); 
%     figure(7); clf;
%     hold on;
%     hist(r1,x);
%     hist(r2,x);
% end;
% 
% conn=0;
% for i=1:m,
%     basi=3*(i-1);
%     submat=zeros(3);
%     for j=1:m,
%         basj=3*(j-1);
%         submat=submat+Hessian(basi+1:basi+3,basj+1:basj+3);
%     end;
%     conn=conn+trace(submat);
%     Hessian(basi+1:basi+3,basi+1:basi+3)=-submat;
% end;
% 
% if conn==0, % ANM is unconnected
%     Hessian=[];
% end;
