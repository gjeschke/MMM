function [rmsd,network] = transition_by_ANM(modes,network0,target)
% function [rmsd,network] = transition_by_ANM(modes,network0,target)
%
% models a coordinate transition by propagation along a few slow modes of
% an anisotropic network model
%

% Input:
% modes     vector of mode numbers for modes used to propagate the network, the numbers
%           refer to the anisotropic network model of the current structure
%           in MMM: model.ANM(model.current_structure).u
% network0  matrix [n,3] of coordinates of network points (knots),
%           restraint variable DEER must refer labelled sites to network
%           knots (C_alpha coordinates of residues), these references are
%           given for the k-th restraint by DEER(k).res1 and DEER(k).res2
%           indexing of residues in network is supposed to conform to
%           model.coarse(model.current_structure).indices
% target    target network coordinates, must be same size as network0 and
%           must correspond
%
% Output:
% rmsd      remaining r.m.s.d. between fitted endpoint of the transition
%           and target structure
% network   the fitted network	for the transition endpoint
%           
% G. Jeschke, 2010

global model

rmsd_threshold=2e-5; % resolution in terms of C_alpha r.m.s.d. normalized to step size
step=3; % maximum step size 1 Å
maxtime=1800; % maximum time in seconds

[m,n]=size(network0); % m number of residues

[mt,nt]=size(target);

if mt~=m,
    add_msg_board('ERROR: Sizes of initial and target model do not agree.');
    rmsd=1e6;
    network=[];
    return
end;

sumcoeff=zeros(size(modes));
rmsd = rmsd_superimpose(target,network0);
rmsd0=1e6;

figure(1); clf;

x=0;
y=rmsd;

h=line(x,y,'Color','k');
axis([0,500,0,1.1*y]);

numit=0;
runtime=0;
tic;
while (rmsd0-rmsd)/step > rmsd_threshold && runtime<maxtime,
    rmsd0=rmsd;
    [grad,sgn,active]=get_gradients(rmsd0,modes,network0,target);
    if active==0, % network is at minimum w.r.t. all modes
        break;
    end;
    coeff=step*sgn.*grad/sum(grad);
    if numit==0,
        coeff0=coeff;
    end;
    numit=numit+1;
    [rmsd,network] = propagate_network(coeff,modes,network0,target);
    model.ANM(model.current_structure).u=reorient_ANM(model.ANM(model.current_structure).u,network0,network);
    sumcoeff=sumcoeff+coeff;
    x=[x numit];
    y=[y rmsd];
    set(h,'XData',x,'YData',y');
    drawnow;
    network0=network;
    runtime=toc;
end;

disp(sprintf('Fitting took %i iterations.',numit));
disp(sprintf('Last gradient sum is: %6.4f with %i active modes',sum(grad),active));
disp(sprintf('Last change in r.m.s.d. was %6.3f',rmsd-rmsd0));
disp(sprintf('Fitting required %5.1f s',runtime));
coeff0=max(sumcoeff)*coeff0/max(coeff0);
figure(2); clf;
plot(modes,sumcoeff,'k.');
hold on;
plot(modes,coeff0,'ro');


function [grad,sgn,active]=get_gradients(rmsd0,modes,network0,target)
% determines gradient vector for figure of merits for selected modes and
% a given network model and DEER restraints 
%

grad=zeros(1,length(modes));
sgn=zeros(1,length(modes));
for k=1:length(modes),
    coeff=1;
    rmsd1 = propagate_network(coeff,modes(k),network0,target);
    rmsd2 = propagate_network(-coeff,modes(k),network0,target);
    grad(k)=rmsd0-min([rmsd1,rmsd2]);
    if rmsd1<rmsd2,
        sgn(k)=coeff;
    else
        sgn(k)=-coeff;
    end;
end;
grad(grad<0)=0;
active=sum(grad>0);

function [rmsd,network]=propagate_network(coeff,modes,network0,target)

global model

[m,n]=size(network0);

% network propagation
network=network0; % initialize propagated network
for k=1:length(coeff),
    cmode=modes(k);
    evec=model.ANM(model.current_structure).u(:,cmode);
    mode=reshape(evec,3,m);
    network=network+coeff(k)*mode';
end;

rmsd = rmsd_superimpose(target,network);
