function ff =get_ff_higher_order(NOpos,time,depth,gridsize,silent)
% function ff =get_ff_higher_order(NOpos,t,depth)
% generates the DEER form factor for given NO midpoint coordinates and
% populations
% it is assumed that the conformation distributions of the labels are
% uncorrelated
% a maximum of five terms of the series expansion are computed
%
% NOpos     cell array of NO midpoint positions and populations, must
%           contain at least 3 spins
% time      time axis
% depth     total modulation depth
% gridsize  optional parameter for the orientation grid size, can be
%           'small','medium','large', or 'very_large', defaults to
%           'medium', also if an invalid grid size is requested
% silent    optional flag to prevent question about long computatioin time,
%           1 no question, 0 user is asked, defaults to 0
%
% based on G. Jeschke, M. Sajid, M. Schulte, A. Godt, Phys. Chem. Chem.
% Phys. 2009, 11, 6580-6591 and a correction of Eq. (5) by Tona von Hagens
%
% G. Jeschke, 2009

nm2MHz=52.04; % conversion factor between nm^(-3) and MHz

libration=0.1; % libration broadening in nm (std. dev.)

threshold_systematic=120; % fictive time in s that defines the limit for systematic computation
                         % above this limit, a Monte carlo computation is
                         % performed
                         
MC=100000; % number of Monte Carlo trials

ra=1;
re=20;
rax=ra:0.05:re;
distr=zeros(size(rax));

if nargin<4,
    gridsize='medium';
end;

if nargin<5,
    silent=0;
end;

[m,n]=size(time);
if m>n, time=time'; end;
ndat=length(time); % number of data points for the pair and three-spin contribution

ff=(1-depth)*ones(1,ndat); % initialize form factor (unmodulated part)

% load orientation grid
switch lower(gridsize),
    case 'small'
        load grid_small 
    case 'medium'
        load grid_medium 
    case 'large'
        load grid_large 
    case 'very_large'
        load grid_very_large 
    otherwise
        load grid_medium
end;
B0=grid(1:3,:);
weights=grid(4,:);
weights=weights'/sum(weights); % renormalize weights
nori=length(weights);

nspin=length(NOpos); % number of spins
lambda=1-(1-depth)^(1/(nspin-1)); % inversion efficiency by Eq. (6)

% Generate coefficients of series expansion, Eq. (5) with correction
coeff=zeros(1,nspin-1);
coeff0=coeff;
fact=1;
for k=1:nspin-1,
    fact=k*fact;
    coeff(k)=(1-lambda)^(nspin-1-k)*lambda^k/fact;
    coeff0(k)=(1-lambda)^(nspin-1-k)*lambda^k*nchoosek(nspin-1,k);
end;
% Determine how many terms need to be evaluated
terms=sum(cumsum(coeff0)<0.99*sum(coeff0))+1; % 99% of all modulation is accounted for
if terms>nspin-1, terms=nspin-1; end;
if terms>5,
    terms=5;
end;

add_msg_board(sprintf('Computing series expansion with %i of %i terms',terms,nspin-1));
add_msg_board(sprintf('Grid size %s with %i orientations.',gridsize,nori));
add_msg_board(sprintf('Total modulation depth of %5.3f assumed.',depth));

% initialize populations and dipolar evolution function arrays
maxrot=0;
maxpop=zeros(1,nspin);
for k=1:nspin,
    [m,n]=size(NOpos{k});
    if m>maxrot, maxrot=m; end;
    pop{k}=NOpos{k}(:,4);
    pop{k}=pop{k}/sum(pop{k});
    maxpop(k)=max(pop{k});
end;
maxpop_pair=0;
for k=1:nspin-1,
    for kk=k+1:nspin,
        if maxpop(k)*maxpop(kk)>maxpop_pair,
            maxpop_pair=maxpop(k)*maxpop(kk);
        end;
    end;
end;
pair_threshold=1e-3*maxpop_pair;
spin1=zeros(1,nspin*(nspin-1)*maxrot*(maxrot-1));
spin2=spin1;
rotamer1=spin1;
rotamer2=spin1;
pairpop=spin1;
poi=0;
total_population=0;
for k1=1:nspin-1, % first spin
    for r1=1:length(pop{k1}), % rotamers of first spin
        for k2=k1+1:nspin, % second spin
            for r2=1:length(pop{k2}), % rotamers of second spin
                cpairpop=pop{k1}(r1)*pop{k2}(r2);
                if cpairpop>pair_threshold,
                    poi=poi+1;
                    spin1(poi)=k1;
                    rotamer1(poi)=r1;
                    spin2(poi)=k2;
                    rotamer2(poi)=r2;
                    pairpop(poi)=cpairpop;
                    total_population=total_population+pairpop(poi);
                end;
            end;
        end;
    end;
end;
pairs=poi;
pairpop=pairpop(1:poi);
pairpop=pairpop/sum(pairpop);
% pairpop=pairpop/(nspin-2); % renormalization
spin1=spin1(1:poi);
spin2=spin2(1:poi);
rotamer1=rotamer1(1:poi);
rotamer2=rotamer2(1:poi);



ftrace=zeros(nspin-1,ndat);
estim=0;
for k=1:terms,
    estim=estim+pairs^k;
end;
update=1+round(5000/estim);
est_time=3*estim*ndat/(350000*nspin);
% add_msg_board(sprintf('Estimated computation time: %6.1fs',est_time));

% if ~silent && est_time>60,
%     answer=questdlg('Do you want really want to compute higher-order correlations?',sprintf('Computation will take very long about %is',round(est_time)),'Compute','Cancel','Cancel');
%     if ~strcmp(answer,'Compute'),
%         ff=[];
%         add_msg_board('Computation of higher-order correlations cancelled by user.');
%         return
%     end;
% end;

allr=zeros(1,pairs);
allrvec=zeros(pairs,3);
for k=1:pairs, % initialize dipolar evolution function arrays
    NO1=NOpos{spin1(k)}(rotamer1(k),1:3)/10; % convert to nm
    NO2=NOpos{spin2(k)}(rotamer2(k),1:3)/10;
    rvec=NO1-NO2;
    r=norm(rvec);
    allr(k)=r;
    allrvec(k,:)=rvec/r;
end;

comp_status=status_figure('Close to stop.');
simcoeff=zeros(1,nspin-1);

tic,
if est_time>threshold_systematic, % Monte Carlo computation
    add_msg_board('Monte Carlo simulation');
    [sortpop,poppoi]=sort(pairpop,'descend');
    popcumsum=cumsum(sortpop);
    dipevo=zeros(nspin-1,ndat); % initialize dipolar evolution traces
    for trial=1:MC,
        if mod(trial,100)==0, 
            comp_status=status_figure((trial-1)/MC,comp_status);
            drawnow;
            if ~comp_status, 
                ff=[];
                add_msg_board('Computation of higher-order correlations cancelled by user.');
                return 
            end;
        end;
        ct=rand;
        st=sqrt(1-ct^2);
        fi=2*pi*rand;
        B0vec=[cos(fi)*st,sin(fi)*st,ct];
        % make a vector of random unique pair selections
        p=(1+sum(popcumsum<rand))*ones(1,nspin-1);
        for k1=2:nspin-1,
            test=1;
            while test,
                test=0;
                p(k1)=(1+sum(popcumsum<rand));
                for k2=1:k1-1,
                    if p(k2)==p(k1); test=1; end; % no pair can occur twice
                end;
                if spin1(k1)~=spin1(1); test=1; end; % pairs must agree in observer spin
            end;
        end;
        for k1=1:nspin-1,
            p(k1)=poppoi(p(k1));
        end;
        for k=1:nspin-1, % initialize dipolar evolution function arrays
            r=allr(p(k))+libration*randn;
            if r<2,
                w1=1-((2-r)/0.5)^3; % dampening of short distances, currently unused
            else
                w1=1;
            end;
            if r<1.5,
                w1=0;
            end;
            revec=allrvec(p(k),:);
            ct=sum(B0vec.*revec); % cosine of angle between B0 and spin-spin vector
            wdd=(3*ct^2-1)*2*pi*nm2MHz/r^3; % dipolar frequency, 1 nm corresponds to 52.04 MHz
            dipevo(k,:)=cos(wdd*time);
            rpoi=1+round((length(rax)-1)*(r-ra)/(re-ra));
            if rpoi>0 && rpoi<=length(rax),
                distr(rpoi)=distr(rpoi)+1;
            end;
        end;
        for t=1:terms, % loop over terms
            ltrace=zeros(1,ndat); % initialize trace
            switch t
                case 1
                    for p1=1:nspin-1,
                        ltrace=ltrace+dipevo(p1,:);
                        simcoeff(1)=simcoeff(1)+1/MC;
                    end;
                case 2
                    for p1=1:nspin-1,
                        for p2=1:nspin-1,
                            if p1~=p2,
                                ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:);
                                simcoeff(2)=simcoeff(2)+1/MC;
                            end;
                        end;
                    end;
                case 3
                    for p1=1:nspin-1,
                        for p2=1:nspin-1,
                            if p1~=p2,
                                for p3=1:nspin-1,
                                    if p3~=p1 && p3~=p2,
                                        ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:).*dipevo(p3,:);
                                        simcoeff(3)=simcoeff(3)+1/MC;
                                    end;
                                end;
                            end;
                        end;
                    end;
                case 4
                    for p1=1:nspin-1,
                        for p2=1:nspin-1,
                            if p1~=p2,
                                for p3=1:nspin-1,
                                    if p3~=p1 && p3~=p2,
                                        for p4=1:nspin-1,
                                            if p4~=p3 && p4~=p2 && p4~=p1,
                                                ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:).*dipevo(p3,:).*dipevo(p4,:);
                                                simcoeff(4)=simcoeff(4)+1/MC;
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
                case 5
                    for p1=1:nspin-1,
                        for p2=1:nspin-1,
                            if p1~=p2,
                                for p3=1:nspin-1,
                                    if p3~=p1 && p3~=p2,
                                        for p4=1:nspin-1,
                                            if p4~=p3 && p4~=p2 && p4~=p1,
                                                for p5=1:nspin-1,
                                                    if p5~=p4 && p5~=p3 && p5~=p2 && p5~=p1,
                                                        ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:).*dipevo(p3,:).*dipevo(p4,:).*dipevo(p5,:);
                                                        simcoeff(5)=simcoeff(5)+1/MC;
                                                    end;
                                                end;
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
            end;
            ftrace(t,:)=ftrace(t,:)+ltrace*coeff0(t)/MC;
        end;
    end;
else % systematic computation
    add_msg_board('Systematic sampling of conformation space');
    dipevo=zeros(pairs,ndat); % initialize dipolar evolution traces
    for ori=1:nori,
        if mod(ori,update)==0, 
            comp_status=status_figure((ori-1)/nori,comp_status);
            drawnow;
            if ~comp_status, 
                ff=[];
                add_msg_board('Computation of higher-order correlations cancelled by user.');
                return 
            end;
        end;
        B0vec=B0(:,ori)';
        w=weights(ori);
        poi=0;
        for k=1:pairs, % initialize dipolar evolution function arrays
            r=allr(k)+libration*randn;
            if r<2,
                w1=1-((2-r)/0.5)^3; % dampening of short distances, currently unused
            else
                w1=1;
            end;
            if r<1.5,
                w1=0;
            end;
            revec=allrvec(k,:);
            ct=sum(B0vec.*revec); % cosine of angle between B0 and spin-spin vector
            wdd=(3*ct^2-1)*2*pi*nm2MHz/r^3; % dipolar frequency, 1 nm corresponds to 52.04 MHz
            dipevo(k,:)=pairpop(k)*cos(wdd*time);
            rpoi=1+round((length(rax)-1)*(r-ra)/(re-ra));
            if rpoi>0 && rpoi<=length(rax),
                distr(rpoi)=distr(rpoi)+pairpop(k);
            end;
        end;
        for t=1:terms, % loop over terms
            ltrace=zeros(1,ndat); % initialize trace
            switch t
                case 1
                    for p1=1:pairs,
                        ltrace=ltrace+dipevo(p1,:);
                        simcoeff(1)=simcoeff(1)+w*pairpop(p1);
                    end;
                case 2
                    for p1=1:pairs,
                        for p2=1:pairs,
                            if p1~=p2 && spin1(p1)==spin1(p2),
                                ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:);
                                simcoeff(2)=simcoeff(2)+w*pairpop(p1)*pairpop(p2);
                            end;
                        end;
                    end;
                case 3
                    for p1=1:pairs,
                        for p2=1:pairs,
                            if p1~=p2 && spin1(p1)==spin1(p2),
                                for p3=1:pairs,
                                    if p3~=p1 && p3~=p2 && spin1(p1)==spin1(p3),
                                        ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:).*dipevo(p3,:);
                                        simcoeff(3)=simcoeff(3)+w*pairpop(p1)*pairpop(p2)*pairpop(p3);
                                    end;
                                end;
                            end;
                        end;
                    end;
                case 4
                    for p1=1:pairs,
                        for p2=1:pairs,
                            if p1~=p2 && spin1(p1)==spin1(p2),
                                for p3=1:pairs,
                                    if p3~=p1 && p3~=p2 && spin1(p1)==spin1(p3),
                                        for p4=1:pairs,
                                            if p4~=p3 && p4~=p2 && p4~=p1 && spin1(p1)==spin1(p4),
                                                ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:).*dipevo(p3,:).*dipevo(p4,:);
                                                simcoeff(4)=simcoeff(4)+w*pairpop(p1)*pairpop(p2)*pairpop(p3)*pairpop(p4);
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
                case 5
                    for p1=1:pairs,
                        for p2=1:pairs,
                            if p1~=p2 && spin1(p1)==spin1(p2),
                                for p3=1:pairs,
                                    if p3~=p1 && p3~=p2 && spin1(p1)==spin1(p3),
                                        for p4=1:pairs,
                                            if p4~=p3 && p4~=p2 && p4~=p1 && spin1(p1)==spin1(p4),
                                                for p5=1:pairs,
                                                    if p5~=p4 && p5~=p3 && p5~=p2 && p5~=p1 && spin1(p1)==spin1(p5),
                                                        ltrace=ltrace+dipevo(p1,:).*dipevo(p2,:).*dipevo(p3,:).*dipevo(p4,:).*dipevo(p5,:);
                                                        simcoeff(5)=simcoeff(5)+w*pairpop(p1)*pairpop(p2)*pairpop(p3)*pairpop(p4)*pairpop(p5);
                                                    end;
                                                end;
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
            end;
            ftrace(t,:)=ftrace(t,:)+w*ltrace*coeff0(t);
        end;
    end;
end;
ttrace=zeros(1,ndat);
% disp(simcoeff);
% figgy=gcf;
% figure(15); clf;
% plot(ftrace(1,:)/simcoeff(1),'k');
% figure(16); clf;
% plot(ftrace(2,:)/simcoeff(2),'r');
% figure(17); clf;
% plot(rax,distr,'k');
% figure(figgy);
for k=1:terms,
    ttrace=ttrace+ftrace(k,:)/simcoeff(k);
end;
if comp_status, status_figure(1); end;
% t0=toc;
ff=ff+ttrace;
ff=ff/max(ff);

% add_msg_board(sprintf('Actual computation time: %6.1fs',t0));
