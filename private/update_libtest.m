function testset = update_libtest(testset,fname,handles,library,no_context)
% Updates a rotamer library test data set by loading experimental
% structures, and computing rotamers and distance distributions
% only <rsim> and <sigr_sim> are updated
% all required PDB structure files must exist on the Matlab path
%
% testset   test data set as obtained by rd_rotlib_test_data
% fname     output file name, log file name is derived from that
% library   (optional) name of rotamer library file, otherwise taken from
%           test data set
% no_context    optional flag that requests context-free computation,
%               defaults to false
%
% G. Jeschke, 2012/2015

global model
global hMain

set(hMain.figure,'Pointer','watch');

[path,name,ext] = fileparts(fname);
name=strcat(name,'.log');
fname2=fullfile(path,name);

if nargin>3 && ~isempty(library),
    testset.library=library;
end;

if ~exist('no_context','var'),
    no_context = false;
end;

fid=fopen(fname2,'wt');
fprintf(fid,'Test logfile for library %s\n\n',testset.library);

cpdb='****';
processed=false;
msg.error=0;
msg.text='No error.';

experiment=zeros(1,length(testset.data));
simulation=zeros(1,length(testset.data));

he_exp=experiment;
he_sim=simulation;
hep=0;

m_exp=experiment;
m_sim=simulation;
mp=0;

cps=zeros(1,100);
exp_ps=zeros(100,length(testset.data));
sim_ps=zeros(100,length(testset.data));
spoi=0;

llist={};
tstart=tic;
for k=1:length(testset.data),
    clear llist
    md=testset.data(k);
    if ~strcmpi(md.pdb,cpdb),
        spoi=spoi+1;
        stags{spoi}=md.pdb;
        cpoi=0;
        processed=false;
        cpdb=md.pdb;
        if md.cmode==-1,
            fname=strcat(cpdb,'.pdb1');
        else
            fname=strcat(cpdb,'.pdb');
        end;
        msg=new_pdb(fname);
        if msg.error,
            add_msg_board(sprintf('Problem reading file %s',fname));
            add_msg_board(msg.text);
        end
    end;
    switch md.cmode
        case {0,-1}
            llist{1}=md.res1;
            llist{2}=md.res2;
        case 1
            poi1=strfind(md.res1,'}');
            if ~isempty(poi1),
                res1=md.res1(poi1+1:end);
            else
                res1=md.res1;
            end;
            llist{1}=strcat('{:}',res1);
            poi2=strfind(md.res2,'}');
            if ~isempty(poi2),
                res2=md.res2(poi1+1:end);
            else
                res2=md.res2;
            end;
            llist{2}=strcat('{:}',res2);
        case 2
            poi1=strfind(md.res1,')');
            if ~isempty(poi1),
                res1=md.res1(poi1+1:end);
            else
                res1=md.res1;
            end;
            lpoi=0;
            chainseq=[md.chainseq md.chainseq(1)];
            for kk=1:length(md.chainseq),
                lpoi=lpoi+1;
                llist{lpoi}=sprintf('(%s)%s',chainseq(kk),res1);
            end;
    end;
    if md.cmode < 3,
        labels=get_labels(llist,testset.library,testset.T,handles);
        lindices=zeros(length(labels),4);
        for kk=1:length(labels),
            cindices=labels(kk).indices;
            if ~isempty(cindices),
                lindices(kk,:)=cindices;
            end;
        end;
    end;
    switch md.cmode
        case {0,-1}
            r=-100;
            sigr=0;
            ind1=resolve_address(md.res1);
            NO1=get_NO(ind1,lindices,labels);
            if NO1,
                NOpos1=labels(NO1).NOpos;
                ind2=resolve_address(md.res2);
                NO2=get_NO(ind2,lindices,labels);
                if NO2,
                    NOpos2=labels(NO2).NOpos;
                end;
                [rax,distr]=get_distribution(NOpos1,NOpos2);
                distr=distr/sum(distr);
                r=sum(rax.*distr);
                diff=rax-r;
                sigr=sqrt(sum(distr.*diff.^2));
            end;
            if r<0,
                fprintf(fid,'ERROR: Distance for pair %s,%s could not be simulated.\n',md.res1,md.res2);
            end;
            md.rsim=r;
            md.ssim=sigr;
        case 1
            ind1=resolve_address(md.res1);
            nn=length(model.structures{ind1(1)}(ind1(2)).xyz);
            poi1=strfind(md.res1,'}');
            if ~isempty(poi1),
                res1=md.res1(poi1+1:end);
            else
                res1=md.res1;
            end;
            poi2=strfind(md.res2,'}');
            if ~isempty(poi2),
                res2=md.res2(poi1+1:end);
            else
                res2=md.res2;
            end;
            rvec=zeros(1,nn);
            svec=zeros(1,nn);
            for kk=1:nn,
                cres1=sprintf('{%i}%s',kk,res1);
                cres2=sprintf('{%i}%s',kk,res2);
                ind1=resolve_address(cres1);
                NO1=get_NO(ind1,lindices,labels);
                if NO1,
                    NOpos1=labels(NO1).NOpos;
                    ind2=resolve_address(cres2);
                    NO2=get_NO(ind2,lindices,labels);
                    if NO2,
                        NOpos2=labels(NO2).NOpos;
                    end;
                    [rax,distr]=get_distribution(NOpos1,NOpos2);
                    distr=distr/sum(distr);
                    r=sum(rax.*distr);
                    diff=rax-r;
                    sigr=sqrt(sum(distr.*diff.^2));
                end;
                if r<0,
                    fprintf(fid,'ERROR: Distance for pair %s,%s could not be simulated.\n',md.res1,md.res2);
                end;
                rvec(kk)=r;
                svec(kk)=sigr;
            end;
            md.rsim=mean(rvec);
            md.ssim=sqrt(sum(svec.^2)/nn);
            sdev=std(rvec);
            fprintf(fid,'Standard deviation over %i models for pair [%s]%s,%s is: %4.2f nm\n',nn,md.pdb,md.res1,md.res2,sdev);
        case 2
            poi1=strfind(md.res1,')');
            if ~isempty(poi1),
                res1=md.res1(poi1+1:end);
            else
                res1=md.res1;
            end;
            lpoi=0;
            chainseq=[md.chainseq md.chainseq(1)];
            rvec=zeros(1,length(md.chainseq));
            svec=zeros(1,length(md.chainseq));
            for kk=1:length(md.chainseq),
                lpoi=lpoi+1;
                cres1=sprintf('(%s)%s',chainseq(kk),res1);
                cres2=sprintf('(%s)%s',chainseq(kk+1),res1);
                ind1=resolve_address(cres1);
                NO1=get_NO(ind1,lindices,labels);
                if NO1,
                    NOpos1=labels(NO1).NOpos;
                    ind2=resolve_address(cres2);
                    NO2=get_NO(ind2,lindices,labels);
                    if NO2,
                        NOpos2=labels(NO2).NOpos;
                    end;
                    [rax,distr]=get_distribution(NOpos1,NOpos2);
                    distr=distr/sum(distr);
                    r=sum(rax.*distr);
                    diff=rax-r;
                    sigr=sqrt(sum(distr.*diff.^2));
                end;
                if r<0,
                    fprintf(fid,'ERROR: Distance for pair %s,%s could not be simulated.\n',md.res1,md.res2);
                end;
                rvec(kk)=r;
                svec(kk)=sigr;
            end;
            md.rsim=mean(rvec);
            md.ssim=sqrt(sum(svec.^2)/length(svec));
            sdev=std(rvec);
            fprintf(fid,'Standard deviation over %i chain pairs for pair [%s]%s,%s is: %4.2f nm\n',length(md.chainseq),md.pdb,md.res1,md.res2,sdev);
        case 3
            [rmean,sigr,rmean_sd,sigr_sd,uc1,uc2] = get_distance_new(md,library,no_context);
            md.rsim=rmean/10;
            md.ssim=sigr/10;
            md.rsim_sd = rmean_sd/10;
            md.ssim_sd = sigr_sd/10;
            md.uc1 = uc1/10;
            md.uc2 = uc2/10;
    end;
    experiment(k)=md.rexp;
    simulation(k)=md.rsim;
    cpoi=cpoi+1;
    cps(spoi)=cpoi;
    exp_ps(spoi,cpoi)=md.rexp;
    sim_ps(spoi,cpoi)=md.rsim;
    if strcmpi(md.sec,'H') && strcmpi(md.exposure,'S'),
        hep=hep+1;
        he_exp(hep)=md.rexp;
        he_sim(hep)=md.rsim;
    end;
    if md.type~=2,
        mp=mp+1;
        m_exp(mp)=md.rexp;
        m_sim(mp)=md.rsim;
    end;
    testset.data(k)=md;
end;

he_exp=he_exp(1:hep);
he_sim=he_sim(1:hep);

m_exp=m_exp(1:mp);
m_sim=m_sim(1:mp);

fprintf(fid,'\n%i test distances were processed\n\n',length(experiment));

[p,s]=polyfit(simulation,experiment,1);
y=polyval(p,[1.5,7.2]);

cc=corrcoef(experiment,simulation);
fprintf(fid,'Correlation coefficient experiment/simulation: %6.4f\n',cc(1,2));
fprintf(fid,'Fit polynomial exp = %5.3f Å + %5.3f sim\n',10*p(2),p(1));

diff=experiment-simulation;
md=mean(diff);
sig=sqrt(sum(diff.^2)/length(diff));
diff2=diff-md;
sig2=sqrt(sum(diff2.^2)/length(diff2));

fprintf(fid,'Mean deviation: %4.1f Å\n',10*md);
fprintf(fid,'Mean absolute deviation: %4.1f Å\n',10*mean(abs(diff)));
fprintf(fid,'Standard deviation: %4.1f Å\n',10*sig);

fprintf(fid,'\nStatistics for %i solvent-exposed helical sites:\n',hep);
[p,s]=polyfit(he_sim,he_exp,1);

cc=corrcoef(he_exp,he_sim);
fprintf(fid,'Correlation coefficient experiment/simulation: %6.4f\n',cc(1,2));
fprintf(fid,'Fit polynomial exp = %5.3f Å + %5.3f sim\n',10*p(2),p(1));

diff=he_exp-he_sim;
md=mean(diff);
sig=sqrt(sum(diff.^2)/length(diff));
diff2=diff-md;
sig2=sqrt(sum(diff2.^2)/length(diff2));

fprintf(fid,'Mean deviation: %4.1f Å\n',10*md);
fprintf(fid,'Mean absolute deviation: %4.1f Å\n',10*mean(abs(diff)));
fprintf(fid,'Standard deviation: %4.1f Å\n',10*sig);

fprintf(fid,'\nStatistics for %i mean distances:\n',mp);
[p,s]=polyfit(m_sim,m_exp,1);

cc=corrcoef(m_exp,m_sim);
fprintf(fid,'Correlation coefficient experiment/simulation: %6.4f\n',cc(1,2));
fprintf(fid,'Fit polynomial exp = %5.3f Å + %5.3f sim\n',10*p(2),p(1));

diff=m_exp-m_sim;
md=mean(diff);
sig=sqrt(sum(diff.^2)/length(diff));
diff2=diff-md;
sig2=sqrt(sum(diff2.^2)/length(diff2));

fprintf(fid,'Mean deviation: %4.1f Å\n',10*md);
fprintf(fid,'Mean absolute deviation: %4.1f Å\n',10*mean(abs(diff)));
fprintf(fid,'Standard deviation: %4.1f Å\n',10*sig);

fprintf(fid,'\nStatistics per structure:\n');

for k=1:spoi,
    nc=cps(k);
    if nc>2,
        m_sim=sim_ps(k,1:nc);
        m_exp=exp_ps(k,1:nc);
        diff=m_exp-m_sim;
        md=mean(diff);
        sig=sqrt(sum(diff.^2)/length(diff));
        fprintf(fid,'%s(%i): sig(dr)= %4.1f, <dsigr>= %4.1f Å, MAE: %4.1f Å\n',stags{k},nc,10*sig,10*md,10*mean(abs(diff)));
    end;
end;

time=toc(tstart);

fprintf(fid,'\nTest computation took %8.1f s\n',time);

fclose(fid);

add_msg_board('Test computation finished');

set(hMain.figure,'Pointer','arrow');

% --------------------------------------------------------------------
function msg=new_pdb(fname)
% load new PDB

global model
global hMain

if exist('model','var')
    model=[];
end;

% initialize display
axes(hMain.axes_model);
cla reset;
axis equal
axis off
set(gca,'Clipping','off');
set(gcf,'Renderer','opengl');
hold on
hMain.camlight=camlight;
guidata(hMain.axes_model,hMain);
hMain.virgin=0;

msg=add_pdb(fname);

function labels=get_labels(llist,library,T,handles,no_context)

global model
global hMain

dynamic=false;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

% check whether sites are already labelled and whether all restraint sites
% do exist, context-free labels are always recomputed
lindices=zeros(length(labels),4);
for k=1:length(labels),
    cindices=labels(k).indices;
    if ~isempty(cindices),
        lindices(k,:)=cindices;
    end;
end;
poi=0;
to_do_list{1}=' ';
for k=1:length(llist),
    adr1=llist{k};
    ind1=resolve_address(adr1);
    if isempty(ind1),
        add_msg_board(sprintf('ERROR: Constraint %i has first label at site %s',k,adr1));
        add_msg_board(sprintf('This site does not exist in current structure %s',mk_address(1)));
    end;
    found=false;
    for l=1:length(labels),
        diff=ind1-lindices(l,:);
        if sum(abs(diff))==0,
            found=true;
        end;
    end;
    if ~found,
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}),
                found=true;
            end;
        end;
        if ~found || no_context,
            poi=poi+1;
            to_do_list{poi}=adr1;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %s %i',to_do_list{k},library,'library',no_context);
        hMain.store_undo=false;
        hMain.dynamic_rotamers=dynamic;
        cmd(handles,command);
    end;
end;

labels=label_information(model.sites);

function labels=label_information(sites)

global model
global label_defs

poi=0;
for k0=1:length(sites),
    for k1=1:length(sites{k0}),
        for k=1:length(sites{k0}(k1).residue),
            poi=poi+1;
            labels(poi).indices=sites{k0}(k1).residue(k).indices;
            id=tag2id(sites{k0}(k1).residue(k).label,label_defs.restags);
            labels(poi).name=label_defs.residues(id).short_name;
            labels(poi).T=sites{k0}(k1).residue(k).T;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            x=sum(NOpos(:,1).*NOpos(:,4));
            y=sum(NOpos(:,2).*NOpos(:,4));
            z=sum(NOpos(:,3).*NOpos(:,4));
            labels(poi).xyz=[x y z];
            labels(poi).rmsd=NOpos_rmsd(NOpos);
            labels(poi).NOpos=NOpos;
        end;
    end;
end;

function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function found=get_NO(ind,lindices,labels)

found=false;
for l=1:length(labels),
    diff=ind-lindices(l,:);
    if sum(abs(diff))==0,
        found=l;
    end;
end;

function [rmean,sigr,rmean_sd,sigr_sd,uc1,uc2] = get_distance_new(md,library,no_context)

rmean_sd = [];
sigr_sd = [];
uc1 = [];
uc2 = [];
poi = strfind(library,'#');
if isempty(poi), mode = 'single'; else mode = 'MC'; end;

rindices1 = resolve_address(md.res1);
rindices2 = resolve_address(md.res2);
switch mode
    case 'single'
        load(library);
        label1 = rotamer_populations(rindices1,rot_lib,-1,0,-1,'','',no_context);
        label2 = rotamer_populations(rindices2,rot_lib,-1,0,-1,'','',no_context);
        if no_context,
            rmean = norm(label1.NOpos(1:3) - label2.NOpos(1:3));
            sigr = 0.05;
        else
            [rax,distr]=get_distribution(label1(1).NOpos,label2(1).NOpos,0.1);
            distr = distr/sum(distr);
            rmean = 10*sum(rax.*distr);
            rdev = 10*rax - rmean;
            sigr = sqrt(sum(rdev.^2.*distr));
        end;
        add_msg_board(sprintf('Mean distance : %5.1f Å',rmean));
        add_msg_board(sprintf('Std. deviation: %5.1f Å',sigr));
    case 'MC'
        poi = strfind(library,'#');
        for k = 1:5,
            clib = sprintf('%s%i%s',library(1:poi-1),k,library(poi+1:end));
            load(clib);
            label1 = rotamer_populations(rindices1,rot_lib,-1,0,-1,'','',no_context);
            label2 = rotamer_populations(rindices2,rot_lib,-1,0,-1,'','',no_context);
            if isempty(label1) || isempty(label2),
                add_msg_board('ERROR: Site is too tight to be labelled with all libraries.');
                return
            end;
            NOpos1{k} = label1.NOpos;
            NOpos2{k} = label2.NOpos;
        end;
        [rax,distr]=get_distribution(NOpos1{1},NOpos2{1},0.1);
        distributions = zeros(25,length(distr));
        rmean = zeros(1,25);
        stddev = zeros(1,25);
        mdistr = zeros(size(distr));
        for k1 = 1:5,
            for k2 = 1:5,
                no = 5*(k1-1)+k2;
                [~,distr]=get_distribution(NOpos1{k1},NOpos2{k2},0.1);
                distr = distr/sum(distr);
                mdistr = mdistr + distr;
                distributions(no,:) = distr;
                rmean(no) = 10*sum(rax.*distr);
                % fprintf(1,'%i: <r> = %8.4f\n',no,rmean(no));
                rdev = 10*rax - rmean(no);
                stddev(no) = sqrt(sum(rdev.^2.*distr));
                % fprintf(1,'%i: sr = %8.4f\n',no,stddev(no));
            end;
        end;
        acoor = zeros(5,3);
        for k1 = 1:5,
            NOpos = NOpos1{k1};
            NO=NOpos(:,1:3); % NO center for k-th frame in the pos. 1
            pop=NOpos(:,4); % weight for k-th frame in the pos. 1    end;
            x=sum(NO(:,1).*pop);
            y=sum(NO(:,2).*pop);
            z=sum(NO(:,3).*pop);
            acoor(k1,:) = [x,y,z];            
        end;
        mcoor = sum(acoor,1)/5;
        dev = acoor - repmat(mcoor,5,1);
        uc1 = sqrt(sum(sum(dev.^2))/5);
        acoor = zeros(5,3);
        for k1 = 1:5,
            NOpos = NOpos2{k1};
            NO=NOpos(:,1:3); % NO center for k-th frame in the pos. 1
            pop=NOpos(:,4); % weight for k-th frame in the pos. 1    end;
            x=sum(NO(:,1).*pop);
            y=sum(NO(:,2).*pop);
            z=sum(NO(:,3).*pop);
            acoor(k1,:) = [x,y,z];            
        end;
        mcoor = sum(acoor,1)/5;
        dev = acoor - repmat(mcoor,5,1);
        uc2 = sqrt(sum(sum(dev.^2))/5);
        mdistr = mdistr/sum(mdistr);
        mrmean = 10*sum(rax.*mdistr);
        mdiff = std(rmean);
        mrdev = 10*rax - mrmean;
        mstddev = sqrt(sum(mrdev.^2.*mdistr));
        msdiff = std(stddev);
        add_msg_board(sprintf('Mean distance : %5.1f Å predicted with uncertainty %5.1f Å',mrmean,mdiff));
        add_msg_board(sprintf('Std. deviation: %5.1f Å predicted with uncertainty %5.1f Å',mstddev,msdiff));
        rmean = mrmean;
        sigr = mstddev;
        rmean_sd = mdiff;
        sigr_sd = msdiff;
end;
