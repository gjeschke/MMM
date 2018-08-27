function fit_core=generate_constraints(fname,transform,pdbid)
% generate distance constraints for a transform of a ten-helix core
% transporter modelled after core TMD movements in Mhp1
%

ofile=sprintf('%s.dat',fname);
icore=sprintf('%s_Ca_core.mat',pdbid);
load(icore);
load(transform);

fit_core=transform_core(Ca_core,coredef,trafo,TMD_trafos);

dist_list=distance_constraints(pdbid,core,fit_core,transform);

function fit_core=transform_core(Ca_core0,coredef,trafo,TMD_trafos)

midpoint=mean(Ca_core0,1);
[m,n]=size(Ca_core0);
origmat=repmat(midpoint,m,1);
Ca_core=Ca_core0-origmat;

bundle=[1,2,6,7];
hash=[3,4,8,9];
bcoor=extract_domain(bundle,coredef,Ca_core);
hcoor=extract_domain(hash,coredef,Ca_core);
mb=mean(bcoor,1);
mh=mean(hcoor,1);
[pc0,v0]=rmsd_line_3D(bcoor);
d1=norm(pc0+v0-Ca_core(1,:));
d2=norm(pc0-v0-Ca_core(1,:));
if d2>d1,
    v0=-v0;
end;
cconns=mh-mb;
cconns=cconns/norm(cconns);
zs=v0/norm(v0);
y0s=cconns;
xs=cross(y0s,zs);
xs=xs/norm(xs);
ys=cross(zs,xs);
ys=ys/norm(ys);

tmatrix=DCM2affine([xs;ys;zs]);

Ca_core1=affine_coor_set(Ca_core,tmatrix);

fit_core=Ca_core1;

for k=1:length(TMD_trafos),
    [TMD,ilist]=extract_domain(k,coredef,Ca_core1);
    for kk=1:length(TMD_trafos{k}),
        trafonum=TMD_trafos{k}(kk);
        origin=Ca_core1(trafo(trafonum).number,:);
        TMD=shifted_transformation(TMD,trafo(trafonum).tmatrix,origin);
    end;
    fit_core(ilist,:)=TMD;
end;

function [Ca_coor,ilist]=extract_domain(helix_list,coredef,Ca_coor0)

domain_list=coredef(helix_list,:);
[m,n]=size(domain_list);
ilist=[];
for k=1:m,
    ilist=[ilist domain_list(k,1):domain_list(k,2)];
end;
Ca_coor=Ca_coor0(ilist,:);

function dist_list=distance_constraints(pdbid,core,Ca_core,transform)

[mc,nc]=size(core);
% make relative core indices
stdcore=core;
poi=0;
for k=1:mc,
    poi=poi+1;
    stdcore(k,1)=poi;
    len=core(k,2)-core(k,1);
    poi=poi+len;
    stdcore(k,2)=poi;
end;
    
sites=zeros(mc,5);
for k=1:mc,
    currsites(1)=stdcore(k,1);
    currsites(5)=stdcore(k,2);
    currsites(3)=round((currsites(1)+currsites(5))/2); % center of TMD, 1/2 of TMD
    currsites(2)=round((currsites(1)+currsites(3))/2); % at about 1/4 of TMD
    currsites(4)=round((currsites(3)+currsites(5))/2); % at about 3/4 of TMD
    periflag=mod(k,2); % decide whether helix starts at periplasmic or cytoplasmic side
    if ~periflag,
        currsites=fliplr(currsites);
    end;
    sites(k,:)=currsites;
end;

% generate real residue number site table, five sites per TMD, relative residue indices in core

[mc,nc]=size(core);
real_sites=zeros(mc,5);
for k=1:mc,
    currsites(1)=core(k,1);
    currsites(5)=core(k,2);
    currsites(3)=round((currsites(1)+currsites(5))/2); % center of TMD, 1/2 of TMD
    currsites(2)=round((currsites(1)+currsites(3))/2); % at about 1/4 of TMD
    currsites(4)=round((currsites(3)+currsites(5))/2); % at about 3/4 of TMD
    periflag=mod(k,2); % decide whether helix starts at periplasmic or cytoplasmic side
    if ~periflag,
        currsites=fliplr(currsites);
    end;
    real_sites(k,:)=currsites;
end;

cg=sprintf('%s_coarse_grained',pdbid);
load(cg);

[mc,nc]=size(sites);
% determine number of site pairs (distances between "membrane-parallel" sites in TMDs are constrained)
pairs=0;
for k=1:mc-1,
    pairs=pairs+nc*(mc-k);
end;
dist_list=zeros(1,pairs);
residue_list=zeros(pairs,2);

fprintf(1,'Determining distance constraints for structure %s (%i site pairs)\n',pdbid,pairs);

% make distance constraint and residue list
poi=0;
for k=1:mc-1, % TMD1 index
    for kk=k+1:mc, % TMD2 index
        for s=1:5, % site index in TMD
            poi=poi+1;
            res1=sites(k,s);
            residue_list(poi,1)=real_sites(k,s);
            res2=sites(kk,s);
            residue_list(poi,2)=real_sites(kk,s);
            xyz1=Ca_core(res1,:); % Cartesian coordinates of Calpha in TMD1
            xyz2=Ca_core(res2,:); % Cartesian coordinates of Calpha in TMD2
            dist_list(poi)=norm(xyz2-xyz1);
        end;
    end;
end;

% make constraint file for fit from template
fname=sprintf('%s_%s_constraints.dat',pdbid,transform);
fid=fopen(fname,'wt');
fprintf(fid,'%%Direct Calpha constraints for transformation %s\n',transform);
fprintf(fid,'%%of structure %s\n',pdbid);
fprintf(fid,'# basis %i\n',pairs);
fprintf(fid,'# pdb %s\n',pdbid);
fprintf(fid,'# direct\n');
for k=1:pairs,
    fprintf(fid,'%3i   %3i   %5.2f   %5.2f\n',residue_list(k,1),residue_list(k,2),dist_list(k)/10,0.05); % this is in nm!
end;
fprintf(fid,'# end');
fclose(fid);
