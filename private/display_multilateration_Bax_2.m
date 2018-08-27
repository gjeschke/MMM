function display_multilateration_Bax_2

global general
global hMain
global model
global Ramachandran
global lovell

load Lovell_Pre_Proline
lovell.prepro = ram_allowed;
load Lovell_Proline
lovell.pro = ram_allowed;
load Lovell_Glycine
lovell.gly =ram_allowed;
load Lovell_general
lovell.general = ram_allowed;

rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs

load([general.Ramachandran 'Ramachandran']);

loc126 = [-10.69,17.12,-4.04];

sequence = 'CTKVPEL';

maxtime = 2*3600;
num_models = 5;
grid_points = 151;
grid_size = 60;
opaqueness = 0.5;
show_locations = true;

locate_H6 = true;

mypath=pwd;

cd(general.restraint_files);

[FileName,PathName] = uigetfile({'*.dat'},'Select constraints file');
if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Loading of constraints canceled by user.');
    message.error=5;
    message.text='Cancelled.';
    return;
end;
reset_user_paths(PathName);
general.restraint_files=PathName;
restraints=fullfile(PathName,FileName);

restraint_data=rd_restraints(restraints);

[ref_points,dist,sigr,info,cancelled]=process_localization_restraints(restraint_data,loc126);

resnum = 0;
switch restraint_data.locate(1).tag,
    case '149R1'
        resnum = 149;
    case '169R1'
        resnum = 169;
    case '193R1'
        resnum = 193;
end;


[np,~] = size(ref_points);

for k = 1:np,
    if strcmpi(info(k).adr,'126'),
        ref_points(k,:) = loc126;
        info(k).xyz = loc126;
    end;
end;

[point, sos, singular] = multilaterate(ref_points,dist);

[m,~]=size(point);
if m==1,
    if singular,
        add_msg_board('Warning: Ill-conditioned matrix encountered in preliminary coordinate estimate');
    end;
    add_msg_board(sprintf('Located mean coordinate: [%5.2f,%5.2f,%5.2f] Å',point));
    add_msg_board(sprintf('Root mean square inconsistency: %5.2f Å',sqrt(sos/np)));    
else
    add_msg_board('Ambiguous localization');
    add_msg_board(sprintf('Located mean coordinate 1: [%5.2f,%5.2f,%5.2f] Å',point(1,:)));
    add_msg_board(sprintf('Located mean coordinate 2: [%5.2f,%5.2f,%5.2f] Å',point(2,:)));
end;

hfig = gcf;
set(hfig,'Pointer','watch');
nx=grid_points;
ny=nx;
nz=nx;
if isfield(model,'density_tags'),
    if isfield(model,'densities'),
        poi=length(model.densities);
        if poi==0;
            model=rmfield(model,'densities');
             model.density_tags=':';
        end;
    end;
else
    poi=0;
    model.density_tags=':';
end;

most_probable = [0,0,0];
max_probability = 0;
add_msg_board('Computing probability density. Please be patient.');
for k=1:m,
    if k>1,
        tag = sprintf('%s_%i',restraint_data.locate(1).tag,k);
    else
        tag=restraint_data.locate(1).tag;
    end;
    col = [1,0,0];
    cube = zeros(ny,nx,nz);
    x = linspace(point(k,1) - grid_size/2, point(k,1) + grid_size/2, nx);
    y = linspace(point(k,2) - grid_size/2, point(k,2) + grid_size/2, ny);
    z = linspace(point(k,3) - grid_size/2, point(k,3) + grid_size/2, ny);
    tic,
    for kx = 1:nx,
        for ky = 1:ny,
            for kz = 1:nz,
                coor = [x(kx),y(ky),z(kz)];
                diff = ref_points - repmat(coor,np,1);
                cdist = sqrt(sum(diff.^2,2));
                earg = (cdist - dist)./sigr;
                prob = exp(-earg.^2);
                cube(ky,kx,kz) = prod(prob);
                if cube(ky,kx,kz) > max_probability,
                    max_probability = cube(ky,kx,kz);
                    most_probable = coor;
                end;
            end;
        end;
    end;
    toc,
    poi=poi+1;
    model.densities{poi}.x=x;
    model.densities{poi}.y=y;
    model.densities{poi}.z=z;
    model.densities{poi}.tag=tag;
    model.densities{poi}.cube=cube;
    model.density_tags=sprintf('%s%s:',model.density_tags,tag);
    sdens=sum(sum(sum(cube)));
    level=restraint_data.locate(1).probability;
    level0 = level;
    for kl=1:99,
        mask=(cube>=kl/100);
        test=sum(sum(sum(mask.*cube)));
        if test<=level0*sdens,
            level=kl/100;
            break;
        end;
    end;
    [xg,yg,zg]=meshgrid(x,y,z);
    p = patch(isosurface(xg,yg,zg,cube,level));
    set(p, 'FaceColor', col, 'EdgeColor', 'none','FaceAlpha',opaqueness,'FaceLighting','gouraud','Clipping','off');
    set(p, 'CDataMapping','direct','AlphaDataMapping','none');
    dg.gobjects=p;
    dg.tag=['LOC_' tag];
    dg.color=col;
    dg.level=level;
    dg.transparency=opaqueness;
    dg.active=true;

    if isfield(model,'surfaces'),
        model.surfaces=[model.surfaces dg];
    else
        model.surfaces=dg;
    end;
    camlookat(hMain.axes_model);
end;

add_msg_board(sprintf('Highest location probability is %5.3f at [%4.1f; %4.1f, %4.1f].',max_probability,most_probable));
m1 = max(max(cube(1,:,:)));
m2 = max(max(cube(end,:,:)));
m3 = max(max(cube(:,1,:)));
m4 = max(max(cube(:,end,:)));
m5 = max(max(cube(:,:,1)));
m6 = max(max(cube(:,:,end)));
mm = max(max(max(cube)));
border_prop = max([m1,m2,m3,m4,m5,m6])/mm;

add_msg_board(sprintf('Maximum relative probability at density cube border is %5.2f%%',100*border_prop));

% ### the following is Bax-specific ###
% it requires that the z axis is the dimer symmetry axis

h6m=0;
h6coor = zeros(2000,3);
for k=131:149,
    adr = sprintf('[HEL6]%i',k);
    [~,rcoor] = get_object(adr,'xyz_paradigm');
    [mc,~] = size(rcoor);
    h6coor(h6m+1:h6m+mc,:) = rcoor;
    h6m=h6m+mc;
end;
h6coor = h6coor(1:h6m,:);

[~,pcoor] = get_object('[BAXC]','xyz');

add_msg_board('Now checking intra-dimer constraint');

most_probable_2 = [0,0,0];
max_probability_2 = 0;

switch resnum
    case 149
        r_intradimer = 23.9; % ### hard-coded constraint ###
        sigr_intradimer = 3.4; % ### hard-coded constraint ###
    case 169
        r_intradimer = 30.0; % ### hard-coded constraint ###
        sigr_intradimer = 3.0; % ### hard-coded constraint ###
    case 193
%         r_intradimer = 21.1; % ### hard-coded constraint ###
%         sigr_intradimer = 4.5; % ### hard-coded constraint ###
        r_intradimer = 44.9; % ### hard-coded constraint ###
        sigr_intradimer = 3.0; % ### hard-coded constraint ###
    otherwise
        fprintf(2,'Residue with unknown intradimer distance constraint. Exiting.\n');
        set(hfig,'Pointer','arrow');
        cd(mypath);
        return
end;

for k=1:m,
    if k>1,
        tag = sprintf('%s_%i_sym',restraint_data.locate(1).tag,k);
    else
        tag=restraint_data.locate(1).tag;
    end;
    col = [0,1,0];
    cube = model.densities{k}.cube; % ### assumes that modeling was started from a "fresh" MMM
    x = linspace(point(k,1) - grid_size/2, point(k,1) + grid_size/2, nx);
    y = linspace(point(k,2) - grid_size/2, point(k,2) + grid_size/2, ny);
    z = linspace(point(k,3) - grid_size/2, point(k,3) + grid_size/2, ny);
    tic,
    for kx = 1:nx,
        for ky = 1:ny,
            for kz = 1:nz,
                coor = [x(kx),y(ky),z(kz)];
                rdimer = 2*sqrt(coor(1)^2+coor(2)^2); % intradimer distance is twice the distance of 149 from the symmetry axis (z)
                earg = (rdimer - r_intradimer)/sigr_intradimer; 
                prob = exp(-earg.^2);
                cube(ky,kx,kz) = prob*cube(ky,kx,kz);
                if cube(ky,kx,kz) > max_probability_2,
                    max_probability_2 = cube(ky,kx,kz);
                    most_probable_2 = coor;
                end;
            end;
        end;
    end;
    toc,
    poi=poi+1;
    model.densities{poi}.x=x;
    model.densities{poi}.y=y;
    model.densities{poi}.z=z;
    model.densities{poi}.tag=tag;
    model.densities{poi}.cube=cube;
    model.density_tags=sprintf('%s%s_SYM:',model.density_tags,tag);
    sdens=sum(sum(sum(cube)));
    level=restraint_data.locate(1).probability;
    level0 = level;
    for kl=1:99,
        mask=(cube>=kl/100);
        test=sum(sum(sum(mask.*cube)));
        if test<=level0*sdens,
            level=kl/100;
            break;
        end;
    end;
    [xg,yg,zg]=meshgrid(x,y,z);
    p = patch(isosurface(xg,yg,zg,cube,level));
    set(p, 'FaceColor', col, 'EdgeColor', 'none','FaceAlpha',opaqueness,'FaceLighting','gouraud','Clipping','off');
    set(p, 'CDataMapping','direct','AlphaDataMapping','none');
    dg.gobjects=p;
    dg.tag=['LOC_' tag];
    dg.color=col;
    dg.level=level;
    dg.transparency=opaqueness;
    dg.active=true;
    
    tic,
    n = 20;
    d = 1.5;
    locations = location_models(x,y,z,cube,n,d);
    time = toc;
    fprintf('Determining %i locations took %5.2f s\n',n,time);
    wr_locations(restraint_data.locate(1).tag,locations);

    if isfield(model,'surfaces'),
        model.surfaces=[model.surfaces dg];
    else
        model.surfaces=dg;
    end;
    if show_locations,
        for l = 1: length(locations.p),
            [xl,yl,zl,tl]=point2trisphere(locations.xyz(l,:),d/2,2);
            obj=trisurf(tl,xl,yl,zl);
            set(obj, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha',locations.p(l),'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none');    
        end;
    end;
    camlookat(hMain.axes_model);
end;

save(dg.tag,'cube','x','y','z');

add_msg_board(sprintf('Highest location probability with intra-dimer constraint is %5.3f at [%4.1f; %4.1f, %4.1f].',max_probability_2,most_probable_2));
add_msg_board(sprintf('%4.1f%% of integral probability within %4.1f%% of maximum probability density.',100*level0,100*level));
add_msg_board(sprintf('This corresponds to an %4.1f Å intra-dimer distance %i/%i.',2*sqrt(most_probable_2(1)^2+most_probable_2(2)^2),resnum,resnum)); 


if resnum == 193,    
    add_msg_board('Computing distance distributions...');
    [~,CB55] = get_object('[BAXC](A)55.CB','coor');
    [~,CB72] = get_object('[BAXC](A)72.CB','coor');
    [~,CB87] = get_object('[BAXC](A)87.CB','coor');
    ra = 10;
    re = 80;
    nr = 651;
    rax = linspace(ra,re,nr);
    distr55 = zeros(1,nr);
    distr72 = zeros(1,nr);
    distr87 = zeros(1,nr);
    distr55m = zeros(1,nr);
    distr72m = zeros(1,nr);
    distr87m = zeros(1,nr);
    nrm = nr -1;
    dr = re -ra;
    for kx = 1:nx,
        for ky = 1:ny,
            for kz = 1:nz,
                coor = [x(kx),y(ky),z(kz)];
                p = cube(ky,kx,kz);
                r55 = norm(coor - CB55);
                n55 = 1 +round(nrm*(r55 - ra)/dr);
                if n55 > 0 && n55 <= nrm
                    distr55(n55) = distr55(n55) + p;
                    distr55m(n55) = distr55m(n55) + mask(ky,kx,kz);
                end;
                r72 = norm(coor - CB72);
                n72 = 1 + round(nrm*(r72 - ra)/dr);
                if n72 > 0 && n72 <= nrm
                    distr72(n72) = distr72(n72) + p;
                    distr72m(n72) = distr72m(n72) + mask(ky,kx,kz);
                end;
                r87 = norm(coor - CB87);
                n87 = 1 + round(nrm*(r87 - ra)/dr);
                if n87 > 0 && n87 <= nrm
                    distr87(n87) = distr87(n87) + p;
                    distr87m(n87) = distr87m(n87) + mask(ky,kx,kz);
                end;
            end;
        end;
    end;
    figure(13); clf;
    plot(rax/10,distr55,'r');
    hold on;
    plot(rax/10,distr72,'b');
    plot(rax/10,distr87,'g');

    figure(14); clf;
    plot(rax/10,distr55m,'r');
    hold on;
    plot(rax/10,distr72m,'b');
    plot(rax/10,distr87m,'g');
end;


if ~strcmpi(restraint_data.locate(1).tag,'149R1') || ~locate_H6,
    set(hfig,'Pointer','arrow');
    cd(mypath);
    return
end;

[~,coorN1]=get_object('[HEL6]149.N1','coor');
[~,coorO1]=get_object('[HEL6]149.O1','coor');
res149t=(coorN1+coorO1)/2;

add_msg_board(sprintf('Template helix 6, residue 149 coordinate at [%4.1f; %4.1f, %4.1f].',res149t));

structures = 0;
runtime=0;
mm = max(max(max(cube)));
level=restraint_data.locate(1).probability;
level0 = level;
for kl=1:99,
    mask=(cube>=kl/100);
    test=sum(sum(sum(mask.*cube/mm)));
    if test<=level0*sdens,
        level=kl/100;
        break;
    end;
end;
% level = 0.065; % ###

too_low = 0;
trials = 0;
c149 = zeros(num_models,3);
c130 = zeros(num_models,3);

[~,coorN] = get_object('[BAXC]125.N','coor');
[~,coorCA] = get_object('[BAXC]125.CA','coor');
[~,coorC] = get_object('[BAXC]125.C','coor');
[~,coorO] = get_object('[BAXC]125.O','coor');
anchorNp = [coorN;coorCA;coorC;coorO];
[~,coorN] = get_object('[BAXC]126.N','coor');
[~,coorCA] = get_object('[BAXC]126.CA','coor');
[~,coorC] = get_object('[BAXC]126.C','coor');
[~,coorO] = get_object('[BAXC]126.O','coor');
anchorN = [coorN;coorCA;coorC;coorO];
[~,coorN] = get_object('[HEL6]131.N','coor');
[~,coorCA] = get_object('[HEL6]131.CA','coor');
[~,coorC] = get_object('[HEL6]131.C','coor');
[~,coorO] = get_object('[HEL6]131.O','coor');
anchorC_0 = [coorN;coorCA;coorC;coorO];

allowed = find(cube/mm > level); % allowed indices into probability grid
fprintf(1,'There are %i grid points with sufficient label 149 probability of %4.1f%%.\n',length(allowed),100*level);
C2_rotmat = affine('rotz',pi);
outside = 0;
clashes = 0;
pclashes = 0;
noloop = 0;
tic;
label_coora = [res149t 1];
label_coora = label_coora';
dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);
best = 0;
while structures < num_models && runtime<maxtime,
    runtime = toc;
    trials = trials + 1;
    [coor,errcode,clashes] = mk_loop(sequence, anchorN, anchorNp, 100);
    if errcode,
        noloop = noloop + 1;
        fprintf(2,'Loop modeling failure.\n');
        continue;
    end;
    anchorC = coor(end-4:end-1,:);
    loopcoor = coor(5:end-5,:);
    [rms,~,transmat]=rmsd_superimpose(anchorC,anchorC_0);
    label_coorc = transmat*label_coora;
    [mx,xp] = min(abs(x- label_coorc(1)));
    [my,yp] = min(abs(y- label_coorc(2)));
    [mz,zp] = min(abs(z- label_coorc(3)));
    if mx > dx || my > dy || dz > dz, % label position is outside density cube
        outside = outside + 1;
        continue
    else
        p = cube(yp,xp,zp)/mm;
        if p > best,
            fprintf(1,'Trial %i, New best probability %8.6f for level %8.6f.\n',trials,p,level);
            best = p;
        end;
        if p < level, % insufficient probability at this grid point
            too_low=too_low+1;
        else
            h6s = h6coor;
            h6c = h6s;
            for k=1:h6m,
                h6s(k,:) = affine_trafo_point(h6s(k,:),transmat);
                h6c(k,:) = h6s(k,:);
                h6c(k,1) = -h6c(k,1);
                h6c(k,2) = -h6c(k,2);
            end;
            loops = loopcoor;
            loopc = loops;
            loopc(:,1) = -loopc(:,1);
            loopc(:,2) = -loopc(:,2);
            clash = check_clash([h6s;loops],[h6c;loopc]);
            if clash,
                clashes = clashes + 1;
            else
                pclash = check_clash(h6s,pcoor);
                if pclash,
                    pclashes = pclashes + 1;
                else
                    structures = structures + 1;
                    fprintf(1,'Structure %i with relative probability %6.4f found.\n',structures,p);
                    c130(structures,:) = coor(6,:);
                    c149(structures,:) = anchorC(2,:);
                    descr(structures).c130 = c130(structures,:);
                    descr(structures).l149 = label_coorc(1:3)';
                    descr(structures).probability = p;
                    all_h6coor{structures} = h6s;
                    all_loopcoor{structures} = loopcoor;
                end;
            end;
        end;
    end;
end;
c130 = c130(1:structures,:);
c149 = c149(1:structures,:);

c130s = c130;
c149s = c149;
for k=1:structures,
    c130s(k,:) = affine_trafo_point(c130(k,:),C2_rotmat);
    c149s(k,:) = affine_trafo_point(c149(k,:),C2_rotmat);
end;

fprintf(1,'Out of %i trials %i were rejected by residue 149 location probability. %5.2f%% were accepted.\n',trials,outside); 
fprintf(1,'Out of %i trials %i were rejected by helix 6/helix 6'' clash. %5.2f%% were accepted.\n',trials-outside,clashes); 
fprintf(1,'Out of %i trials %i were rejected by helix 6/core clash. %5.2f%% were accepted.\n',trials-outside-clashes,pclashes); 
fprintf(1,'Needed %i trials to generate %i structure models in %5.1f s.\n',trials,structures,runtime);



for k=1:structures,
    plot3([c130(k,1),c149(k,1)],[c130(k,2),c149(k,2)],[c130(k,3),c149(k,3)],'b','LineWidth',3);
    plot3(c130(k,1),c130(k,2),c130(k,3),'g.','MarkerSize',14);
    plot3(c149(k,1),c149(k,2),c149(k,3),'r.','MarkerSize',14);
    plot3([c130s(k,1),c149s(k,1)],[c130s(k,2),c149s(k,2)],[c130s(k,3),c149s(k,3)],'m','LineWidth',3);
    plot3(c130s(k,1),c130s(k,2),c130s(k,3),'c.','MarkerSize',14);
    plot3(c149s(k,1),c149s(k,2),c149s(k,3),'y.','MarkerSize',14);
end;

save Bax_H6_NMR_models descr all_h6coor all_loopcoor

% figure(7); clf;
% distances = zeros(1,length(descr));
% probabilities = zeros(1,length(descr));
% d126 = zeros(1,length(descr));
% for k=1:length(descr),
%     distances(k) = descr(k).dist;
%     probabilities(k) = descr(k).probability;
%     d126(k) = descr(k).d126;
% end;
% 
% plot(distances,'k.');
% 
% figure(8); clf;
% plot(probabilities,'r.');
% 
% figure(9); clf;
% plot(d126,'r.');


figure(hfig);


% ### end of Bax-specific code ###


set(hfig,'Pointer','arrow');

cd(mypath);


function [ref_points,dist,sigr,info,cancelled]=process_localization_restraints(restraints,loc126)

global model
global hMain

cancelled=true;
ref_points=[];
dist=[];
sigr=[];
info.probability=0.95;

if ~isfield(restraints,'locate'),
    return;
end;

snum=model.current_structure;

if isfield(restraints,'PDB'),
    if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
        button = questdlg(sprintf('Constraint file specifies template %s, while current template is %s. Do you want to continue?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
        if strcmp(button,'No'),
            cancelled=true;
            return
        end;
    end;
end;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

% check whether sites are already labelled and whether all restraint sites
% do exist
% identity of the label is checked
% labeling temperature is NOT checked
T_list=zeros(1,200);
if ~isempty(labels),
    lindices=zeros(length(labels),4);
    for k=1:length(labels),
        cindices=labels(k).indices;
        if ~isempty(cindices),
            lindices(k,:)=cindices;
        end;
    end;
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.locate),
        adr1=restraints.locate(k).adr;
        ind1=resolve_address(adr1);
        if isempty(ind1),
            add_msg_board(sprintf('ERROR: Constraint %i has label at site %s',k,adr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            add_msg_board('Processing of localization constraints cancelled');
            cancelled=true;
            return;
        end;
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,restraints.locate(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.locate(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}=restraints.locate(k).label;
                T_list(poi)=restraints.locate(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.locate(k).label,adr1));
            end;
        end;
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.locate),
        adr1=restraints.locate(k).adr;
        found=false;
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.locate(k).label),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            label_list{poi}=restraints.locate(k).label;
            T_list(poi)=restraints.locate(k).T;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

labels=label_information(model.sites);

dist=zeros(length(restraints.locate),1);
sigr=dist;
for k=1:length(restraints.locate),
    adr1=restraints.locate(k).adr;
    ind1=resolve_address(adr1);
    dist(k)=10*restraints.locate(k).r;
    sigr(k)=10*restraints.locate(k).sigr;
    info(k).indices=ind1;
    info(k).probability=restraints.locate(k).probability;
    info(k).adr=adr1;
    if strcmpi(adr1,'126'),
        info(k).xyz = loc126;
        info(k).rmsd = 0.3;
    else
        f1=false;
        for l=1:length(labels),
            diff1=ind1-labels(l).indices;
            if sum(abs(diff1))==0,
                f1=true;
                info(k).xyz=labels(l).xyz;
                info(k).rmsd=labels(l).rmsd;
            end;
        end;
        if ~f1,
            add_msg_board('ERROR: Automatic rotamer computation error.');
            add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
            return;
        end;
    end;
end;
ref_points = zeros(length(restraints.locate),3);
for k=1:length(restraints.locate),
    ref_points(k,:) = info(k).xyz;
end;
if strcmpi(restraints.locate(1).tag,'193R1'),
    ref_points(end+1,:) = [-8.2,8.2,-41.1]; % ### hard-coded constraint 149-193
    dist(end+1) = 36.5; % ### hard-coded constraint 149-193; 36.5 or 53.0 Å
    sigr(end+1) = 3.5; % ### hard-coded constraint 149-193
end;
cancelled=false;

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

function clash = check_clash(coor1,coor2)

% clash = check_clash(coor1,coor2)

% checks whether two sets of coordinates produce a mutual clash, based on
% carbon Lennard-Jones potential and a 0.8 forgive factor

global LJ

conv_factor=(4.185*1000); % conversion to SI, - if working per mol

[m1,~] = size(coor1);
[m2,~] = size(coor2);


allowed = 5;
forgive = 0.9; % get softening factor for the LJ potential
full_clash = 0.05; % any pair of atoms closer than 0.05 Å immediately counts as clash

RT = 3*298*8314/2;
max_energy = allowed*RT;

a = coor1;
b = coor2;
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));
if min(abs(pair_dist)) < full_clash,
    clash = true; return;
end

Rmin2_i0 = LJ.Rmin2(6*ones(m1,1));
eps_i = LJ.eps(6*ones(m1,1))*conv_factor;
Rmin2_j0 = LJ.Rmin2(6*ones(m2,1));
eps_j = LJ.eps(6*ones(m2,1))*conv_factor;

eps = sqrt(eps_i(:)*eps_j);
Rmin2_i = Rmin2_i0*forgive; % soften LJ with the forgive factor
Rmin2_j = Rmin2_j0*forgive; % soften LJ with the forgive factor
Rmin = repmat(Rmin2_i(:),1,m2) + repmat(Rmin2_j,m1,1);
            % ATTENTION! Rmin2_i - is actually Rmin_i/2
q = (Rmin./pair_dist).^2;
q = q.*q.*q;
Energy = eps.*(q.^2-2*q);
Energy(pair_dist>10) = 0;
TotalEnergy = sum(sum(Energy));

if TotalEnergy > max_energy,
    clash = true;
else
    clash = false;
end;

function locations = location_models(x,y,z,cube,n,d)
% Represents a probability density cube by n high probability points 
% the points have a minimum distance d from each other
% with this constraint, the next point always corresponds to the highest
% probability point
%
% x,y,x     axes of the probability density cube
% cube      probability density cube
% n         number of points to be generated
% d         minimum distance between generated points
%
% location  structure of high probability locations with fields:
%           .xyz    location coordinate
%           .p      relative probability

nleft = n;

locations.xyz = zeros(n,3);
locations.p = zeros(n,1);

d2 = d^2;
dx = x(2)-x(1);
nx = floor(d/dx);
dy = y(2)-y(1);
ny = floor(d/dy);
dz = z(2)-z(1);
nz = floor(d/dz);
[my,mx,mz] = size(cube);

cube = cube/max(max(max(cube))); % probability normalization

l = 0;
while nleft > 0,
    nleft = nleft - 1;
    l = l+1;
    [p, position] = max(cube(:)); 
    [poiy,poix,poiz] = ind2sub(size(cube),position);
    locations.xyz(l,:) = [x(poix),y(poiy),z(poiz)];
    locations.p(l) = p;
    xa = poix - nx;
    if xa < 1, xa = 1; end;
    xe = poix + nx;
    if xe > mx, xe = mx; end;
    ya = poiy - ny;
    if ya < 1, ya = 1; end;
    ye = poiy + ny;
    if ye > my, ye = my; end;
    za = poiz - nz;
    if za < 1, za = 1; end;
    ze = poiz + nz;
    if ze > mz, ze = mz; end;
    for kx = xa:xe,
        rx2 = ((kx-poix)*dx)^2;
        for ky = ya:ye,
            ry2 = ((ky-poiy)*dy)^2;
            for kz = za:ze,
                rz2 = ((kz-poiz)*dz)^2;
                if rx2 + ry2 + rz2 <= d2, % this point is closer than d to the just selected location
                    cube(ky,kx,kz) = 0; % and is therefore masked
                end;
            end;
        end;
    end;
end;

function wr_locations(tag,locations)

% generate header line
header=sprintf('HEADER    LABEL LOCATIONS');
header=fillstr(header,50);
today=date;
today=[today(1:7) today(10:11)];
header=sprintf('%s%s   %s',header,today,tag);
% state supported format and originating program
format=sprintf('REMARK   4 %s COMPLIES WITH FORMAT V. 3.20, 01-DEC-08',tag);
origin=sprintf('REMARK   5 WRITTEN BY MMM (LABEL LOCATIONS)');

fname = strcat('locations_',tag);
fname = strcat(fname,'.pdb');
fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='File could not be written';
    return;
end;

fprintf(fid,'%s\n',header);
fprintf(fid,'TITLE     BAX SPIN LABEL LOCATIONS %s\n',tag);
fprintf(fid,'REMARK   4\n%s\n',format);
fprintf(fid,'REMARK   5\n%s\n',origin);

for l = 1:length(locations.p),
    atnum=1000; % initialize atom serial number
    if length(locations.p)>1,
        fprintf(fid,'MODEL     %4i\n',l);
    end;
    cid='A';
    rtag='LOC';
    atag='NAL';
    rid=sprintf('%4i ',999);
    tline='HETATM';
    hetflag=1;
    element='NA';
    if length(element)<2 && length(atag)<4,
        atag=[' ' atag];
    end;
    if length(atag)<4,
        atag=fillstr(atag,4);
    end;
    ltag = ' ';
    occupancy = locations.p(l);
    xyz=locations.xyz(l,:);
    Bfactor=0.1;
    atnum=atnum+1;
    fprintf(fid,'%s%5i %4s%s%s %s%4s   %8.3f%8.3f%8.3f',tline,atnum,atag,ltag,rtag,cid,rid,xyz);
    fprintf(fid,'%6.2f%6.2f          %2s\n',occupancy,Bfactor,element);
    if length(locations.p) > 1,
        fprintf(fid,'ENDMDL\n');
    end;
end;
fclose(fid);

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');
% newstring = char(padarray(uint8(string)', newlength-length(string), 32,'post')');
