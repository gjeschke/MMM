function display_multilateration

global general
global hMain
global model

snum = model.current_structure;

maxtime = 30*60;
num_models = 200;
grid_points = 176;
grid_size = 75;
opaqueness = 0.5;

node_radius = 2;
constraint_width = 4;
node_color = [0,0,0.25];
cloud_color = [1,0,0;0,1,0];
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

cd(PathName);

[ref_points,dist,sigr,info,cancelled]=process_localization_restraints(restraint_data);

display_mode = restraint_data.locate(1).display_mode;

[np,~] = size(ref_points);

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

if strcmpi(display_mode,'full') || strcmpi(display_mode,'nodes'),
    for k1 = 1:np,
        coor = ref_points(k1,:);
        [x,y,z,t]=point2trisphere(coor,node_radius,2);
        obj=trisurf(t,x,y,z);
        set(obj, 'FaceColor', node_color, 'EdgeColor', 'none', 'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
        set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
    end;
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


% profile on

most_probable = [0,0,0];
add_msg_board('Computing probability density. Please be patient.');
for k=1:m,
    max_probability = 0;
    if k>1,
        tag = sprintf('%s_%i',restraint_data.locate(1).tag,k);
    else
        tag=restraint_data.locate(1).tag;
    end;
    col = cloud_color(k,:);
    cube = zeros(ny,nx,nz);
    x = linspace(point(k,1) - grid_size/2, point(k,1) + grid_size/2, nx);
    y = linspace(point(k,2) - grid_size/2, point(k,2) + grid_size/2, ny);
    z = linspace(point(k,3) - grid_size/2, point(k,3) + grid_size/2, nz);
    fulldistr = 0;
    for kk = 1:length(info),
        fulldistr = fulldistr + info(kk).fulldistr;
    end;
    if fulldistr, % this code duplication saves a lot of time
        for kx = 1:nx,
            for ky = 1:ny,
                for kz = 1:nz,
                    coor = [x(kx),y(ky),z(kz)];
                    diff = ref_points - reprowvector(coor,np); 
                    cdist = sqrt(sum(diff.^2,2));
                    [~,poi] = min(abs(10*info(k).rax - cdist(k)));
                    prob = get_probabilities(cdist,dist,sigr,info);
                    cube(ky,kx,kz) = prod(prob);
                    if cube(ky,kx,kz) > max_probability,
                        max_probability = cube(ky,kx,kz);
                        most_probable = coor;
                    end;
                end;
            end;
        end;
    else
        for kx = 1:nx,
            for ky = 1:ny,
                for kz = 1:nz,
                    coor = [x(kx),y(ky),z(kz)];
                    diff = ref_points - reprowvector(coor,np); 
                    cdist = sqrt(sum(diff.^2,2));
                    earg = (cdist - dist)./(sqrt(2)*sigr);
                    prob = exp(-earg.^2);
                    cube(ky,kx,kz) = prod(prob);
                    if cube(ky,kx,kz) > max_probability,
                        max_probability = cube(ky,kx,kz);
                        most_probable = coor;
                    end;
                end;
            end;
        end;
    end;
    cube = cube/max_probability;
    poi=poi+1;
    model.densities{poi}.x=x;
    model.densities{poi}.y=y;
    model.densities{poi}.z=z;
    model.densities{poi}.tag=tag;
    model.densities{poi}.cube=cube;
    model.densities{poi}.snum=snum;
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
    cnorm = sum(sum(sum(cube)));
    xm=sum(sum(sum(xg.*cube)))/cnorm;
    ym=sum(sum(sum(yg.*cube)))/cnorm;
    zm=sum(sum(sum(zg.*cube)))/cnorm;
    FV = isosurface(xg,yg,zg,cube,level);
    vert = FV.vertices;
    xmin = min(vert(1,:));
    xmax = max(vert(1,:));
    ymin = min(vert(2,:));
    ymax = max(vert(2,:));
    zmin = min(vert(3,:));
    zmax = max(vert(3,:));
    p = patch(FV);
    set(p, 'FaceColor', col, 'EdgeColor', 'none','FaceAlpha',opaqueness,'FaceLighting','gouraud','Clipping','off');
    set(p, 'CDataMapping','direct','AlphaDataMapping','none');
    dg.gobjects=p;
    dg.tag=['density:' tag];
    dg.color=col;
    dg.level=level;
    dg.transparency=opaqueness;
    dg.active=true;
    dg.cog=[xm,ym,zm];
    dg.mlp=most_probable;
    add_msg_board(sprintf('Center of gravity of density cloud %s is (%4.1f, %4.1f, %4.1f) Å',dg.tag,dg.cog(1),dg.cog(2),dg.cog(3)));
    if isfield(model,'surfaces'),
        model.surfaces=[model.surfaces dg];
    else
        model.surfaces=dg;
    end;
    indices = zeros(1,6);
    indices(1) = -length(model.surfaces);
    record_object(p,indices,[(xmin+xmax)/2,(ymin+ymax)/2,(zmin+zmax)/2,xmax-xmin,ymax-ymin,zmax-zmin]);
    if strcmpi(display_mode,'full') || strcmpi(display_mode,'constraints'),
        for k1 = 1:np,
            coor1 = ref_points(k1,:);
            coor2 = most_probable;
            plot3([coor1(1),coor2(1)],[coor1(2),coor2(2)],[coor1(3),coor2(3)],'Color',0.5*cloud_color(k,:),'LineWidth',constraint_width);
        end;
    end;
    n = restraint_data.ensemble;
    if restraint_data.uncertainty > 0.1,
        d = restraint_data.uncertainty;
    else
        d = 2;
    end;
    locations = location_models(x,y,z,cube,n,d);
    if isfield(model,'surfaces'),
        model.surfaces=[model.surfaces dg];
    else
        model.surfaces=dg;
    end;
    for l = 1: length(locations.p),
        [xl,yl,zl,tl]=point2trisphere(locations.xyz(l,:),d/2,2);
        obj=trisurf(tl,xl,yl,zl);
        set(obj, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha',locations.p(l),'FaceLighting','gouraud','Clipping','off');
        set(obj, 'CDataMapping','direct','AlphaDataMapping','none');    
    end;
    camlookat(hMain.axes_model);
    
    c_tag = restraint_data.locate(1).tag;
    if m > 1,
        c_tag = sprintf('%s_%i',c_tag,k);
    end;

    if isfield(restraint_data,'output'),
        wr_locations(c_tag,locations,restraint_data.output);
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

    camlookat(hMain.axes_model);
end;


set(hfig,'Pointer','arrow');

cd(mypath);


function prob = get_probabilities(cdist,dist,sigr,info)

prob = zeros(1,length(cdist));
for k=1:length(cdist),
    if info(k).fulldistr,
        [~,poi] = min(abs(10*info(k).rax - cdist(k)));
        prob(k) = info(k).distr(poi);
    else
        earg = (cdist(k) - dist(k))/(sqrt(2)*sigr(k));
        prob(k) = exp(-earg^2);
    end;
end;


function [ref_points,dist,sigr,info,cancelled]=process_localization_restraints(restraints)

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

ref_tags = ':';
if isfield(restraints,'reference'),
    nr = length(restraints.reference);
    ref_coor = zeros(nr,3);
    for k = 1:nr,
        ref_tags = [ref_tags restraints.reference(k).tag ':'];
        ref_coor(k,:) = restraints.reference(k).xyz;
    end;
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
        id = tag2id(adr1,ref_tags);
        if isempty(id),
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
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.locate),
        adr1=restraints.locate(k).adr;
        id = tag2id(adr1,ref_tags);
        if isempty(id),
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
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels = [];
end;

dist=zeros(length(restraints.locate),1);
sigr=dist;
for k=1:length(restraints.locate),
    adr1=restraints.locate(k).adr;
    id = tag2id(adr1,ref_tags);
    info(k).indices = [];
    dist(k)=10*restraints.locate(k).r;
    sigr(k)=10*restraints.locate(k).sigr;
    info(k).probability=restraints.locate(k).probability;
    info(k).adr=adr1;
    info(k).fulldistr=restraints.locate(k).fulldistr;
    if restraints.locate(k).fulldistr,
        data = load(restraints.locate(k).file);
        info(k).rax=data(:,1);
        info(k).distr=data(:,2);
        info(k).distr = info(k).distr/max(info(k).distr);
    end;
    if isempty(id),
        ind1=resolve_address(adr1);
        info(k).indices=ind1;
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
    else
        info(k).xyz = restraints.reference(id).xyz;
        info(k).rmsd = restraints.reference(id).rmsd;
    end;
end;
ref_points = zeros(length(restraints.locate),3);
for k=1:length(restraints.locate),
    ref_points(k,:) = info(k).xyz;
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

function wr_locations(tag,locations,oname)

% generate header line
header=sprintf('HEADER    LABEL LOCATIONS');
header=fillstr(header,50);
today=date;
today=[today(1:7) today(10:11)];
header=sprintf('%s%s   %s',header,today,tag);
% state supported format and originating program
format=sprintf('REMARK   4 %s COMPLIES WITH FORMAT V. 3.20, 01-DEC-08',tag);
origin=sprintf('REMARK   5 WRITTEN BY MMM (LABEL LOCATIONS)');

if isempty(oname),
    fname = strcat('locations_',tag);
    fname = strcat(fname,'.pdb');
else
    if strfind(oname,'.pdb'),
        fname = oname;
    else
        fname = [oname '.pdb'];
    end;
end;

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
