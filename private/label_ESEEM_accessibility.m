function fname=label_ESEEM_accessibility(label_data)

global hMain

axes(hMain.axes_model);

set(hMain.MMM,'Pointer','watch');

fname=check_labels(label_data);

set(hMain.MMM,'Pointer','arrow');


function fname=check_labels(label_data)

global model
global general
global DPPC

load(label_data.weight_grid_name); % this defines variable: weight_grid, which
                                   % is the relative modulation depth as a
                                   % function of the square sum of grid
                                   % indices kx^2+ky^2+kz^2 assuming a 0.5 Å
                                   % grid, and rgrid, a distance axis
                                   % corresponding to weight_grid
wgl=length(weight_grid);
maxdist=max(rgrid);
blocking_radius=label_data.blocking_radius;
Pi_factor = 10*label_data.Pi_factor;
b=label_data.bulk;

[msg,context]=get_structure(model.current_structure,'xyz');

[blocked,origin]=blocked_cube(context,blocking_radius);

fname0=[general.tmp_files 'label_accessibility.dat'];
fid=fopen(fname0,'w');
fprintf(fid,'%% MMM spin label ESEEM accessibility for %s\n',label_data.name);
fprintf(fid,'%%   z (Å)          population  z offset (Å) eff. acc. vol. %% residue address\n');

fname=[general.tmp_files 'avg_label_accessibility.dat'];
fid1=fopen(fname,'w');
fprintf(fid1,'%% MMM spin label ESEEM accessibility for %s\n',label_data.name);
fprintf(fid1,'%% eff. acc. vol.  mod. depth     Pi      rotamers %% residue address  p(H2O,<z>)  p(H20,z(CA))\n');

cubeax=-maxdist:0.5:maxdist;
cube_num=length(cubeax);
cmid=(cube_num+1)/2;

poi=0;

if isfield(model,'sites'),
    scans=numel(model.sites);
    for k=1:scans,
        chains = numel(model.sites{k});
        for kc = 1:chains,
            if isfield(model.sites{k}(kc),'residue');
                residues=numel(model.sites{k}(kc).residue);
                for kk=1:residues,
                    poi=poi+1;
                    add_msg_board(sprintf('scanning residue %i of %i residues...',kk,residues));
                    numrot=length(model.sites{k}(kc).residue(kk).rotamers);
                    add_msg_board(sprintf('(%i significant rotamers)',numrot));
                    myindices=model.sites{k}(kc).residue(kk).indices;
                    if myindices(1)~=model.current_structure, % skip residues from other structures
                        continue;
                    end;
                    adr=mk_address(myindices);
                    [~,CAcoor] = get_object(strcat(adr,'.CA'),'coor');
                    z_CA = CAcoor(3);
                    eav=zeros(1,numrot);
                    eav2=eav;
                    show_new=0;
                    % if kk==13, show_new=1; end;
                    maxpop=0;
                    popvec=zeros(1,numrot);
                    zvec=zeros(1,numrot);
                    zvec2=zeros(1,numrot);
                    for kkk=1:numrot,
                        popvec(kkk)=model.sites{k}(kc).residue(kk).rotamers(kkk).pop;
                        if model.sites{k}(kc).residue(kk).rotamers(kkk).pop>maxpop,
                            maxpop=model.sites{k}(kc).residue(kk).rotamers(kkk).pop;
                        end;
                    end;
                    popvec=popvec/sum(popvec);
                    for kkk=1:numrot,
                        mean_z=0;
                        zsum=0;
                        [scoor,elements,conn,NO]=get_virtual_label(model.sites{k}(kc).residue(kk),kkk,false);
                        pop=model.sites{k}(kc).residue(kk).rotamers(kkk).pop;
                        N_label=scoor(NO(1),:);
                        O_label=scoor(NO(2),:);
                        NOcoor=(N_label+O_label)/2;
                        if show_new==1,
                            for sss=1:length(elements),
                                plot3(scoor(sss,1),scoor(sss,2),scoor(sss,3),'k.','MarkerSize',18);
                            end;
                            plot3(NOcoor(1),NOcoor(2),NOcoor(3),'r.','MarkerSize',18);
                            [mc,nc]=size(conn);
                            for s1=1:mc,
                                co1=scoor(s1,:);
                                for s2=1:nc,
                                    if conn(s1,s2)~=0,
                                        co2=scoor(conn(s1,s2),:);
                                        plot3([co1(1),co2(1)],[co1(2),co2(2)],[co1(3),co2(3)],'k','LineWidth',2);
                                    end;
                                end;
                            end;
                        end;
                        cube=set_cube(NOcoor,cubeax,origin,blocked);
                        if show_new==1 && kkk==1,
                            p = patch(isosurface(cubeax+NOcoor(1),cubeax+NOcoor(2),cubeax+NOcoor(3),abs(1-cube)));
                            set(p, 'FaceColor', [0.3,0.3,1], 'EdgeColor', 'none','FaceAlpha',0.5,'FaceLighting','gouraud','Clipping','off');
                            set(p, 'CDataMapping','direct','AlphaDataMapping','none');
                        end;
                        [blocked2,origin2,xb2,yb2,zb2]=blocked_cube(scoor,blocking_radius);
                        if show_new==1 && kkk==1,
                            p = patch(isosurface(xb2,yb2,zb2,blocked2));
                            set(p, 'FaceColor', [0,1,0], 'EdgeColor', 'none','FaceAlpha',0.25,'FaceLighting','gouraud','Clipping','off');
                            set(p, 'CDataMapping','direct','AlphaDataMapping','none');
                        end;
                        cube=set_cube(NOcoor,cubeax,origin2,blocked2,cube);
                        dens_sum=0;
                        dens_sum_2=0;
                        pwz_sum = 0;
                        cnorm=0;
                        for kz=1:cube_num,
                            zc=cubeax(kz);
                            pwz=water_population(zc+NOcoor(3),DPPC);
                            kz2=1+(kz-cmid)^2;
                            for ky=1:cube_num,
                                kyz2=kz2+(ky-cmid)^2;
                                for kx=1:cube_num,
                                    kxyz2=kyz2+(kx-cmid)^2;
                                    if kxyz2<=wgl,
                                        dens1=weight_grid(kxyz2);
                                        cnorm=cnorm+dens1;
                                        if ~cube(ky,kx,kz),
                                            dens_sum=dens_sum+dens1;
                                            dens_sum_2=dens_sum_2+pwz*dens1;
                                            cube(ky,kx,kz)=dens1;
                                            zsum=zsum+dens1;
                                            mean_z=mean_z+dens1*cubeax(kz);
                                        else
                                            cube(ky,kx,kz)=0;
                                        end;
                                    else
                                        cube(ky,kx,kz)=0;
                                    end;
                                end;
                            end;
                        end;
                        eav(kkk)=dens_sum/cnorm;
                        eav2(kkk)=dens_sum_2/cnorm;
    %                    eav(kkk)=sum(sum(sum(cube)))/cnorm;
                        if zsum>0, 
                            mean_z=mean_z/zsum; 
                        else
                            mean_z=0;
                        end;
                        zvec(kkk)=NOcoor(3)+mean_z;
                        zvec2(kkk) = NOcoor(3);
                        if show_new==1,
                            v=0.25*max(max(max(cube)));
                            p = patch(isosurface(cubeax+NOcoor(1),cubeax+NOcoor(2),cubeax+NOcoor(3),cube,v));
                            set(p, 'FaceColor', [1,0,0], 'EdgeColor', 'none','FaceAlpha',0.5,'FaceLighting','gouraud','Clipping','off');
                            set(p, 'CDataMapping','direct','AlphaDataMapping','none');
                            show_new=0;
                        end;
                        if isfield(model.info{model.current_structure},'bilayer') && ~isempty(model.info{model.current_structure}.bilayer),
                            my_eav=eav2(kkk);
                        else
                            my_eav=eav(kkk);
                        end;
                        fprintf(fid,'%10.3f%18.6f%10.3f%15.4f      %% %s\n',NOcoor(3),pop,mean_z,my_eav,adr);
                    end;
                    eav_mean=sum(popvec.*eav);
                    eav_mean_2=sum(popvec.*eav2);
                    z_mean = sum(popvec.*zvec2);
                    if isfield(model.info{model.current_structure},'bilayer') && ~isempty(model.info{model.current_structure}.bilayer),
                        fprintf(fid1,'%11.3f%15.4f%10.3f     %4i     %% %s %8.4f  %8.4f\n',eav_mean,eav_mean_2,Pi_factor*eav_mean_2,numrot,adr,water_population(z_mean),water_population(z_CA));
                    else
                        fprintf(fid1,'%11.3f%15.4f%10.3f     %4i     %% %s %8.4f  %8.4f\n',eav_mean,eav_mean_2,Pi_factor*eav_mean,numrot,adr,water_population(z_mean),water_population(z_CA));
                    end;
                end;
            end;
        end;
    end;
end;
fclose(fid);
fclose(fid1);


function [blocked,xyzmin,xb,yb,zb]=blocked_cube(xyz,blocking_radius)
% Given a list xyz of heavy atom coordinates (in Å), the function generates
% a density cube blocked with ones at all points occupied by an atom and
% zeros at all other points
% the cube encompasses all atoms and is defined at a 0.5 Å raster
% its origin xyzmin and axes xb, yb, zb are also provided
%
% this was tested against solvent-accessible surfaces

if nargin<2,
    d0=2.5; % "blocking length" (Å)
else
    d0=blocking_radius;
end;

dk=round(d0/0.5);
dk2=dk^2;

xyzmin=round(2*min(xyz)-dk-2)/2;
xyzmax=round(2*max(xyz)+dk+2)/2;
xb=linspace(xyzmin(1),xyzmax(1),2*(xyzmax(1)-xyzmin(1))+1);
nx=length(xb);
yb=linspace(xyzmin(2),xyzmax(2),2*(xyzmax(2)-xyzmin(2))+1);
ny=length(yb);
zb=linspace(xyzmin(3),xyzmax(3),2*(xyzmax(3)-xyzmin(3))+1);
nz=length(zb);
blocked=zeros(length(yb),length(xb),length(zb),'int8');
[m,n]=size(xyz);
poi=0;
for k=1:m,
    indices=round(2*(xyz(k,:)-xyzmin))+1;
    for kx=indices(1)-dk:indices(1)+dk,
        dx=kx-indices(1);
        if kx>0 && kx<=nx,
            for ky=indices(2)-dk:indices(2)+dk,
                if ky>0 && ky<=ny,
                    dy=ky-indices(2);
                    for kz=indices(3)-dk:indices(3)+dk,
                        if kz>0 && kz<=nz,
                            dz=kz-indices(3);
                            if dx^2+dy^2+dz^2<=dk2,
                                blocked(ky,kx,kz)=1;
                                poi=poi+1;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;


function cube=set_cube(center,cubeax,origin,blocked,cube)
% sets blocked points in a cube
% if argument cube is not provided, it is initialized based on axis length

if nargin<5,
    cube=zeros(length(cubeax),length(cubeax),length(cubeax));
end;
[ix,iy,iz,minshift]=get_indices(center,cubeax,origin,blocked);
if ~isempty(ix) && ~isempty(iy) && ~isempty(iz),
    i0x=minshift(1)+(1:length(ix));
    i0y=minshift(2)+(1:length(iy));
    i0z=minshift(3)+(1:length(iz));
end;
if nargin<5,
    cube(i0y,i0x,i0z)=blocked(iy,ix,iz);
else
    cube(i0y,i0x,i0z)=cube(i0y,i0x,i0z)+double(blocked(iy,ix,iz));
end;

function [ix,iy,iz,minshift]=get_indices(center,cubeax,origin,blocked)
% returns indices into the cube grid blocked that fall into the cube to be
% tested with center point center and axis cubeax

minshift=zeros(1,3);
point0=center+min(cubeax)*ones(1,3);
point1=center+max(cubeax)*ones(1,3);
[ny,nx,nz]=size(blocked);
nmax=[nx,ny,nz];
indices0=round(2*(point0-origin))+1;
indices1=round(2*(point1-origin))+1;

for k=1:3,
    if indices0(k)<1,
        minshift(k)=1-indices0(k);
        indices0(k)=1; 
    end;
    if indices1(k)>nmax(k),
        indices1(k)=nmax(k); 
    end;
end;
ix=indices0(1):indices1(1);
iy=indices0(2):indices1(2);
iz=indices0(3):indices1(3);

function pwz=water_population(zc,lipid)

global DPPC

if nargin<2,
    lipid=DPPC;
end;

poi=1+round((length(lipid.z)-1)*abs(zc)/max(lipid.z));
if poi<1, poi=1; end;
if poi>length(lipid.PW), poi=length(lipid.PW); end;
pwz=lipid.PW(poi);
