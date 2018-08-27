function bilayer_hole(snum)
% Determines which part of the bilayer plane is taken up by the protein
% this information is added to (or updated) in model.info{snum}.bilayer
%
% snum  structure index
%
% G. Jeschke, 2010

global model

T=293;
forgive=2;
% Carbon, 6
Rmin=forgive*2.275; 
eps=0.020; % CT1 "peptide Alpha carbon"

grid=0.2; % grid spacing (in Å) for computation
rim=5; % rim of the grid outside the protein structure

if ~isfield(model.info{snum},'bilayer') || isempty(model.info{snum}.bilayer),
    return;
end;

add_msg_board('Computing lipid bilayer environment of the structure');
add_msg_board('Please be patient...');
h=gcf;
set(h,'Pointer','watch');
drawnow;

thickness=model.info{snum}.bilayer.thickness;
[msg,xyz]=get_object(snum,'xyz','-nowater');
[msg,elm]=get_object(snum,'elements','-nowater');
[m,n]=size(xyz);
fprintf(1,'%i atoms in structure.\n',m);
poi=find(elm~=1); % exclude protons
xyz=xyz(poi,:);
[m,n]=size(xyz);
fprintf(1,'%i atoms after proton removal.\n',m);
poi=find(abs(xyz(:,3))<=thickness/2); % restrict atoms to bilayer 
xyz=xyz(poi,:);
[m,n]=size(xyz);
fprintf(1,'%i atoms in bilayer range.\n',m);
poi=find(xyz(:,3)>=0); % indices for upper monolayer
uxyz=xyz(poi,:);
[m,n]=size(uxyz);
fprintf(1,'%i atoms in upper monoloayer.\n',m);
poi=find(xyz(:,3)<0); % indices for lower monolayer
lxyz=xyz(poi,:);
[m,n]=size(lxyz);
fprintf(1,'%i atoms in lower monoloayer.\n',m);

minx=min(xyz(:,1));
maxx=max(xyz(:,1));
miny=min(xyz(:,2));
maxy=max(xyz(:,2));

xax=minx-rim:grid:maxx+rim;
nx=length(xax);
yax=miny-rim:grid:maxy+rim;
ny=length(yax);
upper_layer=zeros(nx,ny);
lower_layer=zeros(nx,ny);
fprintf(1,'xy data array has size (%i,%i).\n',nx,ny);

uxyz2=uxyz;
lxyz2=lxyz;
xyz3=uxyz;
dist=uxyz(:,2);
for kx=1:nx,
    x=xax(kx);
    poi=abs(uxyz(:,1)-x)<2*Rmin; % only nearby atoms
    uxyz2=uxyz(poi,:);
    [m,n]=size(uxyz2);
    poi=abs(lxyz(:,1)-x)<2*Rmin; % only nearby atoms
    lxyz2=lxyz(poi,:);
    for ky=1:ny,
        y=yax(ky);
        poi=abs(uxyz2(:,2)-y)<2*Rmin; % only nearby atoms
        xyz3=uxyz2(poi,1:2);
        [m,n]=size(xyz3);
        lipid=repmat([x,y],m,1);
        dist=sqrt(sum((xyz3-lipid).^2,2));
        upper_layer(kx,ky)=repulsion_population(dist,eps,Rmin,T);
        poi=abs(lxyz2(:,2)-y)<2*Rmin; % only nearby atoms
        xyz3=lxyz2(poi,1:2);
        [m,n]=size(xyz3);
        lipid=repmat([x,y],m,1);
        dist=sqrt(sum((xyz3-lipid).^2,2));
        lower_layer(kx,ky)=repulsion_population(dist,eps,Rmin,T);
    end;
end;
figure(1); clf; image(100*(lower_layer>0.02));
figure(2); clf; image(100*(upper_layer>0.02));
set(h,'Pointer','arrow');
fprintf(1,'Now removing inside holes.\n');
% remove inside holes
for kx=1:nx,
    uylow=find(upper_layer(kx,:)<0.02,1);
    if isempty(uylow), 
        uylow=0;
        uyhigh=ny+1;
    else
        uyhigh=find(upper_layer(kx,:)<0.02,1,'last');
    end;
    lylow=find(lower_layer(kx,:)<0.02,1);
    if isempty(lylow), 
        lylow=0; 
        lyhigh=ny+1;
    else
        lyhigh=find(lower_layer(kx,:)<0.02,1,'last');
    end;
    for ky=1:ny,
        if upper_layer(kx,ky)>0.02 && ky>uylow && ky<uyhigh,
            uxlow=find(upper_layer(:,ky)<0.02,1);
            if ~isempty(uxlow),
                uxhigh=find(upper_layer(:,ky)<0.02,1,'last');
                if kx>uxlow && kx<uxhigh,
                    upper_layer(kx,ky)=0;
                end;
            end;
        end;
        if lower_layer(kx,ky)>0.02 && ky>lylow && ky<lyhigh,
            lxlow=find(lower_layer(:,ky)<0.02,1);
            if ~isempty(lxlow),
                lxhigh=find(lower_layer(:,ky)<0.02,1,'last');
                if kx>lxlow && kx<lxhigh,
                    lower_layer(kx,ky)=0;
                end;
            end;
        end;
    end;
end;
figure(3); clf; image(100*(lower_layer>0.02));
figure(4); clf; image(100*(upper_layer>0.02));

model.info{snum}.bilayer.xax=xax;
model.info{snum}.bilayer.yax=yax;
model.info{snum}.bilayer.upper_layer=upper_layer;
model.info{snum}.bilayer.lower_layer=lower_layer;



function pop=repulsion_population(dist,eps,Rmin,T)

gas_un=8.314472;    % universal gas constant in CI (J/(mol*K)    

if isempty(dist),
    pop=1;
    return
end;

arg=Rmin*ones(size(dist))./dist;
energy=eps*sum(arg.^12);
pop=exp(-energy/(gas_un*T));