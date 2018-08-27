function [all_grid,min_entry]=homdimer_grid(constraints,sites,grid) 

% script runs through defined geometrical grid for free 3D rigid-body transformation
% assuming a c2 symmetry of the perfect dimer

% constraints - structure with distance constraints
% sites - structure with simulated sites

% polar angles defining orientation of the dimer c2 axis:
ph=grid.alpha;
th=grid.beta;
% gama=grid.gama; % is set to zero later on when constructing respective rotaion matrix

% set of xyz Cartesian coordinates for translation vector:
x=grid.x;
y=grid.y;
% z=grid.z; % will be set to zero as the translation will take place in the
% dimer plane (perpendicular to the c2 axis)

lph=length(ph);
lth=length(th);
lx=length(x);
ly=length(y);
grid_size=lph*lth*lx*ly;

all_grid0=zeros(grid_size,5);
% all_grid=zeros(length(alpha),length(beta),length(gama),length(x),length(y),length(z)); % initilize grid result variable
% rms_min=10e6;

% keyboard
poi=0; % counter for a number of grids
for k1=1:lph
    for k2=1:lth
        for k3=1:lx
            for k4=1:ly
                poi=poi+1;                  
                v=[ph(k1),th(k2),x(k3),y(k4)];
                 % compute rms
                rms=rms_sites_homdimer(constraints,sites,v);
                all_grid0(poi,:)=[k1 k2 k3 k4 rms];
%                         all_grid(k1,k2,k3,k4)=rms;
%                         if rms<rms_min
%                             rms_min=rms;
%                             min_entry=[k1,k2,k3,k4,rms_min];
%                         end;
                 
            end
        end
    end
end
% keyboard
% [all_grid,ind]=sortrows(all_grid0,7);
% min_entry=all_grid0(ind(1),:);

[all_grid0,~]=sortrows(all_grid0,5);
min_entry=all_grid0(1,:);
if poi>100000
    all_grid=all_grid0(1:100000,:);
else
    all_grid=all_grid0;
end
