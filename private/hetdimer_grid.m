% function [all_grid,min_grid]=hetdimer_grid(constraints,sites,grid) 
% function [all_grid,poi,min_entry]=hetdimer_grid(constraints,sites,grid) % disabled 200613

function [all_grid,min_entry]=hetdimer_grid(constraints,sites,grid) 

% script runs through defined geometrical grid for free 3D rigid-body transformation

% constraints - structure with distance constraints
% sites - structure with simulated sites


% keyboard
% three Euiler angles for a free 3D rotation:
alpha=grid.alpha;
beta=grid.beta;
gama=grid.gama;
% set of xyz Cartesian coordinates for translation vector:
x=grid.x;
y=grid.y;
z=grid.z;

la=length(grid.alpha);
lb=length(grid.beta);
lg=length(grid.gama);
lx=length(grid.x);
ly=length(grid.y);
lz=length(grid.z);
grid_size=la*lb*lg*lx*ly*lz;

all_grid0=zeros(grid_size,7);
% all_grid=zeros(length(alpha),length(beta),length(gama),length(x),length(y),length(z)); % initilize grid result variable
% rms_min=10e6;

% keyboard
poi=0; % counter for a number of grids
for k1=1:length(alpha)
    for k2=1:length(beta)
        for k3=1:length(gama)
            for k4=1:length(x)
                for k5=1:length(y)
                    for k6=1:length(z)
                        poi=poi+1;                  
                        v=[alpha(k1),beta(k2),gama(k3),x(k4),y(k5),z(k6)];
                         % compute rms
                        rms=rms_sites_hetdimer(constraints,sites,v);
                        all_grid0(poi,:)=[k1 k2 k3 k4 k5 k6 rms];
%                         all_grid(k1,k2,k3,k4,k5,k6)=rms;
%                         if rms<rms_min
%                             rms_min=rms;
%                             min_entry=[k1,k2,k3,k4,k5,k6,rms_min];
%                         end;
                    end
                end
            end
        end
    end
end
% keyboard
% [all_grid,ind]=sortrows(all_grid0,7);
% min_entry=all_grid0(ind(1),:);

[all_grid0,~]=sortrows(all_grid0,7);
min_entry=all_grid0(1,:);
if poi>100000
    all_grid=all_grid0(1:100000,:);
else
    all_grid=all_grid0;
end
