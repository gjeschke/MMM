% function [all_grid,min_grid]=hetdimer_grid(constraints,sites,grid) 
function [grid_select,min_entry]=hetdimer_grid_large(constraints,sites,grid,N)

% script runs through defined geometrical grid for free 3D rigid-body
% transformation

% to be used for the same purpose as hetdimer_grid.m function.
% Used when grid size is far too large such it cannot be stored/processed
% during computation. Therefore not a complete is stored as a result, but
% rather it's part containing grids with only the best rmsd's. Amount of
% grids to be returned can be varied.

% constraints - structure with distance constraints
% sites - structure with simulated sites
% grid - structure specifying grid
% N - number of the best grids wanted

% reduced_grids - structure with N best grids and N second best grid
% solutions; two grids are needed to avoid saturation of the selection
% process during computation (see script)

% three Euiler angles for a free 3D rotation:
alpha=grid.alpha;
beta=grid.beta;
gama=grid.gama;
% set of xyz Cartesian coordinates for translation vector:
x=grid.x;
y=grid.y;
z=grid.z;


% alpha=234;
% beta=37.5;
% gama=42;
% x=32.5;
% y=5;
% z=7;

% keyboard
rms_cut=1.1; % k1*rms_min - initial selection criterion to build best_grid array
best_grid=zeros(2*N,7)+1000;  % initilize best grid array
best_grid0=zeros(N,7)+1000;  % initilize first best grid array

% min_entry=zeros(1,7);    % initilize entry array to story grid with lowest rms
rms_min=10e6;

poi=0; % counter number of good grids

for k1=1:length(alpha)
    for k2=1:length(beta)
        for k3=1:length(gama)
            for k4=1:length(x)
                for k5=1:length(y)
                    for k6=1:length(z)
                                               
                        v=[alpha(k1),beta(k2),gama(k3),x(k4),y(k5),z(k6)];
                        rms=rms_sites_hetdimer(constraints,sites,v);
                        
 %                         %  compute rms-------------------------------------
%                         if rms<rms_min
%                             rms_min=rms;
% %                             min_entry=[k1,k2,k3,k4,k5,k6,rms_min];
%                         end;
%                                               
%                         
%                         if rms<rms_cut*rms_min,
%                             poi=poi+1;
%                             best_grid0(poi,:)=[k1,k2,k3,k4,k5,k6,rms];
%                             
%                         end;
%                         
%                         
% %                         % if array has run out of space, save previous 
% %                         % best grid and start over with best_grid0:
% 
%                         if poi==N
%                             best_grid(1:N,:)=best_grid0;
%                             best_grid=sortrows(best_grid,-7);
%                             best_grid0=zeros(N,7)+1000;
%                             poi=0;
%                         end
                        % -------------------------------------
                        if rms<rms_min,
                            poi=poi+1;
                            best_grid0(poi,:)=[k1,k2,k3,k4,k5,k6,rms];
                            
                        end;
                        
                        
%                         % if array has run out of space, save previous 
%                         % best grid and start over with best_grid0:

                        if poi==N
                            best_grid(1:N,:)=best_grid0;
                            best_grid=sortrows(best_grid,-7);
                            best_grid0=zeros(N,7)+1000;
                            poi=0;
                            rms_min=best_grid(N+1,7);
                        end

                        
                       
                    end
                end
            end
        end
    end
end

% keyboard
best_grid(1:N,:)=best_grid0;
best_grid=sortrows(best_grid,-7);
grid_select=best_grid((N+1):2*N,:);
grid_select=sortrows(grid_select,7); % just to see the best grids on top of the array
min_entry=grid_select(1,:);
