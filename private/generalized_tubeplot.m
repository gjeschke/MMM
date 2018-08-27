function [x,y,z]=generalized_tubeplot(backbone,rung,normal,sec,spr)
% Usage: [x,y,z]=generalized_tubeplot(backbone,rung,normal,sec)
% 
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% backbone: [3,N] vector of backbone curve data
% rung    : [3,N] normalized rung vectors
% normal  : [3,N] normal vectors on rungs (normalized)
% sec     : secondary structure type, 0 coil, 1 helix, 2 sheet, 3 DNA/RNA
% spr       segments per residue, optional, if present and sec==2, an arrow
%           is modelled at the C-terminus
%
% The algorithms fails if you have bends beyond 90 degrees.
%
% based on tubeplot by
% Janus H. Wesenberg, july 2004

global graph_settings

% Geometric parameters (in Angstroem)
switch sec
    case 0
        width=graph_settings.loop_radius/2;
        height=width;
    case 1
        width=graph_settings.helix_width/2;
        height=graph_settings.helix_height/2;
    case 2
        width=graph_settings.sheet_width/2;
        height=graph_settings.sheet_height/2;
    case 3
        width=graph_settings.NA_width/2;
        height=graph_settings.NA_height/2;
end;

width0=width;

if sec~=2 || nargin<5,
    arrow=0;
else
    arrow=1;
end;

n=8;

ct=0.5*graph_settings.coil_radius;

  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(backbone,2)-1)
    if norm(backbone(:,k)-backbone(:,npoints))>ct;
      npoints=npoints+1;
      backbone(:,npoints)=backbone(:,k);
    end
  end
  %Always include endpoint
  if norm(backbone(:,end)-backbone(:,npoints))>0
    npoints=npoints+1;
    backbone(:,npoints)=backbone(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=backbone(:,[2:end,end])-backbone(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=zeros([3,n+1,npoints+2]);
  
  %precalculate cos and sin factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
      if arrow && k>=npoints-spr,
        width=width0*graph_settings.sheet_arrow*(npoints-k)/(spr-1);
        if width<1e-4*width0,
            width=0.01*width0;
        end;
      end;
    convec=cross_colvec(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross_colvec(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
%     xyz(:,:,k+1)=repmat(backbone(:,k),[1,n+1])+...
%         cfact.*repmat(r*nvec,[1,n+1])...
%         +sfact.*repmat(r*convec,[1,n+1]);
    xyz(:,:,k+1)=repmat(backbone(:,k),[1,n+1])+...
        cfact.*repmat(width*rung(:,k),[1,n+1])...
        +sfact.*repmat(height*normal(:,k),[1,n+1]);
  end;
  
  %finally, cap the ends:
  xyz(:,:,1)=repmat(backbone(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(backbone(:,end),[1,n+1]);
  
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  
  %... and plot:
  if nargout<3, surf(x,y,z); end;
  
  