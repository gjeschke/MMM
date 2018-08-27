function rms=rms_sites_hetdimer(constraints,sites,v) 

% Compute r.m.s. deviation of site-to-site distances from given constraints
% version for a heterodimer dimer

% v             parameter vector, angles alpha v(1), beta v(2) and gama v(3)
%               (in degrees) Euler rotation of the (second) monomer molecule w.r.t. the
%               fixed (first) monomer molecule; 
%               x v(3), y v(4) and z v(3) (in Angstroem) determine translation of 
%               second molecule w.r.t. first molecule
%
% sites         coordinates of sites
% constrainst   structure with constraints

alpha=v(1);
beta=v(2);
gama=v(3);
x=v(4);
y=v(5);
z=v(6);
 
sites1=sites.NOsite1;
sites2=sites.NOsite2;
pairs_exp=constraints.pairs_exp;
sig=pairs_exp(:,4);
[np,~]=size(pairs_exp);
    
% keyboard-----------------------------------------------------------------
% eulang=[alpha,beta,gama]; % Euler angles to rotate molecule 2
% Rp2=get_rotmat(eulang); % rotation matrix to rotate molecule 2
% t=[x,y,z]; % translation vector to translate molecule 2
% sites2=rot_trans(sites2,Rp2,t);
% rms=0; % current r.m.s. deviation
% for k=1:np
%     co1=sites1(k,:); % coordinates in molecule 1
%     co2=sites2(k,:); % coordinates in molecule 2
%     dist=norm(co1-co2); % distance between labelled sites
% %     if constraints.type(k)==1, % actual test of constraints
% 
% %     rms=rms+(dist-pairs_exp(k,3))^2; % add square deviation (no sigma for exp distribution is considered!)
%     sig0=pairs_exp(k,4); % include sigma of the experimental distance distribution (better weighting of the exp constraints)
%     rms=rms+((dist-pairs_exp(k,3))/sig0)^2; % add square deviation
%     
% %     end;
% %     disp(dist);
% % keyboard
% end;
% rms=sqrt(rms/np); % compute current r.m.s.

%-- no real speedup--------------------------------------------------------
% eulang=[alpha,beta,gama]; % Euler angles to rotate molecule 2
% Rp2=get_rotmat(eulang); % rotation matrix to rotate molecule 2
% t=[x,y,z]'; % translation vector to translate molecule 2
% % sites2=rot_trans(sites2,Rp2,t);
% rms=0; % current r.m.s. deviation
% for k=1:np
%     co1=sites1(k,:); % coordinates in molecule 1
%     co2=sites2(k,:); % coordinates in molecule 2
%     co2=(Rp2*co2'+t)';
%     
%     dist=norm(co1-co2); % distance between labelled sites
% %     if constraints.type(k)==1, % actual test of constraints
% 
% %     rms=rms+(dist-pairs_exp(k,3))^2; % add square deviation (no sigma for exp distribution is considered!)
%     sig0=pairs_exp(k,4); % include sigma of the experimental distance distribution (better weighting of the exp constraints)
%     rms=rms+((dist-pairs_exp(k,3))/sig0)^2; % add square deviation
%     
% %     end;
% %     disp(dist);
% % keyboard
% end;
% rms=sqrt(rms/np); % compute current r.m.s.


% keyboard
%- up to 20% speed-up------------------------------------------------------

eulang=[alpha,beta,gama]; % Euler angles to rotate molecule 2
Rp2=get_rotmat(eulang); % rotation matrix to rotate molecule 2
t=[x,y,z]; % translation vector to translate molecule 2
% sites2=rot_trans(sites2,Rp2,t);
sites2=rot_trans_vectorized(sites2,Rp2,t);

% % keep it for future cases for speed up (here it does not bring much)
% [m,~]=size(sites2);
% sites2=sites2';
% t=t';
% tx = t(:,ones(1, m));
% sites2=rot_trans_fast(sites2,Rp2,tx);
% sites2=sites2';

% rms=0; % current r.m.s. deviation

% keyboard
diff=sites1-sites2;
dist=sqrt(diff(:,1).*diff(:,1)+diff(:,2).*diff(:,2)+diff(:,3).*diff(:,3));
var=((dist-pairs_exp(:,3))./sig).*((dist-pairs_exp(:,3))./sig);
% var=((dist-pairs_exp(:,3))).*((dist-pairs_exp(:,3)));

rms=sqrt(sum(var/np));
%--------------------------------------------------------------------------

