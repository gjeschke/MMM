% Evaluate Cubic Cardinal spline at N+1 values for given four points and tension.
% Uniform parameterization is used.

% fourpoints is a 4xm array of coordinates of the four points, for 3D m=3.
% T is tension.
% N is number of intervals (spline is evaluted at N+1 values).
% slightly modified version for easier callup, G. Jeschke, 2009


function [MatNbyNPlusOne]=cardinal_spline(fourpoints,T,N)

[m0,m]=size(fourpoints);

MatNbyNPlusOne=zeros(N+1,m);

% u vareis b/w 0 and 1.
% at u=0 cardinal spline reduces to P1.
% at u=1 cardinal spline reduces to P2.

u=0;
MatNbyNPlusOne(1,:)=evalcrdnd(fourpoints(1,:),fourpoints(2,:),fourpoints(3,:),fourpoints(4,:),T,u); % MatNbyNPlusOne(:,1)=length(P0)
du=1/N;
for k=1:N
    u=k*du;
      MatNbyNPlusOne(k+1,:)=evalcrdnd(fourpoints(1,:),fourpoints(2,:),fourpoints(3,:),fourpoints(4,:),T,u);
end


% % This program or any other program(s) supplied with it does not provide any
% % warranty direct or implied. This program is free to use/share for
% % non-commerical purpose only. 
% % contact: M A Khan
% % http://www.geocities.com/mak2000sw/
% % Email: khan_goodluck@yahoo.com 


