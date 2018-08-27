function color=color_grade_bwr(x)
% relative color grade in a scheme starting blue and ending red, going via
% white
% x     number between 0 (blue) and 1 (red), with 0.5 corresponding to
%       white, values outside the range result in blue (<0) or red (>1)
sat=0.8;

k=32;
n=2*k+1;
colmap=zeros(n,3);
colmap(1:k+1,1)=sat;
colmap(k+1:n,3)=sat;
colmap(1:k+1,2)=linspace(0,sat,k+1);
colmap(1:k+1,3)=colmap(1:k+1,3);
colmap(k+1:n,1)=linspace(sat,0,k+1);
colmap(k+1:n,2)=colmap(k+1:n,1);

poi=1+round(x*2*k);
if poi<1, poi=1; end;
if poi>n, poi=n; end;

color=colmap(poi,:);