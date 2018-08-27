% Evaluates ND Cubic Cardinal Spline at parameter value u

% INPUT
% P0,P1,P2,P3  are given four points. 
% P1 and P2 are endpoints of curve.
% P0 and P3 are used to calculate the slope of the endpoints (i.e slope of P1 and
% P2).
% T is tension (T=0 for Catmull-Rom type)
% u is parameter at which spline is evaluated


% OUTPUT
% ND-cardinal spline evaluated values at parameter value u

function [Pu] =evalcrdnd(P0,P1,P2,P3,T,u)

Pu=[];

s= (1-T)./2;
% MC is cardinal matrix
MC=[-s     2-s   s-2        s;
    2.*s   s-3   3-(2.*s)   -s;
    -s     0     s          0;
    0      1     0          0];

for i=1:length(P0)
    G(:,i)=[P0(i);   P1(i);   P2(i);   P3(i)];
end

U=[u.^3    u.^2    u    1];

for i=1:length(P0)
    Pu(i)=U*MC*G(:,i);
end

% This program or any other program(s) supplied with it does not provide any
% warranty direct or implied. This program is free to use/share for
% non-commerical purpose only. 
% contact: M A Khan
% Email: khan_goodluck@yahoo.com 
% % http://www.geocities.com/mak2000sw/


