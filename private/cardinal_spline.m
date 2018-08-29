% Cardinal spline interpolation
%
% Inputs:
%   P   four grid points [P1; P2; P3; P4], array size 4x3
%   c   tensor parameter, c=0 corresponds to Catmull-Rom spline
%   N   number of interpolated points minus one
%
% Outputs:
%   M   interpolated points from P2 to P3, array size (N+1)x3

function M = cardinal_spline(P,c,N)

s = (1-c)./2;

% Construct cardinal spline basis matrix
MC = [-s     2-s   s-2    s;
      2*s   s-3   3-2*s  -s;
      -s     0     s      0;
      0      1     0      0];

u = linspace(0,1,N+1).';
U = [u.^3 u.^2 u ones(size(u))];
M = U*MC*P;
