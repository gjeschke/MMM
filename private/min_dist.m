function mdist = min_dist(a,b)
% Returns the minimum distance between a point in coordinate array a and a
% point in coordinate array b
%
% mdist = min_dist(a,b)
%
% G. Jeschke, 21.3.2017

[m1,~] = size(a);
[m2,~] = size(b);
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));
mdist = min(min(pair_dist));