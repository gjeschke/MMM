function B = repcolvector(VEC,N)
%REPROWVECTOR Replicate a column vector VEC to give an array with N columns.
%   B = repcolvector(VEC,N) creates a large matrix A consisting of an 
%   N-by-1 tiling of copies of VEC. 
%
%   Based on built-in REPMAT, but faster for this special case
% 

B = VEC(:,ones(1, N));


