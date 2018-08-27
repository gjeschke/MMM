function [ncoor, rmsd, transmat, reflect] = fit_partial_network(reference, coor, indices)
% Optimum superposition of a subset of a coordinate array onto reference
% points, allowing for rotation, translation, and reflection
% the complete coordinate set ist then transformed and handed back
%
% reference     [n,3] coordinate array of reference points
% coor          [m,3] full coordinate array (m>=n, otherwise empty output)
% indices       indices of the reference points in coor array, defaults to
%               all points
% 
% ncoor         [m,3] transformed full coordinate array
% rmsd          root mean square deviation of reference point coordinates
% transmat      affine [4,4] transformation matrix
% reflect       flag that indicates reflection
%
% G. Jeschke, 2013

ncoor = [];
rmsd = 1e6;
transmat = [];
reflect = false;

[n,~] = size(reference);
[m,~] = size(coor);

if nargin<3,
    indices = 1:m;
end;

if m<n || min(indices)<1 || max(indices) > m || length(indices)~= n,
    return;
end;

icoor = -coor;
[rms1,~,transmat1] = rmsd_superimpose(reference,coor(indices,:));
[rms2,~,transmat2] = rmsd_superimpose(reference,icoor(indices,:));
if rms1<rms2,
    rmsd = rms1; 
    transmat = transmat1; 
    coor0 = coor;
else
    rmsd = rms2;
    transmat = transmat2; 
    coor0 = icoor;
    reflect = true;
end;
ncoor=affine_trafo_coor(coor0,transmat);