function coor = consolidate_coor(coor_cells)
%  coor = consolidate_coor(coor_cells)
%
% Consolidates a cell vector of coordinate arrays to a single coordinate
% array
%
% if one of the arrays does not have dimension (n,3), an empty array is
% returned
%
% G. Jeschke, 24.10.2019

coor = [];

len = length(coor_cells);
nat = 0;
nvec = zeros(1,len);
for k = 1:len
    [n,nd] = size(coor_cells{k});
    if nd ~= 3
        return
    end
    nat = nat + n;
    nvec(k) = n;
end

coor = zeros(nat,3);
poi = 0;
for k = 1:len
    xyz = coor_cells{k};
    coor(poi+1:poi+nvec(k),:) = xyz;
    poi = poi + nvec(k);
end