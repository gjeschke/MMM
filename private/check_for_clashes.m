function clashes = check_for_clashes(coor1,coor2,threshold)

if ~exist('threshold','var') || isempty(threshold)
    threshold = 1.5;
end

clashes.mutual = clash_cost(coor1,coor2,threshold);