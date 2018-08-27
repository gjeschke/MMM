function [library,list] = cull_RNA_library(library,threshold)

l0 = length(library.chains);
rmsd_mat = zeros(length(library.chains));
list = zeros(1,length(library.chains));

for k1 = 1:length(library.chains)-1
    coor1 = library.chains{k1}.xyz{1};
    [m,~] = size(coor1);
    for k2 = k1+1:length(library.chains)
        coor2 = library.chains{k2}.xyz{1};
        rmsd = sqrt(sum(sum((coor1-coor2).^2))/m);
        rmsd_mat(k1,k2) = rmsd;
        rmsd_mat(k2,k1) = rmsd;
    end
end

model_order = zeros(1,length(library.chains));
for k = 1:length(library.chains)
    model_order(k) = get_central_model(rmsd_mat,threshold,model_order);
end
poi = 1;
list(1) = model_order(1);
    
for k = 1:length(library.chains)    
    clist = list(list>0);
    agree = rmsd_mat(model_order(k),clist);
    if min(agree) > threshold
        poi = poi + 1;
        list(poi) = model_order(k);
    end
end

list = list(1:poi);
library.chains = library.chains(list);

fprintf(1,'Library culled from %i entries to %i entries at threshold %5.1f Å\n',l0,poi,threshold)

disp('Aber hallo!');

function cmodel = get_central_model(rmsd_mat,threshold,list)

[m,~] = size(rmsd_mat);
max_matches = -1;
for k = 1:m
    if min(abs(list-k)) > 0
        rmsd_vec = rmsd_mat(k,:);
        for kl = 1:length(list)
            if list(kl) > 0
                rmsd_vec(list(kl)) = 1e6;
            end
        end
        matches = sum(rmsd_vec<threshold);
        if matches > max_matches
            cmodel = k;
            max_matches = matches;
        end
    end
end