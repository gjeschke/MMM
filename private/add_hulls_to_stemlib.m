function add_hulls_to_stemlib(libname)

load(libname,'library');
for k = 1:length(library.chains)
    heavy_coor = eliminate_H(library.chains{k}.xyz{1},library.chains{k}.isotopes);
    faces = convhulln(heavy_coor);
    [na,~] = size(faces);
    [kpa,ar] = reducepatch(faces,heavy_coor,na);
    library.hulls(k).vertices = double(ar);
    library.hulls(k).faces = kpa;
end

save(libname,'library');

function coor = eliminate_H(coor0,isotopes)

[m,~] = size(coor0);
coor = zeros(m,3);
poi = 0;
for k = 1:m
    if isotopes(k,1) ~= 1
        poi = poi + 1;
        coor(poi,:) = coor0(k,:);
    end
end
coor = coor(1:poi,:);