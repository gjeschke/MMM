function [rax,distr] = mk_distance_distribution(adr1,adr2,label)
% [rax,distr] = mk_distance_distribution(adr1,adr2,label)
%
% returns a distance axis and label-to-label distance distribution for two
% sites specified by residue addresses adr1 and adr2
%
% G. Jeschke, 31.10.2018

global rotamer_libraries

rax = [];
distr = [];

rindices1 = resolve_address(adr1);
[m1,n1] = size(rindices1);
if m1 ~= 1 || n1 ~=4
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr1));
    return;
end
rindices2 = resolve_address(adr2);
[m2,n2] = size(rindices2);
if m2 ~= 1 || n2 ~=4
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr2));
    return
end

rot_lib_name = '';

for k = 1:length(rotamer_libraries)
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label)
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T)
            if rotamer_libraries(k).T(kk) == 298
                id = kk;
            end
        end
        if id >0
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
        end
    end
end

load(rot_lib_name);
label1 = rotamer_populations(rindices1,rot_lib);
label2 = rotamer_populations(rindices2,rot_lib);
if ~isempty(label1) && ~isempty(label2)
    [rax,distr]=get_distribution(label1(1).NOpos,label2(1).NOpos,0.1);
    distr = distr/sum(distr);
else
    rax = [];
    distr = [];
end


