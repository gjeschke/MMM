function tinker_path = get_tinker_path
% Path for Tinker input/output files is directory above the Tinker force
% fieled files

tinker_path = which('amber99.prm'); 
if isempty(tinker_path)
    return
end
tinker_path = fileparts(tinker_path);
poi = strfind(tinker_path,'params');
tinker_path = tinker_path(1:poi-1);
