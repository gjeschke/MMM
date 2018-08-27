function reset_user_paths(pname)
% Resets all user directory pathnames to the input path pname if the
% function is called for the first time after program start
%
% G. Jeschke, 2012

global general

if general.virgin_path,
    general.virgin_path=false;
    general.exp_files=pname;
    general.pdb_files=pname;
    general.DEER_files=pname;
    general.reports=pname;
    general.restraint_files=pname;
    general.biblio_files=pname;
    general.reports=pname;
end;