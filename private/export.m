function handles=export(handles,format,extension,resolution,option)
% exports the graphics in the detached model window to a file in the
% specified format and with the specified resolution (optional), 
% one further option for Matlab's print command can be specified

global hModel
global hMain
global general

my_path=pwd;
cd(general.reports);

if ~hMain.detached,
    add_msg_board('ERROR: Visualization export requires a detached model window.');
    return;
end;

if nargin<3,
    add_msg_board('ERROR: Not enough input arguments for visualization export.');
    return
end;

if hMain.atom_graphics_auto,
    switch_it=true;
    adjust_atom_graphics(false);
else
    switch_it=false;
end;
[filename, pathname] = uiputfile(extension, 'Export current visualization to graphics file');
if isequal(filename,0) || isequal(pathname,0)
    add_msg_board('Visualization export cancelled by user');
    return
else
    reset_user_paths(pathname);
    general.reports=pathname;
    fname=fullfile(pathname, filename);
    msg=sprintf('Current visualization saved as file: %s',fname);
end;
add_msg_board('Saving visualization to graphics file...');
switch nargin
    case 3
        print(hModel.figure,format,fname);
    case 4
        print(hModel.figure,format,resolution,fname);
    case 5
        print(hModel.figure,format,resolution,option,fname);
end;
add_msg_board(msg);
if switch_it,
    adjust_atom_graphics(true);
end;

cd(my_path);