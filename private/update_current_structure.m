function update_current_structure
% updates the user interface if the current structure has changed

global model
global hMain

snum=model.current_structure;
stag=model.info{snum}.idCode;
model.current_chain=id2tag(1,model.chain_tags{snum});

if hMain.hierarchy_display,
%     if hMain.large,
%         close(hMain.hierarchy_window_large);
%     else
%         close(hMain.hierarchy_window);
%     end;
    close(hMain.hierarchy_window);
    disp_hierarchy=true;
else
    disp_hierarchy=false;
end;

if isfield(model.info{snum},'resolution') && ~isempty(model.info{snum}.resolution),
    resstring=sprintf('%4.2f Å',model.info{snum}.resolution);
else
    resstring='not specified';
end;
    
set(hMain.MMM,'Name',sprintf('MMM - [%s](%s) Resolution %s',stag,model.current_chain,resstring));
set(hMain.menu_file_save_PDB,'Enable','on');

if disp_hierarchy,
%     if hMain.large,
%         hierarchy_window_large;
%     else
%         hierarchy_window;
%     end;
    hierarchy_window;
end;