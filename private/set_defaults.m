function handles=set_defaults(handles)
% set default values of (minimal) DEER processing parameters

handles.updated=0;
handles.bckg_dim=3;
set(handles.edit_bckg,'String',sprintf('%4.2f',handles.bckg_dim));
handles.mod_depth=0.4;
set(handles.edit_mod_depth,'String',sprintf('%5.3f',handles.mod_depth));
handles.zero_time=0;
set(handles.edit_zero_time,'String',sprintf('%5i',handles.zero_time));
set(handles.checkbox_form_factor,'Enable','off');
