function [v,handback]=edit_update_MMM(handles,hObject,vmin,vmax,vdef,fstr,intgflag)
% Update an edit field
% with check whether value is in proper range (vmin,vmax) and reset to
% default value vdef if necessary, format string fstr for sprintf
% determines appearance
% if intgflag=1, value is an integer and is rounded
% the corrected value is the output parameter
%

if nargin<7
    intgflag=0;
end
v=str2double(get(hObject,'String'));
if intgflag, v=round(v); end
% Protect against wrong inputs
if isnan(v)
    v=vdef;
end
if v<vmin 
    v=vmin;
end
if v>vmax 
    v=vmax;
end
pstr=sprintf(fstr,v);
set(hObject,'String',pstr);
% Update handles structure
guidata(hObject, handles);
handback=handles;
