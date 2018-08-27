function [vin,handback]=disp_update_MMM(handles,hObject,vin,fstr,intgflag)
% Update an edit field
% with check whether value is in proper range (vmin,vmax) and reset to
% default value vdef if necessary, format string fstr for sprintf
% determines appearance
% if intgflag=1, value is an integer and is rounded
% the corrected value is the output parameter
%

if nargin<5,
    intgflag=0;
end;
if intgflag, vin=round(vin); end;
pstr=sprintf(fstr,vin);
set(hObject,'String',pstr);
% Update handles structure
guidata(hObject, handles);
handback=handles;
