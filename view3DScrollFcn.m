% ---------------------------------------------- view3DzoomFcn -------
function view3DScrollFcn(hObject, eventdata, handles)

global hMain
global graph_settings

dy=eventdata.VerticalScrollAmount*eventdata.VerticalScrollCount*graph_settings.scroll_speed;;
camzoom(hMain.axes_model,1-dy)
if isfield(hMain,'camlight'),
    camlight(hMain.camlight);
end;
