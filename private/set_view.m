function set_view(vec,add_light)
% Sets the 3D view direction along the specified vector

global hMain

if nargin<2
    add_light=true;
end

if nargin<1 || isempty(vec)
    cam_pos=get(hMain.axes_model,'CameraPosition');
    cam_tar=get(hMain.axes_model,'CameraTarget');
    vec=cam_pos-cam_tar;
    vec=vec/norm(vec);
elseif length(vec)==3
    vec=vec/norm(vec);
    view(hMain.axes_model,vec);
elseif length(vec)==6
    camup=vec(4:6);
    camup=camup/norm(camup);
    vec=vec(1:3);
    vec=vec/norm(vec);
    view(hMain.axes_model,vec);
    set(hMain.axes_model,'CameraUpVector',camup);
end

mode=1;
if norm(vec-[1,0,0])<1e-3, mode=3; end
if norm(vec-[-1,0,0])<1e-3, mode=4; end
if norm(vec-[0,1,0])<1e-3, mode=5; end
if norm(vec-[0,-1,0])<1e-3, mode=6; end
if norm(vec-[0,0,1])<1e-3, mode=7; end
if norm(vec-[0,0,-1])<1e-3, mode=8; end

if mode>2
    set(hMain.popupmenu_view,'Value',mode);
end

set(hMain.text_view,'String',sprintf('%4.3f,%4.3f,%4.3f',vec));

cam_up=get(hMain.axes_model,'CameraUpVector');
set(hMain.axes_frame,'CameraUpVector',cam_up);
set(hMain.text_camup,'String',sprintf('%4.3f,%4.3f,%4.3f',cam_up));
cam_pos=get(hMain.axes_model,'CameraPosition');
set(hMain.axes_frame,'CameraPosition',cam_pos);
if isfield(hMain,'camlight') && ishandle(hMain.camlight)
    camlight(hMain.camlight);
elseif add_light
    hMain.camlight=camlight;
else
    hMain.camlight=light('Parent',hMain.axes_model);
    camlight(hMain.camlight);
end

if isfield(hMain,'ANM_plot')
    if hMain.ANM_plot
        cam_tar=get(hMain.axes_model,'CameraTarget');
        set(hMain.ANM_axes,'CameraTarget',cam_tar);
        set(hMain.ANM_axes,'CameraPosition',cam_pos);
        set(hMain.ANM_axes,'CameraUpVector',cam_up);
        camlookat(hMain.ANM_axes);
    end
end

if isfield(hMain,'fit_axes')
    if hMain.fit_plot
        if isfield(hMain,'fit_axes') && ishandle(hMain.fit_axes)
            cam_tar=get(hMain.axes_model,'CameraTarget');
            set(hMain.fit_axes,'CameraTarget',cam_tar);
            set(hMain.fit_axes,'CameraPosition',cam_pos);
            set(hMain.fit_axes,'CameraUpVector',cam_up);
            camlookat(hMain.fit_axes);
        else
            hMain.fit_plot=false;
        end
    end
end