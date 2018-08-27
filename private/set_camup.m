function set_camup(vec)
% Sets the camera up vector along the specified vector

global hMain

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
vec0=cam_pos-cam_tar;
vec0=vec0/norm(vec0);

vec=[vec0 vec/norm(vec)];
set_view(vec);
