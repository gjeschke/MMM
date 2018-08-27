function transmat=view_transform(mode)

global hMain

if nargin<1,
    mode='x';
end;

cam_pos=get(hMain.axes_model,'CameraPosition');
cam_tar=get(hMain.axes_model,'CameraTarget');
view_vec=cam_pos-cam_tar;
view_vec=view_vec/norm(view_vec);

cam_up_vec=get(hMain.axes_model,'CameraUpVector');
cam_up_vec=cam_up_vec/norm(cam_up_vec);

switch mode
    case 'x'
        x=view_vec;
        z=cam_up_vec;
        y=cross_colvec(z,x);
        y=y/norm(y);
        z=cross_colvec(x,y);
        z=z/norm(z);
    case 'y'
        y=view_vec;
        z=cam_up_vec;
        x=cross_colvec(y,z);
        x=x/norm(x);
        z=cross_colvec(x,y);
        z=z/norm(z);
    case 'z'
        z=view_vec;
        y=cam_up_vec;
        x=cross_colvec(y,z);
        x=x/norm(x);
        y=cross_colvec(z,x);
        y=y/norm(y);
end;

transmat=eye(4);
transmat(1,1:3)=x;
transmat(2,1:3)=y;
transmat(3,1:3)=z;

% new_view=affine_trafo_vector(view_vec,transmat);
% new_up=affine_trafo_vector(cam_up_vec,transmat);
