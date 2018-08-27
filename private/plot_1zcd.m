function handles=plot_1zcd(handles)

global hMain

whole=9:384;

protein=get_protein('1ZCD');

protein=center_protein(protein);

protein.secondary=zeros(1,length(protein.sequence));

protein=set_secondary(protein,12:30,'helix');
protein=set_secondary(protein,35:41,'helix');
protein=set_secondary(protein,44:50,'sheet');
protein=set_secondary(protein,53:58,'sheet');
protein=set_secondary(protein,59:85,'helix');
protein=set_secondary(protein,95:116,'helix');
protein=set_secondary(protein,121:131,'helix');
protein=set_secondary(protein,134:143,'helix');
protein=set_secondary(protein,150:175,'helix');
protein=set_secondary(protein,181:199,'helix');
protein=set_secondary(protein,205:218,'helix');
protein=set_secondary(protein,223:236,'helix');
protein=set_secondary(protein,247:271,'helix');
protein=set_secondary(protein,290:311,'helix');
protein=set_secondary(protein,327:336,'helix');
protein=set_secondary(protein,340:350,'helix');
protein=set_secondary(protein,357:382,'helix');



axes(handles.axes_model);
cla;
set(gca,'Clipping','off');
hold on

segments=plot_protein(protein,whole);

hMain.camlight=camlight;
guidata(handles.axes_model,hMain);


% shading flat
set(gca,'FontSize',14);
axis equal
axis off
view(0,90);
lighting gouraud 
material shiny
set(handles.button_cam_up,'Enable','off');
set(handles.text_view,'String','view 0 90');
%rotate3D(handles.axes_model);

guidata(handles.axes_model,handles);

axes(handles.axes_frame);
cla;
hold on

plot3([0 1],[0 0],[0 0],'r','LineWidth',2); % x axis
plot3([0 0],[0 1],[0 0],'g','LineWidth',2); % y axis
plot3([0 0],[0 0],[0 1],'b','LineWidth',2); % z axis
view(0,90);
axis([-0.1,1.1,-0.1,1.1,-0.1,1.1]);

guidata(handles.axes_frame,handles);