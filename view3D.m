function view3D(arg,arg2)
% view3D(arg,arg2)
%
% view3D  Interactively rotate, zoom and pan the view of a 3-D plot
% --------------------------------------------------------------------
%
% VIEW3D ROT turns on mouse-based 3-D rotation
% VIEW3D ZOOM turns on mouse-based 3-D zoom and pan
% VIEW3D PAN turns on mouse-based 3-D pan
% VIEW3D OFF turns it off
%
% VIEW3D(FIG,...) works on the figure FIG
%
% Double click to restore the original view
%
% hit "z" key over the figure to switch from ROT or PAN to ZOOM
% hit "r" key over the figure to switch from ZOOM or PAN to ROT
% hit "p" key over the figure to switch from ZOOM or ROT to PAN
%
% in ROT mode:
% press and hold left mouse button to rotate about screen xy axis
% press and hold middle mouse button to rotate about screen z axis
% press and hold left mouse button while pressing shift key to rotate about
% screen z axis
% in ZOOM mode:
% press and hold left mouse button to zoom in and out
% press and hold middle mouse button to move the plot (pan)
% press and hold left mouse button while pressing shift key to move the
% plot (pan)
% in PAN mode:
% press and hold left mouse button to move the plot (pan)
%
% 
% --------------------------------------------------------------------
% inspired from rotate3d by The MathWorks, Inc.
%
% Torsten Vogel 09.04.1999 
% tv.volke@bmw.de
% tested under Matlab 5.2
% --------------------------------------------------------------------
%
% slightly adapted for MMM and tested under Matlab R14a
% Gunnar Jeschke, 18.04.2009

% ---------------------------------------------- inputs --------------
% ---------------------------------------------- inputs --------------
if nargin == 0
   error('not enough inputs')
elseif nargin == 1
   if ishandle(arg)
      error('not enough inputs')
      return
   else
      switch(lower(arg))
      case 'rot'
         viewact(gcf,'rot')
      case 'zoom'
         viewact(gcf,'zoom')
      case 'pan'
         viewact(gcf,'pan')
      case 'off'
         viewact(gcf,'off')
      case 'down'
         view3DDownFcn
      case 'up'
         view3DUpFcn
      case 'keypress'
         view3DkeypressFcn
      case 'view_xy' % rotate via screen xy axis
         view3DxyFcn
      case 'view_z'  % rotate via screen z axis
         view3DzFcn
      case 'view_zoom' % zoom in and out
         view3DzoomFcn
      case 'view_pan'  % move the plot 
         view3DpanFcn
      otherwise
         error('misspelled command argument')
      end
   end
elseif nargin==2
   if ~ishandle(arg)
      error('bad figure handle')
   end
   switch(lower(arg2))
   case 'rot'
      viewact(arg,'rot')
   case 'zoom'
      viewact(arg,'zoom')
   case 'pan'
      viewact(arg,'pan')
   case 'off'
      viewact(arg,'off')
   otherwise
      error('misspelled command argument')
   end
end


% ---------------------------------------------- activation ----------
function viewact(fig,what)


% de-/activates view3D for the given figure


view3DObj = findobj(allchild(fig),'Tag','view3DObj');

if strcmp(what,'rot')
   adjust_atom_graphics(true);
   if isempty(view3DObj)
      view3DObj = makeview3DObj(fig); %the small text box at the lower left corner
   end
   vdata = get(view3DObj,'UserData');
   vdata.what = 'rot';
   set(view3DObj,'UserData',vdata);
elseif strcmp(what,'zoom')
   adjust_atom_graphics(true);
   if isempty(view3DObj)
      view3DObj = makeview3DObj(fig); %the small text box at the lower left corner
   end
   vdata = get(view3DObj,'UserData');
   vdata.what = 'zoom';
   set(view3DObj,'UserData',vdata);
elseif strcmp(what,'pan')
   adjust_atom_graphics(true);
   if isempty(view3DObj)
      view3DObj = makeview3DObj(fig); %the small text box at the lower left corner
   end
   vdata = get(view3DObj,'UserData');
   vdata.what = 'pan';
   set(view3DObj,'UserData',vdata);
elseif strcmp(what,'off')
   adjust_atom_graphics(false);
   if isempty(view3DObj)
      return
   end
   vdata = get(view3DObj,'UserData');
   uirestore(vdata.uistate);
   set(fig,'KeyPressFcn',vdata.oldkeypressfcn)
   set(fig,'WindowScrollWheelFcn',@view3DScrollFcn);
   delete(view3DObj);
end


% ---------------------------------------------- view3DDownFcn -------
function view3DDownFcn

global hMain

view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
mouseclick = get(gcf,'SelectionType');
if isempty(view3DObj)
   return
end 
vdata = get(view3DObj,'UserData');
vdata.oldunits = get(gcf,'Units');
set(gcf,'Units','pixels');
vdata.old_pt = get(0,'PointerLocation');
%  ----------------- store or restore previous view
ViewData = get(get(hMain.axes_model,'zlabel'),'UserData'); 
if isempty(ViewData)
   ViewData = manageViewData('get_from_axes');
   set(get(hMain.axes_model,'zlabel'),'UserData',ViewData)
end
if strcmp(mouseclick,'open')
   manageViewData('set_axes',ViewData);
   set(gcf,'Units',vdata.oldunits)
   camlight(hMain.camlight);
   return
end
%  ----------------- display text box
fig_color = get(gcf,'Color');
c = sum([.3 .6 .1].*fig_color);
set(vdata.textbox,'BackgroundColor',fig_color);
if(c > .5)
   set(vdata.textbox,'ForegroundColor',[0 0 0]);
else
   set(vdata.textbox,'ForegroundColor',[1 1 1]);
end
%  ----------------- what to do?
if strcmp(vdata.what,'rot')
   set(hMain.popupmenu_view,'Value',1);
   if strcmp(mouseclick,'normal')
      set(vdata.textbox,'string','Screen XY Rotation');
      set(gcf,'WindowButtonMotionFcn','view3D(''view_xy'')');
      set(gcf,'Pointer','custom','pointershapecdata',pointershapes('rot'));
   elseif strcmp(mouseclick,'extend')
      set(vdata.textbox,'string','Screen Z Rotation');
      set(gcf,'WindowButtonMotionFcn','view3D(''view_z'')');
      set(gcf,'Pointer','custom','pointershapecdata',pointershapes('rot'));
   end
elseif strcmp(vdata.what,'zoom')
   if strcmp(mouseclick,'normal')
      set(vdata.textbox,'string','Zoom');
      set(gcf,'WindowButtonMotionFcn','view3D(''view_zoom'')');
      set(gcf,'Pointer','custom','pointershapecdata',pointershapes('zoom'));
   elseif strcmp(mouseclick,'extend')
      set(vdata.textbox,'string','Pan');
      set(gcf,'WindowButtonMotionFcn','view3D(''view_pan'')');
      set(gcf,'Pointer','custom','pointershapecdata',pointershapes('pan'));
   end
elseif strcmp(vdata.what,'pan')
      set(vdata.textbox,'string','Pan');
      set(gcf,'WindowButtonMotionFcn','view3D(''view_pan'')');
      set(gcf,'Pointer','custom','pointershapecdata',pointershapes('pan'));
end
set(view3DObj,'UserData',vdata)
set(vdata.textbox,'visi','on')


% ---------------------------------------------- view3DUpFcn ---------
function view3DUpFcn


view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
if isempty(view3DObj)
   return
end
vdata = get(view3DObj,'UserData');
set(gcf,'WindowButtonMotionFcn','','Units',vdata.oldunits,'pointer','arrow')
set(view3DObj,'visi','off')
update_frame_display;

% ---------------------------------------------- view3DkeypressFcn ---
function view3DkeypressFcn

global hMain

view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
if isempty(view3DObj)
   return
end
vdata = get(view3DObj,'UserData');
currchar = lower(get(gcf,'currentchar'));
if strcmp(currchar,'r')
    vdata.what = 'rot';
    set(hMain.uitoggletool_pan,'State','off');
    set(hMain.uitoggletool_zoom,'State','off');
    set(hMain.uitoggletool_select,'State','off');
    set(hMain.uitoggletool_rotate,'State','on');
elseif strcmp(currchar,'z')
    vdata.what = 'zoom';
    set(hMain.uitoggletool_pan,'State','off');
    set(hMain.uitoggletool_rotate,'State','off');
    set(hMain.uitoggletool_select,'State','off');
    set(hMain.uitoggletool_zoom,'State','on');
elseif strcmp(currchar,'p')
    vdata.what = 'pan';
    set(hMain.uitoggletool_rotate,'State','off');
    set(hMain.uitoggletool_zoom,'State','off');
    set(hMain.uitoggletool_select,'State','off');
    set(hMain.uitoggletool_pan,'State','on');
end
set(view3DObj,'UserData',vdata)


% ---------------------------------------------- view3DxyFcn ---------
function view3DxyFcn

global hMain

view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
vdata = get(view3DObj,'UserData');
new_pt = get(0,'PointerLocation');
old_pt = vdata.old_pt;
dx = (new_pt(1) - old_pt(1))*.5;
dy = (new_pt(2) - old_pt(2))*.5;
direction = [0 0 1];
coordsys  = 'camera';
pos  = get(hMain.axes_model,'cameraposition' );
targ = get(hMain.axes_model,'cameratarget'   );
dar  = get(hMain.axes_model,'dataaspectratio');
up   = get(hMain.axes_model,'cameraupvector' );
[newPos newUp] = camrotate(pos,targ,dar,up,-dx,-dy,coordsys,direction);
set(hMain.axes_model,'cameraposition', newPos, 'cameraupvector', newUp);
vdata.old_pt = new_pt;
set(view3DObj,'UserData',vdata)
if isfield(hMain,'camlight') && ishandle(hMain.camlight),
    camlight(hMain.camlight);
else
    hMain.camlight=camlight;
end;



% ---------------------------------------------- view3DzFcn ----------
function view3DzFcn

global hMain

view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
vdata = get(view3DObj,'UserData');
new_pt = get(0,'PointerLocation');
old_pt = vdata.old_pt;
dy = (new_pt(2) - old_pt(2))*.5;
camroll(hMain.axes_model,-dy)
camlight(hMain.camlight);
vdata.old_pt = new_pt;
set(view3DObj,'UserData',vdata)


% ---------------------------------------------- view3DzoomFcn -------
function view3DzoomFcn

global hMain

view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
vdata = get(view3DObj,'UserData');
new_pt = get(0,'PointerLocation');
old_pt = vdata.old_pt;
dy = (new_pt(2) - old_pt(2))/abs(old_pt(2));
camzoom(hMain.axes_model,1-dy)
camlight(hMain.camlight);
vdata.old_pt = new_pt;
set(view3DObj,'UserData',vdata)


% ---------------------------------------------- view3DpanFcn --------
function view3DpanFcn

global hMain

view3DObj  = findobj(allchild(gcf),'Tag','view3DObj');
vdata = get(view3DObj,'UserData');
new_pt = get(0,'PointerLocation');
old_pt = vdata.old_pt;
dx = (new_pt(1) - old_pt(1))/old_pt(1)*4;
dy = (new_pt(2) - old_pt(2))/old_pt(2)*4;
campan(hMain.axes_model,-dx,-dy,'camera')
camlight(hMain.camlight);
vdata.old_pt = new_pt;
set(view3DObj,'UserData',vdata)




% ---------------------------------------------- make view3DObj ------
function view3DObj = makeview3DObj(fig)

global hMain

% save the previous state of the figure window
vdata.uistate  = uisuspend(fig);
% the data structure
vdata.what     = [];
vdata.olp_pt   = [];
vdata.textbox  = [];
vdata.oldunits = [];
vdata.oldkeypressfcn = get(fig,'KeyPressFcn');
vdata.oldscrollwheelfcn = get(fig,'WindowScrollWheelFcn');
% view3DObj
view3DObj = uicontrol('style','text','parent',fig,'Units','Pixels',... 
                      'Position',[2 2 130 20],'Visible','off', ...
                      'HandleVisibility','off','Tag','view3DObj');
vdata.textbox  = view3DObj;
% store current view
ViewData = manageViewData('get_from_axes');
set(get(hMain.axes_model,'zlabel'),'UserData',ViewData);
% functions
set(fig,'WindowButtonDownFcn','view3D(''down'')');
set(fig,'WindowButtonUpFcn','view3D(''up'')');
set(fig,'WindowButtonMotionFcn','');
set(fig,'ButtonDownFcn','');
set(fig,'KeyPressFcn','view3D(''keypress'')');
set(fig,'WindowScrollWheelFcn',@view3DScrollFcn);


set(view3DObj,'UserData',vdata);
% ---------------------------------------------- manage ViewData -----
function ViewData = manageViewData(how,data)

global hMain

if nargin == 1 ; data = [];end
props = {
   'DataAspectRatio'
   'DataAspectRatioMode'
   'CameraPosition'
   'CameraPositionMode'
   'CameraTarget'
   'CameraTargetMode'
   'CameraUpVector'
   'CameraUpVectorMode'
   'CameraViewAngle'
   'CameraViewAngleMode'
   'PlotBoxAspectRatio'
   'PlotBoxAspectRatioMode'
   'Units'
   'Position'
   'View'
   'Projection'
};
if strcmp(how,'get_from_axes')
   ViewData = get(hMain.axes_model,props);
elseif strcmp(how,'get_stored')
   ViewData = get(get(hMain.axes_model,'zlabel'),'UserData');
elseif strcmp(how,'set_axes')
   set(hMain.axes_model,props,data)
   ViewData = [];
end
% -------------------------------------------------------------------------
% Update frame display
function update_frame_display
global hMain

set_view;
% cam_up=get(hMain.axes_model,'CameraUpVector');
% set(hMain.axes_frame,'CameraUpVector',cam_up);
% cam_pos=get(hMain.axes_model,'CameraPosition');
% set(hMain.axes_frame,'CameraPosition',cam_pos);



% -------------------------------------------------------------------------
% get some pointer shapes
function shape = pointershapes(arg)


if strcmp(arg,'zoom')
% -- zoom
shape=[ 2   2   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
        2   1   1   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
        2   1   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
        2   1   2   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
        2   1   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN  ;
        2   1   2   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN  ;
        2   1   2   1   1   1   1   1   2 NaN NaN NaN   2   2   2   2  ;
        2   1   2   1   1   2   1   1   1   2 NaN   2   1   2   1   2  ;
        2   1   2   1   2 NaN   2   1   1   1   2   1   1   2   1   2  ;
        2   2   2   2 NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
      NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   2   1   2  ;
      NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
      NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   2   1   2  ;
      NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   1   2  ;
      NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   1   1   2  ;
      NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   2   2  ];
elseif strcmp(arg,'pan')
% -- pan
shape=[ NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
        NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
        NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
          2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
          2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
        NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
        NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
        NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
        NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ];
elseif strcmp(arg,'rot')
% -- rot
shape=[ NaN NaN NaN   2   2   2   2   2 NaN   2   2 NaN NaN NaN NaN NaN ;
        NaN NaN NaN   1   1   1   1   1   2   1   1   2 NaN NaN NaN NaN ;
        NaN NaN NaN   2   1   1   1   1   2   1   1   1   2 NaN NaN NaN ;
        NaN NaN   2   1   1   1   1   1   2   2   1   1   1   2 NaN NaN ;
        NaN   2   1   1   1   2   1   1   2 NaN NaN   2   1   1   2 NaN ;
        NaN   2   1   1   2 NaN   2   1   2 NaN NaN   2   1   1   2 NaN ;
          2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
          2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
          2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
          2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
        NaN   2   1   1   2 NaN NaN   2   1   2 NaN   2   1   1   2 NaN ;
        NaN   2   1   1   2 NaN NaN   2   1   1   2   1   1   1   2 NaN ;
        NaN NaN   2   1   1   1   2   2   1   1   1   1   1   2 NaN NaN ;
        NaN NaN NaN   2   1   1   1   2   1   1   1   1   2 NaN NaN NaN ;
        NaN NaN NaN NaN   2   1   1   2   1   1   1   1   1 NaN NaN NaN ;
        NaN NaN NaN NaN NaN   2   2 NaN   2   2   2   2   2 NaN NaN NaN ];


end
