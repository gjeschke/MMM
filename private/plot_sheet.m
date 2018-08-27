function sheet=plot_sheet(backbone,rung,normal,col,spr)
% sheet_handles=plot_sheet(backbone,rung,normal,col)
%
% Plots ribbon model for a beta-sheet
% the graphic output is sent to the current axes object
%
% backbone  backbone curve as obtained with get_ribbon
% rung      backbone rungs as obtained with get_ribbon
% normal    normal vectors on ribbon as obtained with get_ribbon
% col       color specifier, defaults to 'g' (green)
% spr       steps per residue (defaults to 10)
%
% sheet     handles to the graphic objects (tube sections) that define the
%           sheet
%
% (c) G. Jeschke, 2009

% Geometric parameters (in Angstroem)
sheet_width=2;
sheet_arrow=1.75;
sheet_height=0.5;
falpha=1; % opaque

if nargin<5,
    spr=10;
end;

if nargin<4,
    col=[0,0,1]; % use blue color if nothing is specified
end;

[m,n]=size(backbone);
sheet=zeros(1,m-1);

for k=1:m,
    s_width=sheet_width;
    if k>=m-spr,
        s_width=sheet_width*sheet_arrow*(m-k)/spr;
    end;
    p1=backbone(k,:)+s_width*rung(k,:)/2-sheet_height*normal(k,:)/2;
    p2=backbone(k,:)+s_width*rung(k,:)/2+sheet_height*normal(k,:)/2;
    p3=backbone(k,:)-s_width*rung(k,:)/2+sheet_height*normal(k,:)/2;
    p4=backbone(k,:)-s_width*rung(k,:)/2-sheet_height*normal(k,:)/2;
    slice2=[p1;p2;p3;p4];
    if k>1,
        [x,y,z]=tube_section(slice1,slice2);
        han=surf(x,y,z);
        set(han,'EdgeColor','none','FaceColor',col,'FaceAlpha',falpha,'FaceLighting','gouraud','Clipping','off');
        sheet(k-1)=han;
    end;
    slice1=slice2;
end;