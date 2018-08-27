function display_angle
% function display_angle
%
% Displays angle between the mean atom coordinates of three objects that
% are currently selected
% if more or less objects are selected, an error message is displayed
%
% G.Jeschke, 2009

global hMain

handles=guidata(hMain.figure);

indices=resolve_address('*');
[m,n]=size(indices);
if m~=3,
    add_msg_board(sprintf('Angle defined only for three selected objects, but %i objects are selected.',m));
else
    indy1=indices(1,:);
    indy1=indy1(indy1>0);
    [msg,xyz1]=get_object(indy1,'xyz');
    xyz1=mean(xyz1,1);
    indy2=indices(2,:);
    indy2=indy2(indy2>0);
    [msg,xyz2]=get_object(indy2,'xyz');
    xyz2=mean(xyz2,1);
    indy3=indices(3,:);
    indy3=indy3(indy3>0);
    [msg,xyz3]=get_object(indy3,'xyz');
    xyz3=mean(xyz3,1);
    if length(indy1)<5 || length(indy2)<5,
        add_msg_board('Angle between mean coordinates of');
    else
        add_msg_board('Angle between');
    end;
    add_msg_board(sprintf('%s,',description(indices(1,:))));
    add_msg_board(sprintf('%s, and',description(indices(2,:))));
    add_msg_board(sprintf('%s',description(indices(3,:))));
    add_msg_board(sprintf('is %5.2f°',180*angle(xyz1,xyz2,xyz3)/pi));
end;
