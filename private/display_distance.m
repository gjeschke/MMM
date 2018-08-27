function display_distance
% function display_distance
%
% Displays distance between the mean atom coordinates of two objects that
% are currently selected
% if more or less objects are selected, an error message is displayed
%
% G.Jeschke, 2009

global hMain

handles=guidata(hMain.figure);

indices=resolve_address('*');
[m,n]=size(indices);
if m~=2,
    add_msg_board(sprintf('Distance defined only for two selected objects, but %i objects are selected.',m));
else
    indy1=indices(1,:);
    indy1=indy1(indy1>0);
    [msg,xyz1]=get_object(indy1,'xyz');
    xyz1=mean(xyz1,1);
    indy2=indices(2,:);
    indy2=indy2(indy2>0);
    [msg,xyz2]=get_object(indy2,'xyz');
    xyz2=mean(xyz2,1);
    if length(indy1)<5 || length(indy2)<5,
        add_msg_board('Distance between mean coordinates of');
    else
        add_msg_board('Distance between');
    end;
    add_msg_board(sprintf('%s and',description(indices(1,:))));
    add_msg_board(sprintf('%s',description(indices(2,:))));
    add_msg_board(sprintf('is %7.4f Å',norm(xyz1-xyz2)));
end;

