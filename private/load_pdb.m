function [msg,snum]=load_pdb(fname)
% function [msg,snum]=load_pdb(fname)
%
% load a PDB file from the local disk
%
% fname     PDB file name with extension
%
% msg       a message record reporting errors
% snum      the number of the imported structure, is empty if the download
%           failed

global hMain
global graph_settings

snum=[];

msg.error=0;
msg.text='No error.';

fname=which(fname);

if isempty(fname),
    msg.error=1;
    msg.text='ERROR: PDB file not found.';
    return
end;

[msg,snum]=add_pdb(fname);

if hMain.virgin,
    hMain.virgin=0;
    % initialize display
    axes(hMain.axes_model);
    cla reset;
    axis equal
    axis off
    set(gca,'Clipping','off');
    set(gcf,'Renderer','opengl');
    hold on
    hMain.camlight=camlight;

    axes(hMain.axes_model);
    axis equal
    axis off
    set(gca,'Clipping','off');
    set(gcf,'Renderer','opengl');
    hold on

    view(graph_settings.az,graph_settings.el);
    lighting gouraud 
    material shiny

    axes(hMain.axes_frame);
    cla reset;
    axis equal
    axis off
    hold on

    plot3([0 1],[0 0],[0 0],'r','LineWidth',2); % x axis
    plot3([0 0],[0 1],[0 0],'g','LineWidth',2); % y axis
    plot3([0 0],[0 0],[0 1],'b','LineWidth',2); % z axis
    axis([-0.1,1.1,-0.1,1.1,-0.1,1.1]);

    guidata(hMain.axes_frame,hMain);
    guidata(hMain.axes_model,hMain);
end;

set(gcf,'Pointer','arrow');
