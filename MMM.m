% MMM Multiscale Modeling of Macromolecules
%
%   MMM
% 
% An open-source program for visualization, inspection, and improvement of
% models of proteins and protein assemblies based on restraints from
% multiple experimental techniques
%
% Inputs:
%    n/a
%
% Outputs:
%    n/a
%                   
% Example: 
%    see http://www.epr.ethz.ch/software/MMM_2011.2_manual.pdf
%
% Other m-files required:   none
%
% Subfunctions:             MMM_prototype.m
%                           MMM (directory)
%
% MAT-files required:       none
%
% References:
%   1)  Y. Polyhach, E. Bordignon and G. Jeschke: Rotamer libraries of
%   spin-labelled cysteins for protein studies. PCCP 13(6) (2011) 2356-2366
%
%   2)  Y. Polyhach and G. Jeschke: Prediction of favourable sites for spin
%   labelling of proteins. Spectroscopy (AMSTERDAM, Netherlands) 2010,
%   24(6), 651-659.

% Author listings
% ===============
% Author:       Gunnar Jeschke
% Work address: Building HCI
%               Wolfgang-Pauli-Str. 10
%               8093 Zurich, Switzerland
% Email:        gjeschke@ethz.ch
% Website:      http://www.epr.ethz.ch/
% 2010;         Last revision: 2010
%
% Author:       Morgan Bye
% Work address: Henry Wellcome Unit for Biological EPR
%               University of East Anglia
%               NORWICH, UK
% Email:        morgan.bye@uea.ac.uk
% Website:      http://www.morganbye.net/eprtoolbox/
% Aug 2012;     Last revision: 29-August-2012

% Version history
% ===============
% Aug 12   MB      Update for dual monitor displays, added header section
%
% 2011     GJ      Update for St Stoll credit
%
% 2010     GJ      Initial release

h=figure;
set(h,'NumberTitle','off');
set(h,'Name','Loading MMM...');
set(h,'Color',[20,43,140]/255);
set(h,'MenuBar','none');
set(h,'ToolBar','none');

% This only gets the total desktop size
% screen_size=get(0,'ScreenSize');
% movegui(h,'center');

% For 2 monitor spanning desktops it looks odd, as the image is stretched
% across 50% of both monitors.

screen_size = get(0,'MonitorPositions');

switch size(screen_size,1)
    case 1
        set(h,'Position',[...
            round(screen_size(3)/4),...
            round(screen_size(4)/4),...
            round(screen_size(3)/2),...
            round(screen_size(4)/2)]);
    otherwise
        % Use monitor 1 (linux)/primary monitor (windows)
        set(h,'Position',[...
            round(screen_size(1,3)/4),...
            round(screen_size(1,4)/4),...
            round(screen_size(1,3)/2),...
            round(screen_size(1,4)/2)]);
end

[img,map]=imread('MMM_title.jpg');
colormap(map);
image(img);
axis off
axis equal
MMM_prototype;
clear img
close(h);

