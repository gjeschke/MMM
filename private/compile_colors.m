function compile_colors
% Compiles a color definition file based on the SVG color definitions of
% the W3C consortium, currently only RGB definitions are given
% use Matlab function rgb2hsv to convert to HSV colors or export option
% CMYK for print/export to convert to CMYK colors
%
% also provides an .m file with an appropriate helper function
% localDefineColors for the Matlab colorPalette GUI
%

global rootdir

defpath=[rootdir 'definitions/'];
outfile=strcat(defpath,'colors.def');
coloutfile=strcat(defpath,'colors.mat');

infile='svgcolors.txt';

colfile=fopen(infile,'rt');
wfile=fopen(outfile,'wt');
fprintf(wfile,'%% Do not edit!\n%% This file exists only for information and as a template for additional definitions.\n');
fprintf(wfile,'%% Actual color definitions are in the binary file colors.mat\n');
tags=':';
defs=textscan(colfile,'%s%s%s');
fclose(colfile);
k=0;
colors=zeros(200,3);
while 1
    k=k+1;
    if strfind(char(defs{1}(k)),'%'), break; end;
    tags=strcat(tags,char(defs{1}(k)));
    tags=strcat(tags,':');
    rgbtriple=str2num(char(defs{3}(k)));
    colors(k,:)=rgbtriple;
    names{k}=char(defs{1}(k));
    fprintf(wfile,'begin %s\n\t rgb %3i%4i%4i\n',char(defs{1}(k)),rgbtriple);
    fprintf(wfile,'end %s\n',char(defs{1}(k)));
end;
fclose(wfile);

n=k-1;
colors=colors(1:n,:)/255;

color_tags=tags;

save(coloutfile,'colors','color_tags','names');

disp(sprintf('%i standard colors defined',n));

