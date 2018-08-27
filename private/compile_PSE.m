function compile_PSE
% uses Molpackage V2 of Gerhard Nieuwenhuiys.

global rootdir

defpath='definitions/';
defpath=strcat(rootdir,defpath);
outfile=strcat(defpath,'PSE.def');
psefile=strcat(defpath,'PSE.mat');
radii=load('radii.dat');

load colors

colfile=fopen('element_colors.txt','rt');
wfile=fopen(outfile,'wt');
fprintf(wfile,'%% Do not edit!\n%% This file exists only for information and as a template for additional definitions.\n');
fprintf(wfile,'%% Actual element definitions are in the binary file PSE.mat\n');
tags=':';
for k=1:92,
    colstr=fgetl(colfile);
    elm=pertable(k);
    tags=strcat(tags,elm.Symbol);
    tags=strcat(tags,':');
    fprintf(wfile,'begin %s\n',elm.Symbol);
    fprintf(wfile,'\t name %s\n',elm.Element);
    fprintf(wfile,'\t Mass %s\n',elm.RAmass);
    ilist=zeros(1,2);
    ilist=wr_isotopes(wfile,elm.NaturalIsotope1,ilist);
    ilist=wr_isotopes(wfile,elm.NaturalIsotope2,ilist);
    ilist=wr_isotopes(wfile,elm.NaturalIsotope3,ilist);
    [m,n]=size(ilist);
    ilist=ilist(2:m,:);
    fprintf(wfile,'\t covalent %4.2f\n',radii(k,1));
    fprintf(wfile,'\t vdW %4.2f\n',radii(k,2));
    fprintf(wfile,'\t color %s\n',colstr);
    fprintf(wfile,'end %s\n',elm.Symbol);
    element.symbol=elm.Symbol;
    element.name=elm.Element;
    element.mass=str2double(elm.RAmass);
    element.isotopes=ilist;
    element.covalent=radii(k,1);
    element.vdW=radii(k,2);
    element.color=colstr;
    colind=tag2id(upper(colstr),upper(color_tags));
    if isempty(colind),
        disp('Error at');
        disp(colstr);
        disp(k);
    end;
    element.rgb=colors(colind,:);
    pse(k)=element;
end;
fclose(wfile);

element_tags=tags;
save(psefile,'pse','element_tags');

function ilist=wr_isotopes(wfile,isotopes,ilist0)
[m,n]=size(ilist0);
ilist=ilist0;
if ~isempty(isotopes),
    isolist=textscan(isotopes,'%s');
    num=length(isolist{1});
    ilist=zeros(m+num,2);
    ilist(1:m,:)=ilist0;
    for k=1:num,
        id=char(isolist{1}(k));
        poi=findstr(id,'%');
        fprintf(wfile,'\t isotope %s %s\n',id(poi+1:end),id(1:poi-1));
        ilist(m+k,1)=str2double(id(poi+1:end));
        ilist(m+k,2)=str2double(id(1:poi-1));
    end;
end;