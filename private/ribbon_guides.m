function [backbone,rung,normal]=ribbon_guides(Ca_coor,C_coor,N_coor,steps_per_residue,sheetflag)
% function [backbone,rung,normal]=ribbon_guides(Ca_coor,C_coor,N_coor,steps_per_residue)
%
% Create coordinates of ribbon guide curves for a protein backbone, using a
% slightly modified version of the Carson-Bugg algorithm (J. Mol. Graph.
% 1986, 4, 121-122)
%
% the guide curves are space by one unit of the input coordinates, usually
% 1 Angstroem
%
% uses: function cardinal_splines
%
% Ca                    coordinates of Calpha atoms [m,3]
% N_coor                coordinates of peptide nitrogen atoms [m,3]
% C_coor                coordinates of carbonyl carbon atoms [m,3]
% steps_per_residue     number of guide points in a segment between two
%                       residues
%
% backbone              coordinates of backbone guide points
% rung                  unit vector from backbone guide points to edge guide points
% normal                unit vector perpendicular to peptide plane
% 
% (c) G. Jeschke, 2009

Tension=0; % assume cardinal splines with zero tension

if nargin<5,
    sheetflag=0;
end;

[m,n]=size(Ca_coor);
edge0=zeros(m,n);
rung0=edge0;
dirsign=-1;

for k=1:m-1,
    a=N_coor(k,:)-Ca_coor(k,:); % x axis in peptide plane
    b=C_coor(k,:)-Ca_coor(k,:);
    c=cross_rowvec(a,b);
    c=c/norm(c); % normal to peptide plane
    d=cross_rowvec(c,a);
    d=d/norm(d); % y axis in peptide plane
    e=Ca_coor(k+1,:)-Ca_coor(k,:);
    f=cross_rowvec(e,c); % vector in peptide plane perpendicular to C_alpha-C_alpha(k+1)
    f=f/norm(f);
    if sheetflag, dirsign=dirsign*-1; end; % flip side, if necessary
    edge0(k,:)=Ca_coor(k,:)+dirsign*f;
    rung0(k,:)=dirsign*f;
end;
edge0(m,:)=Ca_coor(m,:)+c;
rung0(m,:)=dirsign*f;

guide_length=(m-1)*steps_per_residue+1;
backbone=zeros(guide_length,3);
edge=backbone;
normal=backbone;
rung=backbone;

backbone0=[Ca_coor(1,:); Ca_coor; Ca_coor(m,:)];
rung0=[rung0(1,:); rung0; rung0(m,:)];    

for k=1:m-1,
    basnum=(k-1)*steps_per_residue;
    backbone(basnum+1:basnum+steps_per_residue+1,:)=cardinal_spline(backbone0(k:k+3,:),Tension,steps_per_residue);
    rung(basnum+1:basnum+steps_per_residue+1,:)=cardinal_spline(rung0(k:k+3,:),Tension,steps_per_residue);            
end;
for k=1:guide_length-1,
    e=backbone(k+1,:)-backbone(k,:);
    f=rung(k,:);
    g=cross_rowvec(e,f);
    normal(k,:)=g/norm(g);
    rung(k,:)=f/norm(f);
end;
normal(guide_length,:)=normal(guide_length-1,:);
rung(guide_length,:)=rung(guide_length,:)/norm(rung(guide_length,:));

