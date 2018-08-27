function xyz = get_non_interacting_label_position(backbone,rot_lib)
% xyz = get_non_interacting_label_position(backbone,libname)
%
% Provides the mean coordinates xyz of a label that is not interacting with
% backbone or sidechain atoms of the biomacromolecule
% use this function only if structural context is not get fully known
% the coordinates are the mean coordinates of an unrestrained label
% according to the given rotamer library, translated to the local residue
% frame
%
% backbone  (3x3) coordinate array for the attachment frame
%           first line: atom on x axis (backbone N for peptide chains)
%           second line: origin, (Calpha atom for peptide chain)
%           third line: atom in xy plabe (C(=0) atom for peptide chain)
% rot_lib   rotamer library
%
% xyz       coordinates of the mean label position in the supplied 
%           backbone frame
%
% G. Jeschke, 12.10.2015

x= backbone(1,:) - backbone(2,:); % x axis is along C_alpha-N bond
x=x/norm(x);    % unit vector along x
yp=backbone(3,:) - backbone(2,:); % y axis is in the plane spanned by x axis and C-Ca bond
yp=yp/norm(yp);
z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z=z/norm(z);
y=cross_rowvec(z,x); % real (corrected) y axis 
dircos=[x;y;z];
Rp=dircos; % rotation matrix for conversion to standard frame
NO = get_relative_label(rot_lib);
xyz = NO*Rp + backbone(2,:);



function NO = get_relative_label(rot_lib)

midNO = rot_lib.usefull_atoms.midNO;
pops = rot_lib.calibration.pop;
NO = zeros(1,3);
for k = 1:length(rot_lib.library),
    coor = rot_lib.library(k).ecoor;
    NO = NO + pops(k)*(coor(midNO(1),2:4) + coor(midNO(2),2:4))/2;
end;
NO = NO/(sum(pops));