function ecoor = anchor_correction(ecoor,corr)
% function ecoor = anchor_correction(ecoor,corr)
% 
% Correction of an extended coordinate set of an RNA loop for matching the
% final anchor, use mk_RNA_loop_backbone.m and mk_RNA.m to obtain
% corr and ecoor
% the correcting translation and rotation is proportional to the distance
% between the current atom and the initial P atom, divided by the distance
% between the target P atom and the initial P atom
% this distributes the deviation between target and final P atom of the 
% original construct, so that each individual bond length, bond angle, and 
% torsion angle is only weakly affected
%
% ecoor [m,4] coordinate array for m atoms, first column is residue
%       number, columns 2-4 are Cartesian coordinates
% corr  correction parameters, structure   
%       .trans      translation vector
%       .rot        4-element vector, elements 1:3 specify rotation axis,
%                   element 4 specifies rotation angle in degrees
%       .length     length from first atom to anchor point
%       .Pinitial   coordinates of initial P atom       
%
% G. Jeschke, 26.12.2017

[m,~] = size(ecoor);
for k = 1:m
    progress = norm(ecoor(k,2:4)-corr.Pinitial)/corr.length;
    tm1 = affine('translation',progress*corr.trans);
    tm2 = affine('rotn',progress*pi*corr.rot(4)/180,corr.rot(1:3));
    tm = tm1*tm2;
    ccoor = [ecoor(k,2:4) 1]*tm';
    ecoor(k,2:4) = ccoor(1:3);
end
