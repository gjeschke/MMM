function delta=delta_inv_T1(d,d0,lambda,mindel,maxdel)
%
%
% Relaxation enhancement function Delta(1/T1) according to:
% D. Marsh. 2001. Polarity and permeation profiles in lipid membranes. 
% Proc. Natl. Acad. Sci. USA, 98:7777-7782.
% 
% see also: B G. Dzikovski, V. A. Livshits, D. Marsh. 2003. 
%           Biophys. J. 85: 1005-102. 
% Input:
% d         distance from the membrane central plane
% d0        distance of the point of maximum gradient from the membrane
%           central plane
% lambda    decay constant; d, d0, lambda must have the same units
% mindel    minimum relaxation enhancement (outside membrane, water)
% maxdel    maximum relaxation enhancement (membrane center)
%
% Output:
% delta     Delta(1/T1)

d=15-abs(d); % redefinition, so that d is the distance from the headgroup

delta=maxdel+(mindel-maxdel)/(1+exp((d-d0)/lambda));