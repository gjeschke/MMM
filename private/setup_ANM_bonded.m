function Hessian=setup_ANM_bonded(Ca_coor)
% function Hessian=setup_ANM_bonded(Ca_coor)
%
% returns the Hessian for a fully connected anisotropic network model (ANM)
% when given the Calpha coordinates
% the force constants have inverse 6th power dependence on interresidue
% distance, except for next neighbors (i,i+1) and second neighbors (i,i+2)
% along the polypetide chain
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% and:
% L. Yang, G. Song, R L. Jernigan, Proc. Natl. Acad. Sci. 106 (2009)
% 12347–12352.
% and:
% L. Orellana, M. Rueda, C. Ferrer-Costa, J. R. Lopez-Blanco, P. Chacon, M.
% Orozco, J. Chem. Theory Comput. 2010, 6, 2910-2923
%
% parameters are defined in
% initialize_MMM.m  and stored in global variable ENM_param:
%      .p_ANM               exponent for distance-dependent force constant
%      .cutoff_const_ANM    constant contribution to cutoff distance
%      .cutoff_log_ANM      logarithmic scaling of cutoff with number of
%                           residues
%      .connect             number of connected neighbors with larger force
%                           constants
%      .C_connect           connected force constant for direct neighbor in
%                           kcal/Å^2
%      .p_connect           exponent for inverse scaling of connected force
%                           constant with distance in residue sequence
%      .C_cart              unconnected force constant for direct neighbor 
%                           in kcal/Å^2
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param
global model

Hessian=setup_ANM_bonded_parallel(Ca_coor,model,ENM_param);

