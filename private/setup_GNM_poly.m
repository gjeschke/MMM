function Gamma=setup_GNM_poly(Ca_coor)
% function Gamma=setup_GNM_poly(Ca_coor)
%
% returns the Kirchhoff matrix GAMMA for a Gaussian network model (GNM)
% with inverse polynomial dependence of the force constant on interresidue
% distance when given the Calpha coordinates
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% and:
% L. Yang, G. Song, R L. Jernigan, Proc. Natl. Acad. Sci. 106 (2009)
% 12347–12352.
%
% the basic parameters of the GNM (exponent ENM_param.p_GNM and generic
% force constant ENM_param.p_GNM_gamma) are defined in initialize_MMM.m and
% stored in global variable ENM_param
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Gamma=zeros(m);
for i=1:m-1,
    for j=i+1:m,
        r=norm(Ca_coor(i,:)-Ca_coor(j,:)); % distance of Calpha atoms
        cgamma=-ENM_param.p_GNM_gamma*r^(-ENM_param.p_GNM);
        Gamma(i,j)=cgamma;
        Gamma(j,i)=cgamma;
    end;
end;

conn=0;
for i=1:m,
    restore=sum(Gamma(i,:));
    conn=conn+restore;
    Gamma(i,i)=-restore;
end;

if conn==0, % GNM is unconnected
    Gamma=[];
end;
