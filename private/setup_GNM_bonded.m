function Gamma=setup_GNM_bonded(Ca_coor)
% function Gamma=setup_GNM_bonded(Ca_coor)
%
% returns the Kirchhoff matrix GAMMA for a Gaussian network model (GNM)
% with inverse power 6 dependence of the force constant on interresidue
% distance and enhanced force constants for first and second neighbors 
% along the chain when given the Calpha coordinates
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
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param
global model

bond_force_1=80*3.8^(-6); % 1200
bond_force_2=80*6^(-6); % 80
exponent=6;

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Gamma=zeros(m);
for i=1:m-1,
    for j=i+1:m,
        r=norm(Ca_coor(i,:)-Ca_coor(j,:)); % distance of Calpha atoms
        cind0=model.coarse(model.current_structure).indices(i,:);
        cind1=model.coarse(model.current_structure).indices(j,:);
        if sum(abs(cind1(1:3)-cind0(1:3)))==0, % both C_alpha in the same chain
            diff=cind1(4)-cind0(4);
            switch diff
                case 1
                    cgamma=-bond_force_1*ENM_param.p_GNM_gamma;
                case 2
                    cgamma=-bond_force_2*ENM_param.p_GNM_gamma;
                otherwise
                    cgamma=-ENM_param.p_GNM_gamma*r^(-exponent);
            end;
        else
            cgamma=-ENM_param.p_GNM_gamma*r^(-exponent);
        end;
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
