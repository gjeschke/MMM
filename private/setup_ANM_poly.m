function Hessian=setup_ANM_poly(Ca_coor)
% function Hessian=setup_ANM_poly(Ca_coor)
%
% returns the Hessian for a fully connected anisotropic network model (ANM)
% when given the Calpha coordinates
% the force constants have inverse polynomial dependence on interresidue
% distance with order of the polynomial given in global variable
% ENM_param.p_ANM (defined in initialize_MMM.m)
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
% the force constant at a distance of 1 Å p_ANM_gamma is defined in
% initialize_MMM.m  and stored in global variable ENM_param
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param

Ca_dist=3.8; % typical C_alpha-C_alpha distance in Å
C_limit=1000;

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Hessian=zeros(3*m);
contacts=zeros(1,m);
for i=1:m-1,
    basi=3*(i-1);
    for j=i+1:m,
        basj=3*(j-1);
        rvec=Ca_coor(i,:)-Ca_coor(j,:); % distance vector of Calpha atoms
        r=norm(rvec);
        if r<C_limit,
            cgamma=-ENM_param.p_ANM_gamma*r^(-ENM_param.p_ANM);
            submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
            Hessian(basi+1:basi+3,basj+1:basj+3)=cgamma*submat;
            Hessian(basj+1:basj+3,basi+1:basi+3)=cgamma*submat;
        end;
    end;
end;

conn=0;
for i=1:m,
    basi=3*(i-1);
    submat=zeros(3);
    for j=1:m,
        basj=3*(j-1);
        submat=submat+Hessian(basi+1:basi+3,basj+1:basj+3);
    end;
    conn=conn+trace(submat);
    Hessian(basi+1:basi+3,basi+1:basi+3)=-submat;
end;

if conn==0, % ANM is unconnected
    Hessian=[];
end;
