function Hessian=setup_ANM_Hinsen(Ca_coor)
% function Hessian=setup_ANM_Hinsen(Ca_coor)
%
% returns the Hessian for a fully connected anisotropic network model (ANM)
% when given the Calpha coordinates
% the force constants have a distance dependence as suggested by Hinsen et
% al.
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% and:
% K. Hinsen, A. Petrescu, S. Dellerue, M. Bellissent-Funel, G. Kneller,
% Chem. Phys. 261 (2000) 25.
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010


r_local=4;
RT=2.4776e3;

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Hessian=zeros(3*m);
contacts=zeros(1,m);
for i=1:m-1,
    basi=3*(i-1);
    for j=i+1:m,
        basj=3*(j-1);
        rvec=Ca_coor(i,:)-Ca_coor(j,:); % distance vector of Calpha atoms
        r=norm(rvec);
        if r<r_local,
            cgamma=(-8.6e4*r-2.39e5)/RT;
        else
            cgamma=-(1.28e8/r^6)/RT;
        end;   
        submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
        Hessian(basi+1:basi+3,basj+1:basj+3)=cgamma*submat;
        Hessian(basj+1:basj+3,basi+1:basi+3)=cgamma*submat;
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
