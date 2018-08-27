function [Hessian,contacts]=setup_ANM(Ca_coor,restypes,rindices)
% function [Hessian,contacts]=setup_ANM(Ca_coor,restypes,rindices)
%
% returns the Hessian for an anisotropic network model (ANM)
% when given the Calpha coordinates
% a vector contacts of contact numbers for each residue is also returned 
%
% restypes  optional argument that provides the amino acid types
% rindices  residue indices, if this argument is given, intra- and
%           interchain force constants are treated separately
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% note that there is a sign error in the expression for Gamma_ii in Eq.
% (19) of that paper
%
% basic parameters of the GNM (cutoff distance ENM_param.rc and generic
% force constant ENM_param.gamma) are defined in initialize_MMM.m and
% stored in global variable ENM_param
%
% if argument restypes is given, amino acid specific force constants are
% used according to:
% K. Hamacher, J. A. McCammon, J. Chem. Theory Comput. 2006, 2, 873-878.
% in our hands, this leads to worse results than the simpler model,
% NOT RECOMMENDED
%
% intrachain force constants are derived from:
% S. Miyazawa, R. L. Jernigan, J. Mol. Biol. 1996, 256, 623-644
% interchain force constant are derived from:
% O. Keskin, I. Bahar, A.Y. Badretdinov, O.B. Ptitsyn, R.L. Jernigan,
% Protein Sci.  1998, 7: 2578-2586
% both potentials are defined in file contact_potentials
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param

if nargin>1,
    load contact_potentials;
    potentials=contact_energies;
    [mp,np]=size(potentials);
    nzelm=mp*(mp-1)/2+mp;
end;

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Hessian=zeros(3*m);
contacts=zeros(1,m);
for i=1:m-1,
    basi=3*(i-1);
    for j=i+1:m,
        basj=3*(j-1);
        rvec=Ca_coor(i,:)-Ca_coor(j,:); % distance vector of Calpha atoms
        r=norm(rvec);
        if r>ENM_param.rc_ANM, % set force constant for pair of residues
            cgamma=0;
            submat=zeros(3);
        elseif nargin<2
            contacts(i)=contacts(i)+1;
            contacts(j)=contacts(j)+1;
            cgamma=-ENM_param.gamma;
            submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
        else
            contacts(i)=contacts(i)+1;
            contacts(j)=contacts(j)+1;
            submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
            type1=restypes(i);
            type2=restypes(j);
            if type2<type1, swap=type2; type2=type1; type1=swap; end;
            if nargin>2,
                if rindices(i,2)==rindices(j,2),
                    potentials=contact_energies;
                else
                    potentials=solvent_mediated;
                end;
            end;
            if type1==0,
                if type2==0,
                    cgamma=sum(sum(potentials))/nzelm;
                else
                    cgamma=(sum(potentials(type2,:))+sum(potentials(:,type2))-potentials(type2,type2))/mp;
                end;
            else
                cgamma=potentials(type1,type2);
            end;
            if abs(cgamma)<eps,
                disp('Unexpectedly low force constant');
            end;
        end;
        Hessian(basi+1:basi+3,basj+1:basj+3)=cgamma*submat;
        Hessian(basj+1:basj+3,basi+1:basi+3)=cgamma*submat;
    end;
end;

for i=1:m,
    basi=3*(i-1);
    submat=zeros(3);
    for j=1:m,
        basj=3*(j-1);
        submat=submat+Hessian(basi+1:basi+3,basj+1:basj+3);
    end;
    Hessian(basi+1:basi+3,basi+1:basi+3)=-submat;
end;

if sum(contacts)==0, % ANM is unconnected
    Hessian=[];
    contacts=[];
end;
