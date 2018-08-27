function coor=make_peptide_backbone(fi,psi,omega),
%
% Creates peptide backbone coordinates from given dihedral angles assuming
% fixed bond lengths and bond angles
%
% fi, psi, omega     are vectors of dihedral angles of the same length,
%                    which is the number n of residues, psi(n) and omega(n) 
%                    are irrelevant for the backbone coordinates, but must
%                    be numbers
% coor   are Cartesian coordinates, starting with N of the first residue at
%        (0,0,0) and ending with C of the last residue
%
% based on: H. Sugeta, T. Miyazawa, Biopolymers, 1967, 5, 673-679.
%
% (c) G. Jeschke, 2007

% Default bond lengths
b1=1.416; % N-C_alpha
b1vec=[b1,0,0];
b2=1.529; % C_alpha-C
b2vec=[b2,0,0];
b3=1.320; % C-N
b3vec=[b3,0,0];

% Default bond angles
fi1=110*pi/180; % N-C_alpha-C
fi2=114.2*pi/180; % C_alpha-C-N
fi3=121*pi/180; % C-N-C_alpha

acoor=[0;0;0];
coor=zeros(3*length(fi),3);
A=eye(3);

for k=1:length(fi),
    coor(3*k-2,:)=acoor'; % coordinates of N
    acoor=acoor+A*b1vec'; % coordinates of C_alpha
    coor(3*k-1,:)=acoor';
    cfi=cos(fi1);
    sfi=sin(fi1);
    ctau=cos(fi(k));
    stau=sin(fi(k));
    A12=[-cfi,-sfi,0;sfi*ctau,-cfi*ctau,-stau;sfi*stau,-cfi*stau,ctau];
    A=A*A12;
    acoor=acoor+A*b2vec'; % coordinates of C
    coor(3*k,:)=acoor';
    cfi=cos(fi2);
    sfi=sin(fi2);
    ctau=cos(psi(k));
    stau=sin(psi(k));
    A23=[-cfi,-sfi,0;sfi*ctau,-cfi*ctau,-stau;sfi*stau,-cfi*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of next N
    cfi=cos(fi3);
    sfi=sin(fi3);
    ctau=cos(omega(k));
    stau=sin(omega(k));
    A31=[-cfi,-sfi,0;sfi*ctau,-cfi*ctau,-stau;sfi*stau,-cfi*stau,ctau];
    A=A*A31;    
end;