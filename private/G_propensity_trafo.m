function [G,msg]=G_propensity_trafo(indices,zm,thickness,access,propens,theta,phi)
% function [G,msg]=G_propensity_trafo(indices,zm,thickness,access,propens,theta,phi)
%
% hydrophobic energy contribution of residues an insertion into a bilayer
% the bilayer normal is assumed to be the z axis
% based on: L. Adamian, V. Nanda, W. F. De Grado, J. Liang,
%           Proteins, 2005, 59:496-509
%
% G         free energy contribution in kcal/mol, empty variable is returned
%           if none of the addressed residues is an amino acid
% msg       structure that reports warnings or errors
%
% indices   array of indices that address amino acid residues, if any other
%           object is addressed, warning 2 is raised, MTSL (three-letter
%           code R1A) is treated as isoleucin-like
% zm        coordinate of the bilayer midplane
% thickness width of a lipid bilayer
% access    relative accessibility
% propens   propensities TMLIP-H and TMLIP-P from Adamian et al.
% theta     Euler angle for rotation about y'
% phi       Euler angle for rotation about z''

% kcal2kJ=4.1868; % conversion kcal/mol -> kJ/mol, unused

transmat=affine('Euler',[phi,theta,0]);

RT=8.314*(273.15+27)/4.1868; % RT at 27 °C in kcal/mol

G=[];
msg.error=1;
msg.text='No amino acid addressed by indices.';

if nargin<1 || isempty(indices),
    return;
end;

non_amino=0;

[m,n]=size(indices);


succ=0;
G=0;
for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    if length(cindices)~=4,
        non_amino=non_amino+1;
    else
%        [message,xyz]=get_residue(cindices,'xyz');
        [message,xyz]=get_residue(cindices,'Calpha');
        if ~isempty(xyz),
            coor=mean(xyz,1);
            coor=affine_trafo_point(coor,transmat);
            coor(3)=coor(3)-zm;
            position=abs(coor(3))/thickness;
            propensity=0;
            if position<=0.5, % headgroup layer
                propensity=propens(k,1); 
            end;
            if position<0.25, % core
                propensity=propens(k,2);
            end;
            G=G+access(k)*propensity;        
            succ=1;
        end;
    end;
end;

if ~succ, 
    G=[]; 
else
    G=RT*G;
    msg.error=0;
    msg.text='No error.';
end;