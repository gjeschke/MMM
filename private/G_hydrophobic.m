function [G,msg]=G_hydrophobic(indices,z0,zm,eta,access,hydrophob)
% function [G,msg]=G_hydrophobic(indices,z0,zm,eta)
%
% hydrophobic energy contribution of residues an insertion into a bilayer
% the bilayer normal is assumed to be the z axis
% based on: A. Kessel, D. Shental-Bechor, Tu. Haliloglu, N. Ben-Tal,
%           Biophys. J. 2003, 85, 3431-344
% only the sidegroup contribution is considered, note that there is an
% error in Eq. (8) of the paper, see, e.g. Biophys. J. 2002, 82, 244-263
%
% G         free energy contribution in kcal/mol, empty variable is returned
%           if none of the addressed residues is an amino acid
% msg       structure that reports warnings or errors
%
% indices   array of indices that address amino acid residues, if any other
%           object is addressed, warning 2 is raised, MTSL (three-letter
%           code R1A) is treated as isoleucin-like
% z0        coordinate of the bilayer midplane, defaults to 0
% zm        width of a lipid monolayer (see Eq. (8) in Kessel, 2003, this
%           is wrongly labelled z0 in Fig. 2 of the paper), defaults to 15
%           Å
% eta       characteristic length parameter for membrane-water interface,
%           defaults to 0.5
% access    relative accessibility (optional, defaults to 1)
% hydrophob hydrophobicity values for the residues (optional, are extrcated, 
%           if not provided) 

global residue_defs

% kcal2kJ=4.1868; % conversion kcal/mol -> kJ/mol, unused

G=[];
msg.error=1;
msg.text='No amino acid addressed by indices.';

if nargin<1 || isempty(indices),
    return;
end;

if nargin<2,
    z0=0;
end;

if nargin<3,
    zm=15;
end;

if nargin<4,
    eta=0.5;
end;

non_amino=0;

[m,n]=size(indices);

if nargin<5,
    access=ones(1,m);
end;

succ=0;
G=0;
for k=1:m,
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    if length(cindices)~=4,
        non_amino=non_amino+1;
    else
        if nargin<6, % extract hydrophobicity
            [message,tag]=get_residue(cindices,'name');
            if strcmpi(tag,'R1A'),
                tag='ILE';
            end;
            resnum=tag2id(tag,upper(residue_defs.restags));
            if ~isempty(resnum),
                info=residue_defs.residues(resnum);
                [message,xyz]=get_residue(cindices,'xyz');
                coor=mean(xyz,1);
                pr=1/(1+exp(eta*(abs(coor(3)-z0)-zm)));
                G=G+access(k)*info.hydrophobicity*pr;
                succ=1;
            end;
        else
            [message,xyz]=get_residue(cindices,'xyz');
            coor=mean(xyz,1);
            pr=1/(1+exp(eta*(abs(coor(3)-z0)-zm)));
            G=G+access(k)*hydrophob(k)*pr;        
            succ=1;
        end;
    end;
end;

if ~succ, 
    G=[]; 
else
    G=10*G;
    msg.error=0;
    msg.text='No error.';
end;