function ENM_param=set_ANM(ENM_param,ANM_choice)
% Sets parameters for the anisotropic elastic network model according to
% a user-selected parametrization
% energy unit for force constants is kJ, distance unit is Å (sorry)
%
% ENM_param     structure with ealstic network model parameters, the
%               parameters defined below in ANM_choice can also be supplied
%               as fields of ENM_param
% ANM_choice    user selection of model, structure with fields
%               .parametrization    type of ANM, can be
%                                   'cutoff10'  uniform force constants
%                                               with uniform cutoff at 10 Å
%                                   'cutoff13'  uniform force constants
%                                               with uniform cutoff at 13 Å
%                                   'Jeschke'   r^-6 distance scaling with
%                                               direct and next neighbor
%                                               bonded force constants
%                                               enhanced by factor 10000
%                                   'Hinsen'    r^-6 distance scaling with
%                                               direct neighbor bonded 
%                                               force constant enhanced
%                                               similar to K. Hinsen, A. 
%                                               Petrescu, S. Dellerue, M. 
%                                               Bellissent-Funel, G. Kneller,
%                                               Chem. Phys. 261 (2000) 25.
%                                   'ed-ENM'    the essential-dynamics
%                                               refined ealstic network
%                                               model of Orellana et al.,
%                                               JCTC, 2010, 6, 2910-2923
%                                   'ed-ENM-p'  like ed-ENM, but without
%                                               cutoff
%               .imANM              optional flag, if true, the implicit
%                                   membrane ANM of Lezon & Bahar, Biophys.
%                                   J. 2012, 102, 1331-1340 (imANM) is used
%               .mass_weighting     residue-mass weighting of force
%                                   constants and coordinates, optional
%                                   flag, defaults to false, UNTESTED for
%                                   .mass_weighting = true
%
% G. Jeschke, 2012

kcal2kJ=4.1868; % kcal to kJ
nm2A=10; % nanometer to Å
RT=8.314e-3*298; % 8.314 J mol^(-1) K^(-1)· 298 K
rCaCa=3.8; % mean Calpha-Calpha distance for direct neighbor residues

if nargin<2,
    ANM_choice=ENM_param;
end;

if isfield(ANM_choice,'imANM'),
    ENM_param.imANM=ANM_choice.imANM;
else
    ENM_param.imANM=false;
end;

if isfield(ANM_choice,'mass_weighting'),
    ENM_param.mass_weighting=ANM_choice.mass_weighting;
else
    ENM_param.mass_weighting=false;
end;

if isfield(ANM_choice,'parametrization'),
    ENM_param.parametrization=ANM_choice.parametrization;
    switch ANM_choice.parametrization
        case 'cutoff10'
            ENM_param.p_ANM=0; %  exponent for distance-dependent force constants ANM
            ENM_param.cutoff_const_ANM=10;
            ENM_param.cutoff_log_ANM=0;
            ENM_param.connect=0;
            ENM_param.C_connect=0;
            ENM_param.p_connect=0;
            ENM_param.C_cart=6;
        case 'cutoff13'
            ENM_param.p_ANM=0; %  exponent for distance-dependent force constants ANM
            ENM_param.cutoff_const_ANM=13;
            ENM_param.cutoff_log_ANM=0;
            ENM_param.connect=0;
            ENM_param.C_connect=0;
            ENM_param.p_connect=0;
            ENM_param.C_cart=6;
        case 'Hinsen'
            ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
            ENM_param.cutoff_const_ANM=1e6;
            ENM_param.cutoff_log_ANM=0;
            ENM_param.connect=1;
            ENM_param.C_connect=(8.6e5*nm2A^(-3)*rCaCa-2.39e5*nm2A^(-2))/RT;
            ENM_param.p_connect=100;
            ENM_param.C_cart=(128*nm2A^(4))/RT;
        case 'Jeschke'
            ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
            ENM_param.cutoff_const_ANM=1e6;
            ENM_param.cutoff_log_ANM=0;
            ENM_param.connect=2; 
            ENM_param.C_connect=19.62*kcal2kJ;
            ENM_param.p_connect=1.789;
            ENM_param.C_cart=6*kcal2kJ;
        case 'ed-ENM'
            ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
            ENM_param.cutoff_const_ANM=-2.9;
            ENM_param.cutoff_log_ANM=2.9;
            ENM_param.connect=3;
            ENM_param.C_connect=(2.5e4*nm2A^(-2))/RT;
            ENM_param.p_connect=2;
            ENM_param.C_cart=(19.534*nm2A^(4))/RT;
        case 'ed-ENM-p'
            ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
            ENM_param.cutoff_const_ANM=1e6;
            ENM_param.cutoff_log_ANM=0;
            ENM_param.connect=3;
            ENM_param.C_connect=60*kcal2kJ;
            ENM_param.p_connect=2;
            ENM_param.C_cart=6*kcal2kJ;
    end;
end;
