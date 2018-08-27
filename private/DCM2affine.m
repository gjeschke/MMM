function transmat=DCM2affine(DCM)
% function transmat=DCM2affine(DCM)
%
% Converts a direction cosine matrix, obtained, e.g., by SpinCalc or
% Euler2DCM, to an affine transformation matrix for use with affine_trafo,
% affine_trafo_point, and affine_trafo_vector
%
% G. Jeschke, 2009

transmat=eye(4); % initialize identity transformation

transmat(1:3,1:3)=DCM; % substitute "rotation block"
