function transform_model(transform)
% function transform_model(transform)
%
% Perform a coordinate transform on all structures in the model
% transform is an affine transformation matrix (4x4) or a cell array of such
% matrices
%
% transforms atom coordinates, label coordinates, and solvent accessible
% surfaces
%
% after the transform, the graphics may no longer be conssitent with
% coordinates
% the client is responsible for redisplay

global model

if isempty(model.structures),
    return;
end;

for snum=1:length(model.structures),
    transform_structure(snum,transform);
end;

