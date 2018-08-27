function [msg,snum]=symmetry_transform
% function msg=symmetry_transform
%
% checks whether selected objects define a symmetry or pseudosymmetry axis
% if so, transforms the structure to a frame where this axis is the z axis
% and where the center of the structure is on this axis
%
% the number of selected objects determines the multiplicity of the axis
% all selected objects must belong to the same structure, this structure is
% transformed
% all selected objects must have the same number of atoms (sorry)
%
% G. Jeschke, 2009

global model

msg.error=0;
msg.text='No error.';

snum=[];
if isfield(model,'selections')
    if model.selections>=2,
        for k=1:model.selections,
            cindices=model.selected{k};
            if isempty(snum),
                snum=cindices(1);
            elseif snum~=cindices(1),
                snum=-1;
                break;
            end;
        end;
    end;
end;

if isempty(snum),
    msg.error=1;
    msg.text='Selected objects do not define a symmetry axis.';
    return
end;

if snum<=0,
    msg.error=2;
    msg.text='For symmetry definition, all selected objects must belong to the same structure.';
    return
end;

[p0,v,msg]=symmetry_axis;
if ~isempty(p0) && ~isempty(v),
    v=v/norm(v);
    th=acos(v(3));
    if norm(v(1:2))>1e-6,
        v=v(1:2)/norm(v(1:2));
    else
        v=[0,0];
    end;
    phi=atan2(v(2),v(1));
    transmat1=affine('translation',-p0);
    transmat2=affine('Euler',[-phi,-th,0]);
    transform_structure(snum,{transmat1,transmat2});
else
    msg.error=3;
    msg.text='Selected objects have different numbers of atoms.';
end;

