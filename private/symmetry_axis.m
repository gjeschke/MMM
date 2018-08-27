function [p0,symaxis,msg]=symmetry_axis
% function axis=symmetry_axis
%
% Determines an n-fold symmetry axis, assuming that the n selected objects 
% are related by rotational symmetry
%
% empty matrices are returned if the number of atoms/atom locations in the
% selected objects is inconsistent or if less than two objects are selected
%
% symaxis   unit vector along the symmetry axis direction
% p0        point on the symmetry axis, near the center of the structure
%
% G. Jeschke, 2009

global model

min_ext=2;  % minimum extension of the point cloud in axis direction for 
            % considering the axis as well defined

msg.error=1;
msg.text='Selected objects do not define symmetry axis';

symaxis=[];
p0=[];

xyz=[];

poi=0;
if isfield(model,'selections')
    if model.selections>=2,
        for k=1:model.selections,
            cindices=model.selected{k};
            cindices=cindices(cindices>0);
            [msg,xyza]=get_object(cindices,'xyz','-nowater');
            if isempty(xyz),
                xyz=xyza;
                poi=1;
            else
                try
                    xyz=xyz+xyza;
                    poi=poi+1;
                catch my_exception
                    return
                end;
            end;
        end;
        if poi>1 && ~isempty(xyz),
            xyz=xyz/poi;
            [m,n]=size(xyz);
            if m>=4
                [p0,v]=rmsd_line_3D(xyz);
                symaxis=v/norm(v);
                if norm(v)>min_ext,
                    msg.error=0;
                    msg.text='No error.';
                else
                    msg.error=999;
                    msg.text='Warning. Axis may be ill defined.';
                end;
            end;
        end;
    end;
end;

