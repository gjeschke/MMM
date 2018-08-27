function transform_structure_ensemble(snum,transform,bilayer)
% function transform_structure(snum,transform,bilayer)
%
% Perform a coordinate transform on the structure with number snum with
% individual rotations/translations for the individual models in an
% ensemble
% transform is a cell array of affine transformation matrices (4x4)
% with the members corresponding to the individual models
%
% transforms only atom coordinates, spin labels, and precomputed, but not
% attached spin label rotamers
%
% after the transform, the graphics may no longer be consistent with
% coordinates
% the client is responsible for redisplay
%
% bilayer   optional flag, if present and true, bilayer is not invalidated

global model

if nargin<3,
    bilayer=false;
end;

if ~bilayer && isfield(model.info{snum},'bilayer') && ~isempty(model.info{snum}.bilayer),
    if ~isempty(model.info{snum}.bilayer.graphics),
        for k=1:3,
            if ishandle(model.info{snum}.bilayer.graphics(k)),
                delete(model.info{snum}.bilayer.graphics(k));
            end;
        end;
    end;
    model.info{snum}=rmfield(model.info{snum},'bilayer');
    add_msg_board('Warning: Bilayer became invalid on coordinate transformation and was removed');
end;

if isfield(model.info{snum},'cryst'),
    model.info{snum}=rmfield(model.info{snum},'cryst');
    add_msg_board('Warning: Unit cell definition became invalid on coordinate transformation and was removed');
end;


matrix=transform{1};
p0=matrix(1:3,4)';

% set_object(indices,'transform',matrix);
if length(model.structures)>=snum,
    for k=1:length(model.structures{snum}),
        for kk=1:length(model.structures{snum}(k).xyz),    
            matrix = transform{kk};
            xyz=model.structures{snum}(k).xyz{kk};
            [mm,nn]=size(xyz);
            xyz=[xyz ones(mm,1)];
            xyz=matrix*xyz';
            xyz=xyz';
            model.structures{snum}(k).xyz{kk}=xyz(:,1:3);
        end;
    end;
end;

model.info{snum}.center=model.info{snum}.center+p0;
if isfield(model,'labels'),
    for k=1:length(model.labels),
        tag=strtok(model.labels(k).adr,']');
        tags=model.structure_tags;
        ids=model.structure_ids;
        id=tag2id(tag(2:end),tags,ids);
        if id==snum,
            lindices = resolve_address(model.labels(k).adr);
            matrix = transform{lindices(3)};
            xyz=model.labels(k).NOpos(:,1:3);
            [mm,nn]=size(xyz);
            xyz=[xyz ones(mm,1)];
            xyz=matrix*xyz';
            xyz=xyz';
%             [mm,nn]=size(xyz);
%             for kk=1:mm,
%                 xyz(kk,:)=affine_trafo_point(xyz(kk,:),matrix);
%             end;
            model.labels(k).NOpos(:,1:3)=xyz(:,1:3);
        end;
    end;
end;
if isfield(model,'sites'),
    scans=numel(model.sites);
    for k=1:scans,
        for k0=1:length(model.sites{k}),
            if isfield(model.sites{k}(k0),'residue');
                residues=numel(model.sites{k}(k0).residue);
                for kk=1:residues,
                    if model.sites{k}(k0).residue(kk).indices(1)==snum, % transform only residues in current structure
                        matrix = transform{model.sites{k}(k0).residue(kk).indices(3)};
                        xyz=model.sites{k}(k0).residue(kk).NOpos(:,1:3);
                        [mm,nn]=size(xyz);
                        xyz=[xyz ones(mm,1)];
                        xyz=matrix*xyz';
                        xyz=xyz';
                        model.sites{k}(k0).residue(kk).NOpos(:,1:3)=xyz(:,1:3);
                    end;
                end;
            end;
        end;
    end;
end;
