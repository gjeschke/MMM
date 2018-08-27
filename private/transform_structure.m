function transform_structure(snum,transform,bilayer)
% function transform_structure(snum,transform,bilayer)
%
% Perform a coordinate transform on the structure with number snum
% transform is an affine transformation matrix (4x4) or a cell array of such
% matrices
%
% transforms atom coordinates, label coordinates, and solvent accessible
% surfaces
%
% after the transform, the graphics may no longer be conssitent with
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

indices=snum;
if iscell(transform),
    matrix=eye(4);
    for k=1:length(transform),
        matrix=transform{k}*matrix;
    end;
else
    matrix=transform;
end;
p0=matrix(1:3,4)';

% set_object(indices,'transform',matrix);
if length(model.structures)>=snum,
    for k=1:length(model.structures{snum}),
        for kk=1:length(model.structures{snum}(k).xyz),            
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
if isfield(model,'SAS') && ~isempty(model.SAS),
    for k=1:length(model.SAS),
        if model.SAS(k).snum==snum,
            coor=[model.SAS(k).x' model.SAS(k).y' model.SAS(k).z' ones(size(model.SAS(k).x'))];
            coor=matrix*coor';
            model.SAS(k).x=coor(1,:);
            model.SAS(k).y=coor(2,:);
            model.SAS(k).z=coor(3,:);
        end;
    end;
end;
if isfield(model,'densities') && ~isempty(model.densities),
    for k=1:length(model.densities),
        if ~isfield(model.densities{k},'snum') || model.densities{k}.snum==snum,
            matrix2=matrix;
            matrix2(1:3,4)=[0;0;0];
            x=model.densities{k}.x;
            nx=length(x);
            dx=(max(x)-min(x))/nx;
            y=model.densities{k}.y; 
            ny=length(y);
            dy=(max(y)-min(y))/ny;
            z=model.densities{k}.z;
            nz=length(z);
            dz=(max(z)-min(z))/nz;
            centroid=[mean(x);mean(y);mean(z);1];
            centroid=matrix*centroid;
            add_msg_board('Transforming a density cube');
            add_msg_board('(this may take quite a lot of time)');
            drawnow;
            [cube,new_matrix] = affine_cube(model.densities{k}.cube, inv(matrix2));
            [ny,nx,nz]=size(cube);
            model.densities{k}.x=linspace(centroid(1)-(nx-1)*dx/2,centroid(1)+(nx-1)*dx/2,nx);
            model.densities{k}.y=linspace(centroid(2)-(ny-1)*dy/2,centroid(2)+(ny-1)*dy/2,ny);
            model.densities{k}.z=linspace(centroid(3)-(nz-1)*dz/2,centroid(3)+(nz-1)*dz/2,nz);
            model.densities{k}.cube=cube;
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
if isfield(model,'motion'),
    for k=1:length(model.motion),
        id=model.motion(k).snum;
        if id==snum,
            if model.motion(k).transform,
                xyz=model.motion(k).template;
                [mm,nn]=size(xyz);
                xyz=[xyz ones(mm,1)];
                xyz=matrix*xyz';
                xyz=xyz';
                model.motion(k).template=xyz(:,1:3);
                xyz=model.motion(k).target;
                [mm,nn]=size(xyz);
                xyz=[xyz ones(mm,1)];
                xyz=matrix*xyz';
                xyz=xyz';
                model.motion(k).target=xyz(:,1:3);
            else
                model.motion(k).transform=1;
            end;
        end;
    end;
end;
