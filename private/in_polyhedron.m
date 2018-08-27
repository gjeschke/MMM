function in_flag=in_polyhedron(point,tri,vertices,tnorm)
% in_flag=in_polyhedron(point,tri,vertices,tnorm)
%
% Determines whether the point with coordinates coor is inside the
% polyhedron with triangular faces indexed by tri
%
% requires intersect_triangle.m
%
% Input:
% point     1*3 xyz coordinates of the point
% tri       mt*3 list of index triples that specifies the triangular faces
%           of the polyhedron
% vertices  mv*3 list of vertex xyz coordinates of the polyhedron
% tnorm     mt*3 triangle face normal vectors, not necessarily normalized,
%           provided for speed, for degenerate triangles, the face normal
%           vector must be [0,0,0] (calling routine is responsible)
%
% Output:
% in_flag   flag that is 1 for a point inside and 0 for a point outside the
%           polyhedron
%
% G. Jeschke, 2010

[mt,nt]=size(tri);
outside=1.1*max(vertices); % point that is outside polyhedron

in_hit=zeros(mt,3); %list of intersections inside triangle
ihpoi=0;
vert_hit=zeros(mt,3); % list of vertex hits
hit_vert=zeros(1,mt);
vhpoi=0; % number of vertex hits
edge_hit=zeros(mt,3); % list of edge hits
hit_edge=zeros(mt,2);
ehpoi=0; % number of edge hits
full_hit=zeros(mt,3); % list of full hits (ray inside plane)
fhpoi=0; % number of full hits
ray=[point;outside]; % definition of the ray that is cast (to infinity)
degenerate=sum(abs(tnorm),2);

for k=1:mt, % loop over all triangles
    if degenerate(k)>eps, % exclude degenerate triangles
        triind=tri(k,:);
        triangle=vertices(triind,:); % triangle vertex coordinates
        [intersect,hit]=intersect_triangle(ray,triangle,tnorm(k,:)); % test whether ray intersects triangle
        if intersect==1 || intersect==2,
            switch numel(hit)
                case 0 % inside hit
                    ihpoi=ihpoi+1;
                    in_hit(ihpoi,:)=triind;
                case 1 % vertex hit
                    vhpoi=vhpoi+1;
                    vert_hit(vhpoi,:)=triind;
                    hit_vert(vhpoi)=hit;
                case 2 % edge hit
                    ehpoi=ehpoi+1;
                    edge_hit(ehpoi,:)=triind;
                    hit_edge(ehpoi,:)=hit;
                case 3 % full hit
                    fhpoi=fhpoi+1;
                    full_hit(fhpoi,:)=triind;
            end;
        end;
    end;
end;

% avoid multiple counting of full, vertex, and edge hits
% this is done very defensively, in principle some check (vertex vs. edge,
% inside vs. vertex and edge) should never return a hit (except for
% numerical imprecision)

intersections=fhpoi; % number of intersections is at least as large as number of full hits
if ehpoi>0, % there are edge hits
    for k=1:ehpoi,
        triind=edge_hit(k,:);
        hit=false;
        if fhpoi>0,
            for kk=1:fhpoi, % check all full hits
                for kkk=1:3, % check all vertices of current triangle
                    if sum(full_hit(kk,:)==triind(kkk)), % if it belongs to a fully hit triangle,
                        hit=true; % avoid double counting
                        break
                    end;
                end;
                if hit, break; end;
            end;
        end;
        if ~hit,
            intersections=intersections+1;
        end;
    end;
end;
if vhpoi>0, % there are vertex hits
    for k=1:vhpoi,
        triind=vert_hit(k,:);
        hit=false;
        if fhpoi>0,
            for kk=1:fhpoi, % check all full hits
                for kkk=1:3, % check all vertices of current triangle
                    if sum(full_hit(kk,:)==triind(kkk)), % if it belongs to a fully hit triangle,
                        hit=true; % avoid double counting
                        break
                    end;
                end;
                if hit, break; end;
            end;
        end;
        if ~hit && ehpoi>0,
            for kk=1:ehpoi, % check all edge hits
                for kkk=1:3, % check all vertices of current triangle
                    if sum(hit_edge(kk,:)==triind(kkk)), % if it belongs to a fully hit triangle,
                        hit=true; % avoid double counting
                        break
                    end;
                end;
                if hit, break; end;
            end;
        end;        
        if ~hit,
            intersections=intersections+1;
        end;
    end;
end;
if ihpoi>0, % there are inside hits
    for k=1:ihpoi,
        triind=in_hit(k,:);
        hit=false;
        if fhpoi>0,
            for kk=1:fhpoi, % check all full hits
                for kkk=1:3, % check all vertices of current triangle
                    if sum(full_hit(kk,:)==triind(kkk)), % if it belongs to a fully hit triangle,
                        hit=true; % avoid double counting
                        break
                    end;
                end;
                if hit, break; end;
            end;
        end;
        if ~hit && ehpoi>0,
            for kk=1:ehpoi, % check all edge hits
                for kkk=1:3, % check all vertices of current triangle
                    if sum(hit_edge(kk,:)==triind(kkk)), % if it belongs to a fully hit triangle,
                        hit=true; % avoid double counting
                        break
                    end;
                end;
                if hit, break; end;
            end;
        end;        
        if ~hit && vhpoi>0,
            for kk=1:vhpoi, % check all vertex hits
                if sum(triind==hit_vert(kk)), % if hit vertex belongs to current triangle,
                    hit=true; % avoid double counting
                end;
                if hit, break; end;
            end;
        end;        
        if ~hit,
            intersections=intersections+1;
        end;
    end;
end;

in_flag=mod(intersections,2); % odd number of intersections => inside, even number => outside
    