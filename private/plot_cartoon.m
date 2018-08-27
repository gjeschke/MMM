function [gobjects,x,y,z]=plot_cartoon(indices)
% function gobjects=plot_cartoon(indices)
%
% Plots cartoon model for a range of residues and returns
% array of graphics objects handles
%
% the handle graphics objects and graphics settings are stored in a
% structure variable residue_graphics for that residue with fields
% .objects  handles to Matlab graphics objects
% .mode     atom graphics mode, 0 none, 1 tube, 2 helical, 3 sheet-like
% .color    RGB color vector
% .opaque   opaqueness (1-transparency) 
%
% indices   indices to identify structure, coordinate set, chain, and residue
%
% gobjects  vector of handles to graphics objects of all alternate locations
%
% G. Jeschke, 2014

global model
global graph_settings

Tension = 0;

gobjects=[];


[m,~] = size(indices);
nindices = indices;
np = 0;

% remove non-amino acid residues
for k = 1:m,
    info=model.structures{indices(k,1)}(indices(k,2)).residues{indices(k,3)}.info(indices(k,4));
    if info.type == 1,
        np = np+1;
        nindices(np,:) = indices(k,:);
    end;
end;
indices = nindices(1:np,:);
m = np;

snum = indices(1,1);
cnum = indices(1,2);

secpoints = zeros(m,4);
pairs = zeros(2,m);
ppoi = 0;
spoi = 0;
% Helices
if isfield(model.structures{snum}(cnum),'helix_defs'),
    helnum=length(model.structures{snum}(cnum).helix_defs);
    for k=1:helnum, 
        hrange=model.structures{snum}(cnum).helix_defs{k}.range;
        if min(indices(:,4))> hrange(1),
            hrange(1) = min(indices(:,4));
        end;
        if max(indices(:,4))< hrange(2),
            hrange(2) = max(indices(:,4));
        end;
        if hrange(2) > hrange(1),
            [objects,x,y,z,CA]=plot_sec_element(indices(k,1:3),hrange,-1);
            minx=min(min(x));
            maxx=max(max(x));
            miny=min(min(y));
            maxy=max(max(y));
            minz=min(min(z));
            maxz=max(max(z));
            rindices=indices(1:3);
            record_cartoon_object(objects,rindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            secpoints(spoi+1,1) = hrange(1);
            secpoints(spoi+1,2:4) = CA(1,:);
            secpoints(spoi+2,1) = hrange(2);
            secpoints(spoi+2,2:4) = CA(end,:);
            spoi = spoi + 2;
            pairs(:,ppoi+1) = hrange';
            ppoi = ppoi +1;
            gobjects=[gobjects objects];
        end;
    end;
end;
% Strands
if isfield(model.structures{snum}(cnum),'sheet_defs'),
    helnum=length(model.structures{snum}(cnum).sheet_defs);
    for k=1:helnum, 
        hrange=model.structures{snum}(cnum).sheet_defs{k}.range;
        if min(indices(:,4))> hrange(1),
            hrange(1) = min(indices(:,4));
        end;
        if max(indices(:,4))< hrange(2),
            hrange(2) = max(indices(:,4));
        end;
        if hrange(2) > hrange(1),
            [objects,x,y,z,CA]=plot_sec_element(indices(k,1:3),hrange,2);
            minx=min(min(x));
            maxx=max(max(x));
            miny=min(min(y));
            maxy=max(max(y));
            minz=min(min(z));
            maxz=max(max(z));
            rindices=indices(1:3);
            record_cartoon_object(objects,rindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            secpoints(spoi+1,1) = hrange(1);
            secpoints(spoi+1,2:4) = CA(1,:);
            secpoints(spoi+2,1) = hrange(2);
            secpoints(spoi+2,2:4) = CA(end,:);
            spoi = spoi + 2;
            pairs(:,ppoi+1) = hrange';
            ppoi = ppoi +1;
            gobjects=[gobjects objects];
        end;
    end;
end;
secpoints = secpoints(1:spoi,:);
[sorted,poi] = sort(secpoints(:,1));
secpoints = secpoints(poi,:);
pairs = pairs(:,1:ppoi);



% Loops
[ra,pa] = min(indices(:,4));
[re,pe] = max(indices(:,4));
if secpoints(1,1) ~= ra,
    CA = get_CA_coor(indices(pa,:));
    secpoints = [[ra,CA];secpoints];
end;
if secpoints(end,1) ~= re,
    CA = get_CA_coor(indices(pe,:));
    secpoints = [secpoints;[re,CA]];
end;

for k = 1:spoi-1,
    cpair = secpoints(k:k+1,1);
    if cpair(2)- cpair(1) > 1,
        diff = sum(abs(pairs - repmat(cpair,1,ppoi)));
        if min(diff>0),
            % fprintf(1,'Loop to be plotted between residues %i and %i.\n',secpoints(k,1),secpoints(k+1,1));
            CA = zeros(secpoints(k+1,1)-secpoints(k,1)+1,3);
            CA(1,:) = secpoints(k,2:4);
            CA(end,:) = secpoints(k+1,2:4);
            pp = 1;
            for kk = cpair(1)+1 : cpair(2)-1,
                CA(pp+1,:) = get_CA_coor([indices(1,1:3) kk]);
                pp = pp + 1;
            end;
            CA  = [CA(1,:); CA; CA(end,:)];
            [mc,~] = size(CA);
            guide_length=(mc-3)*graph_settings.spr+1;
            backbone=zeros(guide_length,3);
            for kk = 1:mc-3,
                basnum=(kk-1)*graph_settings.spr;
                backbone(basnum+1:basnum+graph_settings.spr+1,:)=cardinal_spline(CA(kk:kk+3,:),Tension,graph_settings.spr);
            end;
%             oldax = 1:mc;
%             newax = linspace(1,mc,guide_length);
%             px = polyfit(oldax,CA(:,1)',3);
%             xp = polyval(px,newax);
%             py = polyfit(oldax,CA(:,2)',3);
%             yp = polyval(py,newax);
%             pz = polyfit(oldax,CA(:,3)',3);
%             zp = polyval(pz,newax);
%             backbone = [xp;yp;zp];
            [x,y,z] = tubeplot(backbone',0.75*graph_settings.coil_radius);
            objects = surf(x,y,z);
            minx=min(min(x));
            maxx=max(max(x));
            miny=min(min(y));
            maxy=max(max(y));
            minz=min(min(z));
            maxz=max(max(z));
            rindices=indices(1:3);
            record_cartoon_object(objects,rindices,[(minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2,maxx-minx,maxy-miny,maxz-minz]);
            set(objects,'EdgeColor','none','FaceColor',graph_settings.coil_color,'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
            set(objects, 'CDataMapping','direct','AlphaDataMapping','none');
        end;
    end;
end;
    
if isempty(gobjects),
    message.error=1;
    message.text='Nothing to plot.';
end;

function CA_coor = get_CA_coor(indices)

cmadr = mk_address(indices(1:3));
resadr = sprintf('%s%i',cmadr,indices(4));
[~,CA_coor] = get_object([resadr '.CA'],'coor');
    
