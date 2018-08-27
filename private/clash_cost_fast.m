function cost = clash_cost_fast(a,b,tol,display)

cost = 0;

if ~exist('tol','var') || isempty(tol)
    tol = 1.5;
end

if ~exist('display','var') || isempty(display)
    display = false;
end

ka = convhulln(a);
kb = convhulln(b);

[na,~] = size(ka);
[kpa,ar] = reducepatch(ka,a,na);
[nb,~] = size(kb);
[kpb,br] = reducepatch(kb,b,nb);
    
inb = inhull(ar,br,kpb,tol);
p1 = ar(inb,:);
ina = inhull(br,ar,kpa,tol);
p2 = br(ina,:);
p = [p1;p2];

if display
    figure(1); clf;
    h1=trisurf(ka,a(:,1),a(:,2),a(:,3));
    set(h1, 'FaceColor', [0,0,1], 'EdgeColor', 'none', 'FaceAlpha',0.5,'FaceLighting','gouraud','Clipping','off');
    set(h1, 'CDataMapping','direct','AlphaDataMapping','none');
    axis off
    hold on
    h2=trisurf(kb,b(:,1),b(:,2),b(:,3));
    set(h2, 'FaceColor', [0,1,0], 'EdgeColor', 'none', 'FaceAlpha',0.5,'FaceLighting','gouraud','Clipping','off');
    set(h2, 'CDataMapping','direct','AlphaDataMapping','none');
end

[np,~] = size(p);
if np > 3
    if display
        [kp,cost] = convhulln(p);
        h3=trisurf(kp,p(:,1),p(:,2),p(:,3));
        set(h3, 'FaceColor', [1,0,0], 'EdgeColor', 'none', 'FaceAlpha',0.5,'FaceLighting','gouraud','Clipping','off');
        set(h3, 'CDataMapping','direct','AlphaDataMapping','none');
    else
        [~,cost] = convhulln(p);
    end
end
