figure(2);
clf;
hold on;
t=linspace(0,2*pi,50);
[x,y,z]=tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
ha=surf(x(:,1:2),y(:,1:2),z(:,1:2));
set(ha,'FaceColor','g','EdgeColor','none','FaceLighting','gouraud');
he=surf(x(:,51:52),y(:,51:52),z(:,51:52));
set(he,'FaceColor','g','EdgeColor','none','FaceLighting','gouraud');
for k=1:10,
    bas=(k-1)*5;
    h=surf(x(:,bas+2:bas+7),y(:,bas+2:bas+7),z(:,bas+2:bas+7));
    set(h,'FaceColor','g','EdgeColor','none','FaceLighting','gouraud');
end;
camlight
axis equal