figure(2); clf;
daspect([1 1 1])
view(62,-26)
hold on;
set(gca,'FontSize',24);
hold on;

tic,
set_object('[NMR](A)','show',{'coil'});
set_object('[NMR](A){:}"Ile"','show',{'stick'});
toc,

camlight 
lighting gouraud
material dull
axis off