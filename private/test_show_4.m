figure(3); clf;
daspect([1 1 1])
view(-11,12)
hold on;
set(gca,'FontSize',24);
hold on;

tic,
set_object('[LHCIIb](A)','show',{'coil'});
toc,

camlight 
lighting gouraud
material dull
axis off