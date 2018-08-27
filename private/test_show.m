

figure(1); clf;
daspect([1 1 1])
view(49,26)
hold on;
set(gca,'FontSize',24);
hold on;

tic,
set_object('[NhaA](A)','show',{'ribbon'});
toc,

tic,
set_object('[NhaA](A)<H.IVa>','show',{'ball&stick'});
toc,

camlight 
lighting gouraud
material dull
