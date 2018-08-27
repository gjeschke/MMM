figure(1); clf;
daspect([1 1 1])
view(49,26)
hold on;
set(gca,'FontSize',24);
hold on;

tic,
set_object('[T4L](A){:}','show',{'ribbon'});
toc,

camlight 
lighting gouraud
material dull
axis off

tic,
set_object('[T4L](A)"AZI"','show',{'ball&stick'});
set_object('[T4L](A)"HED"','show',{'ball&stick'});
set_object('[T4L](A)."Cl"','show',{'space-filling'});
set_object('[T4L]"R1A"','show',{'stick'});
set_object('[T4L]"R1A".N1','show',{'space-filling'});
set_object('[T4L]"R1A".O1:A','show',{'space-filling'});
toc,

view3d(gcf,'rot');
