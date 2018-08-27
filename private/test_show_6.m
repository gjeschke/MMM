figure(2); clf;
daspect([1 1 1])
view(3)
hold on;
set(gca,'FontSize',24);
hold on;

address='[2W8A](A)';

argin{1}='ribbon';
tic,
set_object('[2W8A](A)','show',{'ribbon'});
toc,

camlight 
lighting gouraud
material dull
axis off