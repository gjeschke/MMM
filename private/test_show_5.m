figure(1); clf;
daspect([1 1 1])
view(3)
hold on;
set(gca,'FontSize',24);
hold on;

address='[1QJP]';
% [indices,message]=resolve_address(address);
% 
% [m,n]=size(indices);

argin{1}='ribbon';
tic,
set_object(address,'show',argin);
toc,

camlight 
lighting gouraud
material dull
axis off