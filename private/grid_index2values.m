function grid_selectVal=grid_index2values(grid_select,grid,homhet)
% homhet - flag: 1 corresponds to the homodeimer grid with symmetry axis
%                2 corresponds to heterodimer
 

if nargin<3
    homhet=2; % make heterodimer a default case
end

[n,m]=size(grid_select);
grid_selectVal=zeros(n,m);

switch homhet
    case 1 % homodimer
         for ee=1:n
            numbers=[grid_select(ee,1),grid_select(ee,2),grid_select(ee,3),grid_select(ee,4),grid_select(ee,5)];
            values=[grid.alpha(numbers(1)),grid.beta(numbers(2)),grid.x(numbers(3)),grid.y(numbers(4)),numbers(5)];
            grid_selectVal(ee,:)=values;
        %     disp(grid_selectVal);
        %     keyboard
        end
    case 2  % heterodimer
        for ee=1:n
            numbers=[grid_select(ee,1),grid_select(ee,2),grid_select(ee,3),grid_select(ee,4),grid_select(ee,5),grid_select(ee,6),grid_select(ee,7)];
            values=[grid.alpha(numbers(1)),grid.beta(numbers(2)),grid.gama(numbers(3)),grid.x(numbers(4)),grid.y(numbers(5)),grid.z(numbers(6)),numbers(7)];
            grid_selectVal(ee,:)=values;
        %     disp(grid_selectVal);
        %     keyboard
        end
end

        