function wr_docking_session(fid,dockStat)

% writes docking statistics to ASCII result file

fitORdock=dockStat.fitORdock;
acttime=dockStat.acttime;
constrFile=dockStat.constrFile;
dockcase=dockStat.dockcase;
templatesStr=dockStat.pdbtemplates;
modelQual(:,1)=dockStat.distBest;
modelQual(:,2:4)=dockStat.constraints.pairs_exp(:,3:5);
modelQualtitle=sprintf('%s%s%s%s', '   Model  ','  Exp ','  sig(exp)  ',' num'); % carefull: spaces are matched with the numerical output below!!!!

if fitORdock
    dockproc1='Docking was performed solely by fitting (Nelder-Mead simplex method, Matlab default).';
    dockproc2='This procedure is proone to get stuck in a local minimum especially if the error surface is rough.';
    dockproc3='The latter is practically always the case while docking pretine fragments';
    dockproc4='Independent estimation of a global minimum is higly recommented (with grid search for example).';
    dockprocstr=[dockproc1 dockproc2 dockproc3 dockproc4];
    fitst=dockStat.fitst;
    bestfitonly=dockStat.bestfitonly;
else
    dockproc1='Docking was performed using a grid search algorithm. ';
    if isfield(dockStat,'stackSize')
        dockstack=sprintf('The grid was considered large, therefore a stackind mode with stack size of %i was chosen,',dockStat.stackSize);
        dockproc1=[dockproc1 dockstack];
    end
    dockprocstr=dockproc1;
    gridin=dockStat.gridin;
    gridout=dockStat.gridout;
    if isfield(dockStat,'fitpostdock')
        dockproc2='Subsequently, fitting was performed using the best grid search model on input.';
        dockprocstr=[dockproc1 dockproc2];
        fitpostdock=dockStat.fitpostdock;
    end
end

switch dockcase
    case '1  1'
        dockcase_str1='The case was treated as homooligomer during docking (4 degrees of freedom. ';
        dockcase_str2='Coordinates of the first molecule are obtained by applying Euler rotation with [alpha, beta, 0] and no translation ([0 0 0]) to the initial structure. ';
        dockcase_str3='Coordinates of the second molecule (dimer) or all the other molecules in the complex are obtained either by applying Euler rotation [alpha beta 360/multiplicity] to the initial structure followed by the [x y 0] translation ';
        dockcase_str4='or by applying Euler rotation with [360/multiplicity 0 0] to the already transformed first molecule in the complex followed by the [x y 0] translation. ';
        dockcase_str5='Atom coordinates of the initial structure are meant to be centered first.';
        dockcasestr=[dockcase_str1 dockcase_str2 dockcase_str3 dockcase_str4 dockcase_str5];
        modtitle=sprintf('%s%s%s%s%s%s%s', '   Alpha  ',' Beta','  Gama(FIXED) ',' x    ','   y  ',' z(FIXED)','  rms'); % carefull: spaces are matched with the numerical output below!!!!
    case '1  0'
        dockcase_str1='This genuinely homodimer case was treated as heterodimer duting docking (6 degrees of freedom).';
        dockcase_str2='Complex is obtained by keeping the first molecule fixed and applying Euler rotation [alpha beta gama] to the second molecule followed by translation [x y z].';
        dockcase_str3='Atom coordinates of both molecules are meant to be centered first.';
        dockcasestr=[dockcase_str1 dockcase_str2 dockcase_str3];
        modtitle=sprintf('%s%s%s%s%s%s%s', '   Alpha ','  Beta ','   Gama    ','   x    ','  y  ','     z   ','   rms'); % carefull: spaces are matched with the numerical output below!!!!
    case '2  0'
        dockcase_str1='The case was treated as a heterodimer during docking (6 degrees of freedom).';
        dockcase_str2='Complex is obtained by keeping the first molecule fixed and applying Euler rotation [alpha beta gama] to the second molecule followed by translation [x y z].';
        dockcase_str3='Atom coordinates of both molecules are meant to be centered first.';
        dockcasestr=[dockcase_str1 dockcase_str2 dockcase_str3];
        modtitle=sprintf('%s%s%s%s%s%s%s', ' Alpha  ','   Beta ','   Gama    ','   x    ','  y  ','     z   ','   rms'); % carefull: spaces are matched with the numerical output below!!!!
end


fprintf(fid,'# Docking session performed on %s\n',acttime);
fprintf(fid,'# Docked templates: %s\n', templatesStr);
fprintf(fid,'# Docking constraints file: %s\n', constrFile); fprintf(fid,' \n');
fprintf(fid,'# Docking case: %s\n',dockcasestr); fprintf(fid,' \n');
fprintf(fid,'# Docking procedure: %s\n',dockprocstr); fprintf(fid,' \n');
% put grid here if grid search was performed

if fitORdock
    fprintf(fid,'# Initial fit values:\n'); fprintf(fid,' \n');
    fprintf(fid,'%s\n',modtitle);
    fprintf(fid,'%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n',fitst);
    fprintf(fid,'# Bestfit docking model: \n');
    fprintf(fid,'%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.4f',bestfitonly(1),bestfitonly(2),bestfitonly(3),bestfitonly(4),bestfitonly(5),bestfitonly(6),bestfitonly(7)); fprintf(fid,' \n');
    fprintf(fid,' \n');
    fprintf(fid,'# Quality of the best model, mean distances:\n'); % make distances here!!!
    fprintf(fid,'%s\n',modelQualtitle);
    for ii=1:length(modelQual(:,1))
        fprintf(fid,'%8.2f%8.2f%8.2f%8.2f\n',modelQual(ii,1),modelQual(ii,2),modelQual(ii,3),modelQual(ii,4));
    end
else
    fprintf(fid,'# Initial grid: \n'); % with gridin
    infields=fieldnames(gridin);
    for ii=1:numel(infields)
        if ~strcmp(infields{ii},'grSize')
            fieldvals=getfield(gridin,infields{ii});
            fprintf(fid,'%s: from  %5.2f to %5.2f with %i steps; \n',infields{ii},fieldvals(1),fieldvals(length(fieldvals)),numel(fieldvals));
        else
            fieldvals=getfield(gridin,infields{ii});
            fprintf(fid,'Grid size: %i;\n',fieldvals);
        end
    end
    
    fprintf(fid,' \n'); 
    fprintf(fid,'# Best grid model: \n');
    fprintf(fid,'%s\n',modtitle);
    fprintf(fid,'%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.4f\n', gridout(1,1),gridout(1,2),gridout(1,3),gridout(1,4),gridout(1,5),gridout(1,6),gridout(1,7));
    if isfield(dockStat,'fitpostdock')
        fprintf(fid,'# Best grid model refined by fitting: \n');
        fprintf(fid,'%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.4f\n',fitpostdock(1),fitpostdock(2),fitpostdock(3),fitpostdock(4),fitpostdock(5),fitpostdock(6),fitpostdock(7));
        fprintf(fid,'# Quality of the best model, mean distances:\n'); % make distances here!!!
    end
    fprintf(fid,' \n');
    fprintf(fid,'# Quality of the best model, mean distances:\n'); % make distances here!!!
    fprintf(fid,'%s\n',modelQualtitle);
    for ii=1:length(modelQual(:,1))
        fprintf(fid,'%8.2f%8.2f%8.2f%8.2f\n',modelQual(ii,1),modelQual(ii,2),modelQual(ii,3),modelQual(ii,4));
    end
    fprintf(fid,' \n'); 
    fprintf(fid,'# Best 20 grid search models: \n');
    for ii=1:20
        fprintf(fid,'%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.4f\n',gridout(ii,1),gridout(ii,2),gridout(ii,3),gridout(ii,4),gridout(ii,5),gridout(ii,6),gridout(ii,7));
    end
end
