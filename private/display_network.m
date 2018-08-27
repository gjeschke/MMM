function display_network(handles)

global general
global hMain
global model

if hMain.virgin,
    hMain.virgin=0;
    % initialize display
    axes(handles.axes_model);
    cla reset;
    axis equal
    axis off
    set(gca,'Clipping','off');
    set(gcf,'Renderer','opengl');
    hold on
    hMain.camlight=camlight;
    guidata(handles.axes_model,hMain);
end;

node_radius = 2;
constraint_width = 4;
node_color = [0,0,0.25];
cloud_color = [0.5,0,0];
maxtime = 10*60*60;
mcn = 1000; % number of Monte Carlo trials for sparse distance matrices 
grid_points = 151;
grid_size = 60;
opaqueness = 0.5;

mypath=pwd;

cd(general.restraint_files);

[FileName,PathName] = uigetfile({'*.dat'},'Select constraints file');
if isequal(FileName,0) || isequal(PathName,0),
    add_msg_board('Loading of constraints canceled by user.');
    message.error=5;
    message.text='Cancelled.';
    return;
end;
reset_user_paths(PathName);
general.restraint_files=PathName;
restraints=fullfile(PathName,FileName);

restraint_data=rd_restraints(restraints);

cd(PathName);

n = restraint_data.ensemble;
if restraint_data.uncertainty > 0.5,
    d = restraint_data.uncertainty;
else
    d = 2;
end;

[dmat, sigrmat, ub, lb, types, assignment, ref_points, node_tags, nodes_known, info, cancelled] = process_network_restraints(restraint_data);

if cancelled,
    add_msg_board('ERROR. Processing of constraints cancelled.');
    return;
end;


display_mode = restraint_data.network(1).display_mode;

[np,~] = size(ref_points);
[nn,~] = size(dmat);

df = 3*nn-6;

% test if there is sufficient information
constrained = 0;
restrained = 0;
for k1 = 1:nn-1,
    for k2 = k1+1:nn,
        if types(k1,k2)==1,
            constrained = constrained + 1;
        end;
        if types(k1,k2)==4,
            restrained = restrained + 1;
        end;
    end
end;
add_msg_board(sprintf('There are %i degrees of freedom (%i nodes).',df,nn));
add_msg_board(sprintf('%i distances are constrained by a distribution.',constrained));
add_msg_board(sprintf('%i distances are restrained by lower and upper bounds.',restrained));
if constrained+restrained < df,
    add_msg_board('ERROR. The problem is underdetermined. Quitting.');
    return;
end;
if constrained == nn*(nn-1)/2,
    add_msg_board('Problem is fully determined by constraints. Straight solution.');
else
    add_msg_board('Problem is not fully determined by constraints. Monte Carlo solution.');
end;


% index vector of reference points in full node vector
indvec = zeros(1,np);
for k = 1:nn,
    if nodes_known(k),
        indvec(abs(nodes_known(k))) = k;
    end;
end;

all_known = ~sum(sum(abs(types - ones(nn)))); % are all distances known (type 1 constraints)?

tic;
if all_known,
    coor = dmat2coor(dmat);
    if isempty(coor),
        add_msg_board('ERROR: Distance matrix could not be embedded.');
        return;
    end;
    [coor,err]=bound_refiner(coor,lb,ub);
    if err,
        add_msg_board('ERROR: Network is inconsistent with constraints.');
        return;
    end;
else
    hfig = gcf;
    set(hfig,'Pointer','watch');
    drawnow;
    mcoor = zeros(nn,3);
    k = 0;
    trials = 0;
    runtime = 0;
    eerr = 0;
    berr = 0;
    merr = 0;
    while k < mcn && runtime < maxtime,
        runtime = toc;
        trials = trials + 1;
        [cdmat,err] = network_metrize(lb,ub,dmat,sigrmat,types);
        if ~err,
            coor = dmat2coor(cdmat);
            if ~isempty(coor),
                [coor,err]=bound_refiner(coor,lb,ub);
                if ~err,
                    if k > 1,
                        coor = fit_partial_network(coor0, coor);
                    else
                        coor0 = coor;
                    end;
                    mcoor = mcoor + coor;
                    k = k + 1;
                else
                    berr = berr + 1;
                    % fprintf(1,'Trial %i with failed bound fulfilment.\n',trials);
                end;
            else
                eerr = eerr +1;
                % fprintf(1,'Trial %i with failed embedding.\n',trials);
            end;
            % fprintf(1,'Trial %i successful. %i models.\n',trials,k);
        else
            merr = merr +1;
            % fprintf(1,'Trial %i with failed metrization.\n',trials);
        end;
    end;
    if k == 0,
        add_msg_board(sprintf('ERROR: Distance matrix could not be embedded in any of the %i trials.',trials));
        return;
    else
        coor = mcoor/k;
        runtime = toc;
        add_msg_board(sprintf('Needed %i trials to generate %i network coordinate sets (%i failed metrizations, %i failed embeddings, %i structures inconstent with bounds).',trials,k,merr,eerr,berr));
        add_msg_board(sprintf('Time needed for distance matrix geometry: %5.1f s',runtime));
    end;
end;
drawnow;

if np>=3,
    [ncoor, rmsd, transmat, reflect] = fit_partial_network(ref_points, coor, indvec);
    add_msg_board(sprintf('Coordinate rmsd: %6.4f Å',rmsd));
    if reflect,
        add_msg_board('Network points were reflected.');
    end;
    % update reference point coordinates
    for k = 1:np,
        ref_points(k,:) = ncoor(indvec(k),:);
    end;
    reference_mode = true;
else
    ncoor = coor;
    reference_mode = false;
end;

n_ref = np;

if strcmpi(display_mode,'full') || strcmpi(display_mode,'nodes'),
    for k1 = 1:nn,
        if nodes_known(k1),
            coor = ncoor(k1,:);
            [x,y,z,t]=point2trisphere(coor,node_radius,2);
            obj=trisurf(t,x,y,z);
            xyz=[coor(1),coor(2),coor(3),2*node_radius,2*node_radius,2*node_radius];
    %         indices = 0;
    %         record_object(obj,indices,xyz);
            set(obj, 'FaceColor', node_color, 'EdgeColor', 'none', 'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
        end;
    end;
end;

if strcmpi(display_mode,'full') || strcmpi(display_mode,'constraints'),
    for k1 = 1:nn-1,
        for k2= k1+1:nn,
            if nodes_known(k1) && nodes_known(k2) && types(k1,k2)==1,
                coor1 = ncoor(k1,:);
                coor2 = ncoor(k2,:);
                plot3([coor1(1),coor2(1)],[coor1(2),coor2(2)],[coor1(3),coor2(3)],'Color',node_color,'LineWidth',constraint_width);
            end;
        end;
    end;
end;

for k = 1:length(nodes_known),
    if ~nodes_known(k),
        tag = id2tag(k,node_tags);
        xyz = ncoor(k,:);
        add_msg_board(sprintf('Mean position of %s = (%4.2f, %4.2f, %4.2f) Å', tag, xyz(1), xyz(2), xyz(3)));
    end;
end;

hfig = gcf;
set(hfig,'Pointer','watch');
nx=grid_points;
ny=nx;
nz=nx;
if isfield(model,'density_tags'),
    if isfield(model,'densities'),
        poi=length(model.densities);
        if poi==0;
            model=rmfield(model,'densities');
             model.density_tags=':';
        end;
    end;
else
    poi=0;
    model.density_tags=':';
end;

alpoi = 0;
all_kn = zeros(1,nn);
for kn= 1:nn,
    if ~nodes_known(kn),
        rn = 0;
        rb = 0;
        ref_points = zeros(nn,3);
        bound_points = zeros(nn,3);
        dist = zeros(nn,1);
        sigr = dist;
        rlb = dist;
        rub = dist;
        assign_vec = dist;
        for kr = 1:nn,
            if kr ~= kn 
                if types(kn,kr) == 1,
                    rn = rn+1;
                    ref_points(rn,:) = ncoor(kr,:);
                    dist(rn) = dmat(kn,kr);
                    sigr(rn) = sigrmat(kn,kr);
                    assign_vec(rn) = assignment(kn,kr);
                else
                    rb = rb + 1;
                    bound_points(rb,:) = ncoor(kr,:);
                    rlb(rb) = lb(kn,kr);
                    rub(rb) = ub(kn,kr);
                end;
            end;
        end;
        ref_points = ref_points(1:rn,:);
        dist = dist(1:rn);
        sigr = sigr(1:rn);
        assign_vec = assign_vec(1:rn);
        bound_points = ref_points(1:rb,:);
        rlb = rlb(1:rb);
        rub = rub(1:rb);
        np = rn;
        tag=id2tag(kn,node_tags);
        if rn < 3,
            add_msg_board(sprintf('Warning: Node %s specified by only %i distances.',tag,rn));
            add_msg_board('Computation of probability density will be skipped for this node.');
            continue;
        end;
        most_probable = [0,0,0];
        max_probability = 0;
        add_msg_board(sprintf('Computing probability density for node %s. Please be patient.',tag));
        xyz = ncoor(kn,:);
        cube = zeros(ny,nx,nz);
        x = linspace(xyz(1) - grid_size/2, xyz(1) + grid_size/2, nx);
        y = linspace(xyz(2) - grid_size/2, xyz(2) + grid_size/2, ny);
        z = linspace(xyz(3) - grid_size/2, xyz(3) + grid_size/2, nz);
        fulldistr = 0;
        for kk = 1:length(info),
            fulldistr = fulldistr + info(kk).fulldistr;
        end;
        if fulldistr, % this code duplication saves a lot of time
            if rb>0,
                for kx = 1:nx,
                    for ky = 1:ny,
                        for kz = 1:nz,
                            coor = [x(kx),y(ky),z(kz)];
                            diff = ref_points - reprowvector(coor,np);
                            cdist = sqrt(sum(diff.^2,2));
                            prob = get_probabilities(cdist,dist,sigr,info,assign_vec);
                            diff = bound_points - reprowvector(coor,rb);
                            cdist = sqrt(sum(diff.^2,2));
                            if sum(cdist<rlb), prob = 0; end;
                            if sum(cdist>rub), prob = 0; end;
                            cube(ky,kx,kz) = prod(prob);
                            if cube(ky,kx,kz) > max_probability,
                                max_probability = cube(ky,kx,kz);
                                most_probable = coor;
                            end;
                        end;
                    end;
                end;
            else
                for kx = 1:nx,
                    for ky = 1:ny,
                        for kz = 1:nz,
                            coor = [x(kx),y(ky),z(kz)];
                            diff = ref_points - reprowvector(coor,np);
                            cdist = sqrt(sum(diff.^2,2));
                            prob = get_probabilities(cdist,dist,sigr,info,assign_vec);
                            cube(ky,kx,kz) = prod(prob);
                            if cube(ky,kx,kz) > max_probability,
                                max_probability = cube(ky,kx,kz);
                                most_probable = coor;
                            end;
                        end;
                    end;
                end;
            end;
        elseif rb>0
            for kx = 1:nx,
                for ky = 1:ny,
                    for kz = 1:nz,
                        coor = [x(kx),y(ky),z(kz)];
                        diff = ref_points - reprowvector(coor,np); 
                        cdist = sqrt(sum(diff.^2,2));
                        earg = (cdist - dist)./(sqrt(2)*sigr);
                        prob = exp(-earg.^2);
                        diff = bound_points - reprowvector(coor,rb);
                        cdist = sqrt(sum(diff.^2,2));
                        if sum(cdist<rlb), prob = 0; end;
                        if sum(cdist>rub), prob = 0; end;
                        cube(ky,kx,kz) = prod(prob);
                        if cube(ky,kx,kz) > max_probability,
                            max_probability = cube(ky,kx,kz);
                            most_probable = coor;
                        end;
                    end;
                end;
            end;
        else
            for kx = 1:nx,
                for ky = 1:ny,
                    for kz = 1:nz,
                        coor = [x(kx),y(ky),z(kz)];
                        diff = ref_points - reprowvector(coor,np); 
                        cdist = sqrt(sum(diff.^2,2));
                        earg = (cdist - dist)./(sqrt(2)*sigr);
                        prob = exp(-earg.^2);
                        cube(ky,kx,kz) = prod(prob);
                        if cube(ky,kx,kz) > max_probability,
                            max_probability = cube(ky,kx,kz);
                            most_probable = coor;
                        end;
                    end;
                end;
            end;
        end;
        cube = cube/max_probability;
        poi=poi+1;
        model.densities{poi}.x=x;
        model.densities{poi}.y=y;
        model.densities{poi}.z=z;
        model.densities{poi}.tag=tag;
        model.densities{poi}.cube=cube;
        model.density_tags=sprintf('%s%s:',model.density_tags,tag);
        sdens=sum(sum(sum(cube)));
        level=restraint_data.network(1).probability;
        level0 = level;
        for kl=1:99,
            mask=(cube>=kl/100);
            test=sum(sum(sum(mask.*cube)));
            if test<=level0*sdens,
                level=kl/100;
                break;
            end;
        end;
        [xg,yg,zg]=meshgrid(x,y,z);
        cnorm = sum(sum(sum(cube)));
        xm=sum(sum(sum(xg.*cube)))/cnorm;
        ym=sum(sum(sum(yg.*cube)))/cnorm;
        zm=sum(sum(sum(zg.*cube)))/cnorm;
        FV = isosurface(xg,yg,zg,cube,level);
        vert = FV.vertices;
        xmin = min(vert(1,:));
        xmax = max(vert(1,:));
        ymin = min(vert(2,:));
        ymax = max(vert(2,:));
        zmin = min(vert(3,:));
        zmax = max(vert(3,:));
        p = patch(FV);
        set(p, 'FaceColor', cloud_color, 'EdgeColor', 'none','FaceAlpha',opaqueness,'FaceLighting','gouraud','Clipping','off');
        set(p, 'CDataMapping','direct','AlphaDataMapping','none');
        dg.gobjects=p;
        dg.tag=['DMG_' tag];
        dg.color=cloud_color;
        dg.level=level;
        dg.cog=[xm,ym,zm];
        dg.transparency=opaqueness;
        dg.active=true;
        dg.mlp=most_probable;
        add_msg_board(sprintf('Center of gravity of density cloud %s is (%4.1f, %4.1f, %4.1f) Å',dg.tag,dg.cog(1),dg.cog(2),dg.cog(3)));
        if isfield(model,'surfaces'),
            model.surfaces=[model.surfaces dg];
        else
            model.surfaces=dg;
        end;
        indices = zeros(1,6);
        indices(1) = -length(model.surfaces);
        record_object(p,indices,[(xmin+xmax)/2,(ymin+ymax)/2,(zmin+zmax)/2,xmax-xmin,ymax-ymin,zmax-zmin]);
        add_msg_board(sprintf('Highest location probability is %5.3f at [%4.1f; %4.1f, %4.1f].',max_probability,most_probable));
        m1 = max(max(cube(1,:,:)));
        m2 = max(max(cube(end,:,:)));
        m3 = max(max(cube(:,1,:)));
        m4 = max(max(cube(:,end,:)));
        m5 = max(max(cube(:,:,1)));
        m6 = max(max(cube(:,:,end)));
        mm = max(max(max(cube)));

        border_prop = max([m1,m2,m3,m4,m5,m6])/mm;

        add_msg_board(sprintf('Maximum relative probability at density cube border is %5.2f%%',100*border_prop));
        if strcmpi(display_mode,'full') || strcmpi(display_mode,'constraints'),
            for k2= 1:nn,
                if reference_mode,
                    if k2~=kn && nodes_known(k2) && types(kn,k2)==1,
                        coor1 = most_probable;
                        coor2 = ncoor(k2,:);
                        plot3([coor1(1),coor2(1)],[coor1(2),coor2(2)],[coor1(3),coor2(3)],'Color',node_color,'LineWidth',constraint_width);
                    end;
                    if k2~=kn && ~nodes_known(k2) && types(kn,k2)==1,
                        coor1 = most_probable;
                        coor2 = ncoor(k2,:);
                        plot3([coor1(1),coor2(1)],[coor1(2),coor2(2)],[coor1(3),coor2(3)],'Color',cloud_color,'LineWidth',constraint_width);
                    end;
                else
                    if k2~=kn && types(kn,k2) == 1,
                        coor1 = most_probable;
                        coor2 = ncoor(k2,:);
                        plot3([coor1(1),coor2(1)],[coor1(2),coor2(2)],[coor1(3),coor2(3)],'Color',cloud_color,'LineWidth',constraint_width);
                    end;
                end;
            end;
        end;
        locations = location_models(x,y,z,cube,n,d);
        if isfield(model,'surfaces'),
            model.surfaces=[model.surfaces dg];
        else
            model.surfaces=dg;
        end;
        objects = zeros(1,length(locations.p));
        for l = 1: length(locations.p),
            [xl,yl,zl,tl]=point2trisphere(locations.xyz(l,:),d/2,2);
            obj=trisurf(tl,xl,yl,zl);
            objects(l) = obj;
            set(obj, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'FaceAlpha',locations.p(l),'FaceLighting','gouraud','Clipping','off');
            set(obj, 'CDataMapping','direct','AlphaDataMapping','none');    
        end;
        alpoi = alpoi+1;
        all_locations{alpoi} = locations;
        all_objects{alpoi} = objects;
        all_tags{alpoi} = tag;
        all_kn(alpoi) = kn;
    end;
end;

camlookat(hMain.axes_model);

if isfield(restraint_data,'output') && alpoi>0,
    wr_locations(restraint_data.network(1).pid,all_locations,restraint_data.output,all_tags);
    if n_ref < 4,
        wr_locations(restraint_data.network(1).pid,all_locations,restraint_data.output,all_tags,true);
    end;
end;
set(hfig,'Pointer','arrow');

cd(mypath);

function prob = get_probabilities(cdist,dist,sigr,info,assign_vec)

prob = zeros(1,length(cdist));
for k=1:length(cdist),
    if info(assign_vec(k)).fulldistr,
        [~,poi] = min(abs(10*info(k).rax - cdist(k)));
        prob(k) = info(k).distr(poi);
    else
        earg = (cdist(k) - dist(k))/(sqrt(2)*sigr(k));
        prob(k) = exp(-earg^2);
    end;
end;


function [dist, sigr, ub, lb, types, assignment, ref_points, node_tags, nodes_known, info, cancelled] = process_network_restraints(restraints)

global model
global hMain

dist = [];
sigr = [];
ub = [];
lb = [];
types = [];
assignment = [];
ref_points = [];

min_lower = 3; % default lower bound of node-node distance
max_upper = 100; % default upper bound of node-node distance

if isfield(restraints,'lower'),
    min_lower = 10*restraints.lower; % nm to Å
end;
if isfield(restraints,'upper'),
    max_upper = 10*restraints.upper; % nm to Å
end;

cancelled=true;

info.probability=0.95;
node_tags = ':';
label_adr = ':';
label_indices = zeros(1000,4);
nodes_known = zeros(1,1000);
nodes_assign = zeros(1,1000);
n_nodes = 0;
n_labels = 0;

if ~isfield(restraints,'network'),
    return;
end;

ref_tags = ':';
if isfield(restraints,'reference'),
    nr = length(restraints.reference);
    ref_coor = zeros(nr,3);
    for k = 1:nr,
        ref_tags = [ref_tags restraints.reference(k).tag ':'];
        ref_coor(k,:) = restraints.reference(k).xyz;
    end;
end;

if isfield(model,'current_structure'),
    snum = model.current_structure;
    template = true;
    if isfield(restraints,'PDB'),
        if ~strcmpi(model.info{snum}.idCode,restraints.PDB),
            button = questdlg(sprintf('Constraint file specifies template %s, while current template is %s. Do you want to continue?',restraints.PDB,model.info{snum}.idCode),'Mismatch between templates','Yes','No','No');
            if strcmp(button,'No'),
                cancelled=true;
                return
            end;
        end;
    end;
else
    template = false;
end;

n_ref2 = 0;
% compile node and label list
for k=1:length(restraints.network),
    adr1=restraints.network(k).adr1;
    ni = tag2id(adr1,node_tags);
    if isempty(ni),
        n_nodes = n_nodes + 1;
        node_tags = [node_tags adr1 ':'];
        ind1 = resolve_address(adr1);
        r1 = tag2id(adr1,ref_tags);
        if ~isempty(r1),
            n_ref2 = n_ref2 + 1;
            nodes_known(n_nodes) = -n_ref2;
            nodes_assign(n_nodes) = r1;
        elseif ~isempty(ind1) && length(ind1)==4, % valid residue address
            n_labels = n_labels + 1;
            nodes_known(n_nodes) = n_labels;
            label_indices(n_labels,:) = ind1;
            label_adr = [label_adr adr1 ':'];
        end;
    end;
    adr2=restraints.network(k).adr2;
    ni = tag2id(adr2,node_tags);
    if isempty(ni),
        n_nodes = n_nodes + 1;
        node_tags = [node_tags adr2 ':'];
        ind2 = resolve_address(adr2);
        r2 = tag2id(adr2,ref_tags);
        if ~isempty(r2),
            n_ref2 = n_ref2 + 1;
            nodes_known(n_nodes) = -n_ref2;
            nodes_assign(n_nodes) = r2;
        elseif ~isempty(ind2) && length(ind2)==4, % valid residue address
            n_labels = n_labels + 1;
            nodes_known(n_nodes) = n_labels;
            label_indices(n_labels,:) = ind2;
            label_adr = [label_adr adr2 ':'];
        end;
    end;
end;

label_indices = label_indices(1:n_labels,:);
nodes_known = nodes_known(1:n_nodes);
nodes_assign = nodes_assign(1:n_nodes);

types = eye(n_nodes);
dist = zeros(n_nodes);
sigr = zeros(n_nodes);
assignment = zeros(n_nodes);
lb = min_lower*(ones(n_nodes) - eye(n_nodes));
ub = max_upper*(ones(n_nodes) - eye(n_nodes));
ref_points = zeros(n_labels + n_ref2,3);

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

% check whether label sites are already labelled
% identity of the label is checked
% labeling temperature is NOT checked
T_list=zeros(1,200);
if ~isempty(labels),
    lindices=zeros(length(labels),4);
    for k=1:length(labels),
        cindices=labels(k).indices;
        if ~isempty(cindices),
            lindices(k,:)=cindices;
        end;
    end;
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:n_labels,
        adr1 = id2tag(k,label_adr);
        ind1 = label_indices(k,:);
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,restraints.network(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.network(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}=restraints.network(k).label;
                T_list(poi)=restraints.network(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.locate(k).label,adr1));
            end;
        end;
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:n_labels,
        adr1 = id2tag(k,label_adr);
        found = false;
        for l=1:length(to_do_list),
            if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.network(k).label),
                found=true;
            end;
        end;
        if ~found,
            poi=poi+1;
            to_do_list{poi}=adr1;
            label_list{poi}=restraints.network(k).label;
            T_list(poi)=restraints.network(k).T;
            add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
        end;
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

if isfield(model,'sites'),
    labels = label_information(model.sites);
else
    labels = {};
end;

for k=1:length(restraints.network),
    adr1=restraints.network(k).adr1;
    ind1=resolve_address(adr1);
    k1 = tag2id(adr1,node_tags);
    r1 = tag2id(adr1,ref_tags);
    adr2=restraints.network(k).adr2;
    ind2=resolve_address(adr2);
    k2 = tag2id(adr2,node_tags);
    r2 = tag2id(adr2,ref_tags);
    r = 10*restraints.network(k).r;
    sr = 10*restraints.network(k).sigr;
    type = 1;
    if r < 0,
        if sr > 0
            type = 2;
        else
            type = 4;
        end;
    elseif sr < 0
        type = 3;
    end;
    types(k1,k2) = type;
    types(k2,k1) = type;
    assignment(k1,k2) = k;
    assignment(k2,k1) = k;
    switch type
        case 1 % proper constraint r, sigr
            dist(k1,k2) = r;
            dist(k2,k1) = r;
            lb(k1,k2) = r - 2*sr;
            if lb(k1,k2) < min_lower,
                lb(k1,k2) = min_lower;
            end;
            lb(k2,k1) = lb(k1,k2);
            ub(k1,k2) = r + 2*sr;
            if ub(k1,k2) > max_upper,
                ub(k1,k2) = max_upper;
            end;
            ub(k2,k1) = ub(k1,k2);
            sigr(k1,k2) = sr;
            sigr(k2,k1) = sr;
        case 2 % lower bound
            lb(k1,k2) = -r;
            lb(k2,k1) = -r;
            ub(k1,k2) = max_upper;
            ub(k2,k1) = max_upper;
        case 3 % upper bound
            lb(k1,k2) = min_lower;
            lb(k2,k1) = min_lower;
            ub(k1,k2) = -sr;
            ub(k2,k1) = -sr;
        case 4 % lower and upper bound
            lb(k1,k2) = -r;
            lb(k2,k1) = -r;
            ub(k1,k2) = -sr;
            ub(k2,k1) = -sr;
    end;          
    info(k).ind1=ind1;
    info(k).ind2=ind2;
    info(k).probability=restraints.network(k).probability;
    info(k).adr1=adr1;
    info(k).adr2=adr2;
    info(k).fulldistr=restraints.network(k).fulldistr;
    if restraints.network(k).fulldistr,
        data = load(restraints.network(k).file);
        info(k).rax = data(:,1);
        info(k).distr = data(:,2);
        info(k).distr = info(k).distr/max(info(k).distr);
    end;
    f1 = true;
    f2 = true;
    if ~isempty(r1),
        info(k).xyz1 = restraints.reference(r1).xyz;
        info(k).rmsd1 = restraints.reference(r1).rmsd;
        f1 = true;
    elseif ~isempty(ind1),
        f1=false;
        for l=1:n_labels,
            diff1=ind1-labels(l).indices;
            if sum(abs(diff1))==0,
                f1=true;
                info(k).xyz1=labels(l).xyz;
                info(k).rmsd1=labels(l).rmsd;
            end;
        end;
    end;
    if ~isempty(r2),
        info(k).xyz2 = restraints.reference(r2).xyz;
        info(k).rmsd2 = restraints.reference(r2).rmsd;
        f2 = true;
    elseif ~isempty(ind2),
        f2=false;
        for l=1:n_labels,
            diff2=ind2-labels(l).indices;
            if sum(abs(diff2))==0,
                f2=true;
                info(k).xyz2=labels(l).xyz;
                info(k).rmsd2=labels(l).rmsd;
            end;
        end;
    end;
    if ~f1 || ~f2,
        add_msg_board('ERROR: Automatic rotamer computation error.');
        add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
        return;
    end;
end;
poi = 0;
for k=1:n_nodes,
    if nodes_known(k)>0,
        poi = poi+1;
        ref_points(poi,:) = labels(nodes_known(k)).xyz;
    elseif nodes_known(k)<0,
        poi = poi+1;
        ref_points(poi,:) = restraints.reference(nodes_assign(k)).xyz;
    end;
end;
cancelled=false;

function labels=label_information(sites)

global model
global label_defs

poi=0;
for k0=1:length(sites),
    for k1=1:length(sites{k0}),
        for k=1:length(sites{k0}(k1).residue),
            poi=poi+1;
            labels(poi).indices=sites{k0}(k1).residue(k).indices;
            id=tag2id(sites{k0}(k1).residue(k).label,label_defs.restags);
            labels(poi).name=label_defs.residues(id).short_name;
            labels(poi).T=sites{k0}(k1).residue(k).T;
            NOpos=model.sites{k0}(k1).residue(k).NOpos;
            x=sum(NOpos(:,1).*NOpos(:,4));
            y=sum(NOpos(:,2).*NOpos(:,4));
            z=sum(NOpos(:,3).*NOpos(:,4));
            labels(poi).xyz=[x y z];
            labels(poi).rmsd=NOpos_rmsd(NOpos);
        end;
    end;
end;

function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function clash = check_clash(coor1,coor2)

% clash = check_clash(coor1,coor2)

% checks whether two sets of coordinates produce a mutual clash, based on
% carbon Lennard-Jones potential and a 0.8 forgive factor

global LJ

conv_factor=(4.185*1000); % conversion to SI, - if working per mol

[m1,~] = size(coor1);
[m2,~] = size(coor2);


allowed = 5;
forgive = 0.9; % get softening factor for the LJ potential
full_clash = 0.05; % any pair of atoms closer than 0.05 Å immediately counts as clash

RT = 3*298*8314/2;
max_energy = allowed*RT;

a = coor1;
b = coor2;
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));
if min(abs(pair_dist)) < full_clash,
    clash = true; return;
end

Rmin2_i0 = LJ.Rmin2(6*ones(m1,1));
eps_i = LJ.eps(6*ones(m1,1))*conv_factor;
Rmin2_j0 = LJ.Rmin2(6*ones(m2,1));
eps_j = LJ.eps(6*ones(m2,1))*conv_factor;

eps = sqrt(eps_i(:)*eps_j);
Rmin2_i = Rmin2_i0*forgive; % soften LJ with the forgive factor
Rmin2_j = Rmin2_j0*forgive; % soften LJ with the forgive factor
Rmin = repmat(Rmin2_i(:),1,m2) + repmat(Rmin2_j,m1,1);
            % ATTENTION! Rmin2_i - is actually Rmin_i/2
q = (Rmin./pair_dist).^2;
q = q.*q.*q;
Energy = eps.*(q.^2-2*q);
Energy(pair_dist>10) = 0;
TotalEnergy = sum(sum(Energy));

if TotalEnergy > max_energy,
    clash = true;
else
    clash = false;
end;

function locations = location_models(x,y,z,cube,n,d)
% Represents a probability density cube by n high probability points 
% the points have a minimum distance d from each other
% with this constraint, the next point always corresponds to the highest
% probability point
%
% x,y,x     axes of the probability density cube
% cube      probability density cube
% n         number of points to be generated
% d         minimum distance between generated points
%
% location  structure of high probability locations with fields:
%           .xyz    location coordinate
%           .p      relative probability

nleft = n;

locations.xyz = zeros(n,3);
locations.p = zeros(n,1);

d2 = d^2;
dx = x(2)-x(1);
nx = floor(d/dx);
dy = y(2)-y(1);
ny = floor(d/dy);
dz = z(2)-z(1);
nz = floor(d/dz);
[my,mx,mz] = size(cube);

cube = cube/max(max(max(cube))); % probability normalization

l = 0;
while nleft > 0,
    nleft = nleft - 1;
    l = l+1;
    [p, position] = max(cube(:)); 
    [poiy,poix,poiz] = ind2sub(size(cube),position);
    locations.xyz(l,:) = [x(poix),y(poiy),z(poiz)];
    locations.p(l) = p;
    xa = poix - nx;
    if xa < 1, xa = 1; end;
    xe = poix + nx;
    if xe > mx, xe = mx; end;
    ya = poiy - ny;
    if ya < 1, ya = 1; end;
    ye = poiy + ny;
    if ye > my, ye = my; end;
    za = poiz - nz;
    if za < 1, za = 1; end;
    ze = poiz + nz;
    if ze > mz, ze = mz; end;
    for kx = xa:xe,
        rx2 = ((kx-poix)*dx)^2;
        for ky = ya:ye,
            ry2 = ((ky-poiy)*dy)^2;
            for kz = za:ze,
                rz2 = ((kz-poiz)*dz)^2;
                if rx2 + ry2 + rz2 <= d2, % this point is closer than d to the just selected location
                    cube(ky,kx,kz) = 0; % and is therefore masked
                end;
            end;
        end;
    end;
end;

function wr_locations(pid,all_locations,oname,all_tags,mirror_image)

if ~exist('mirror_image','var'),
    mirror_image = false;
end;

% generate header line
header=sprintf('HEADER    SITE NETWORK');
header=fillstr(header,50);
today=date;
today=[today(1:7) today(10:11)];
header=sprintf('%s%s   %s',header,today,pid);
% state supported format and originating program
format=sprintf('REMARK   4 %s COMPLIES WITH FORMAT V. 3.20, 01-DEC-08',pid);
origin=sprintf('REMARK   5 WRITTEN BY MMM (SITE NETWORK)');

if isempty(oname),
    fname = strcat('network_',pid);
    if mirror_image,
        fname = strcat(fname,'_inverted');
    end;
    fname = strcat(fname,'.pdb');
else
    if strfind(oname,'.pdb'),
        fname = oname;
        if mirror_image,
            fname = strcat(fname,'_inverted');
        end;
    else
        fname = [oname '.pdb'];
        if mirror_image,
            fname = [oname '_inverted.pdb'];
        end;
    end;
end;

fid=fopen(fname,'wt');
if fid==-1,
    message.error=1;
    message.text='File could not be written';
    return;
end;

fprintf(fid,'%s\n',header);
fprintf(fid,'TITLE     BAX SPIN LABEL LOCATIONS %s\n',pid);
fprintf(fid,'REMARK   4\n%s\n',format);
fprintf(fid,'REMARK   5\n%s\n',origin);

for k = 1: length(all_tags),
    tag = all_tags{k};
    rtags{k} = tag;
    if length(tag) > 3, rtags{k} = tag(1:3); end;
    while length(rtags{k})<3,
        rtags{k} = [rtags{k} 'Z'];
    end;
end;

for l = 1:length(all_locations{1}.p),
    atnum = 0; % initialize atom serial number
    if length(all_locations{1}.p)>1,
        fprintf(fid,'MODEL     %4i\n',l);
    end;
    cid='A';
    atag='NAL';
    tline='HETATM';
    hetflag=1;
    element='NA';
    if length(element)<2 && length(atag)<4,
        atag=[' ' atag];
    end;
    if length(atag)<4,
        atag=fillstr(atag,4);
    end;
    ltag = ' ';
    Bfactor=0.1;
    for rnum = 1:length(all_tags),
        locations = all_locations{rnum};
        rtag = rtags{rnum};
        rid=sprintf('%4i ',rnum);
        occupancy = locations.p(l);
        xyz=locations.xyz(l,:);
        if mirror_image,
            xyz = -xyz;
        end;
        atnum=atnum+1;
        fprintf(fid,'%s%5i %4s%s%s %s%4s   %8.3f%8.3f%8.3f',tline,atnum,atag,ltag,rtag,cid,rid,xyz);
        fprintf(fid,'%6.2f%6.2f          %2s\n',occupancy,Bfactor,element);
    end;
    if length(locations.p) > 1,
        fprintf(fid,'ENDMDL\n');
    end;
end;
fclose(fid);

function newstring=fillstr(string,newlength)
% pads a string with spaces

format=sprintf('%%s%%%is',newlength-length(string));
newstring=sprintf(format,string,' ');
% newstring = char(padarray(uint8(string)', newlength-length(string), 32,'post')');
