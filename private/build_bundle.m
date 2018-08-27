function build_bundle

% prototype version, according to the White paper 'Considerations on
% modeling helix arrangement in Bax'
% G. Jeschke, 2012

global general
global model

requested_cpu=5;

diagnostics=false; % set true for test purposes

very_long=30; % a very long upper distance bound for lower bound constraints (nm)
max_trials=1000000000;
max_time=5*3600;
max_trans=50; % maximum translation of a segment in any direction (Å)
tolerance=2.0;

% load test_oriented_6
% plot_ensemble(constraints,constraint_list,segment_list,geometry);
% keyboard

add_msg_board('### Helix bundle assembler ###');

stag=tag2id(model.structure_tags,model.current_structure);

if isfield(model,'sites'), % check, whether spin label rotamers were computed
    labels=label_information(model.sites);
    if ~isempty(labels), % create index table of known labels
        lindices=zeros(length(labels),4);
        for k=1:length(labels),
            cindices=labels(k).indices;
            if ~isempty(cindices),
                lindices(k,:)=cindices;
            end;
        end;
    else
        lindices=[];
    end;
else
    labels=[];
    lindices=[];
end;

% Get constraint file name and read constraint file, see Section 5.2 of
% White paper
[fname,pname]=uigetfile('*.dat','Load constraints from file');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Constraint loading cancelled by user');
    return
else
    reset_user_paths(pname);
    general.restraint_files=pname;
    constraints=rd_restraints(fullfile(pname,fname));
end;

if isfield(constraints,'initial_model') && ~isempty(constraints.initial_model),
    initial_model=get_initial_model(constraints.initial_model);
else
    initial_model.segments=[];
end;

if isfield(constraints,'max_time'),
    max_time=constraints.max_time*3600;
end;
if isfield(constraints,'tolerance'),
    tolerance=constraints.tolerance;
end;

% Get result file name
[fname,pname]=uiputfile('*.mat','Save results in file (.mat and .log)');
if isequal(fname,0) || isequal(pname,0)
    add_msg_board('Modeling cancelled by user');
    return
else
    filename=fullfile(pname,fname);
    [pathstr, name] = fileparts(filename);
    fname_mat=fullfile(pathstr,sprintf('%s.mat',name));
    fname_log=fullfile(pathstr,sprintf('%s.log',name));
end;

% prepare parallel computing
if requested_cpu>1 && exist('matlabpool'),
    % prepare parallel loop
    sz=matlabpool('size');
    if sz~=requested_cpu,
        if sz~=0,
            matlabpool('CLOSE');
        end;
        matlabpool('local',requested_cpu);
    end;
end;


% check consistency of PARALLEL and PERPENDICULAR statements

if isfield(constraints,'perpendicular') && isfield(constraints,'parallel'),
    for k=1:length(constraints.perpendicular),
        for kk=1:length(constraints.parallel),
            if constraints.perpendicular(k)==constraints.parallel(kk),
                add_msg_board(sprintf('ERROR: Segment %i is defined as both parallel and perpendicular to symmetry axis.',constraints.perpendicular(k)));
                return
            end;
        end;
    end;
end;

% determine the reference segment

if isfield(constraints,'bundle'),
    possible=constraints.bundle;
    if isfield(constraints,'parallel'), % block parallel helices from consideration as reference helices
        for k=1:length(constraints.parallel),
            for kk=1:length(possible),
                if possible(kk)==constraints.parallel(k),
                    possible(kk)=0;
                end;
            end;
        end;
    end;
    preferred=zeros(size(constraints.bundle));
    if isfield(constraints,'perpendicular'), % block parallel helices from consideration as reference helices
        for k=1:length(constraints.perpendicular),
            for kk=1:length(constraints.bundle),
                if constraints.bundle(kk)~=constraints.perpendicular(k),
                    preferred(kk)=1;
                end;
            end;
        end;
    end;
    preferred=preferred.*possible;
    poi=find(preferred);
    if ~isempty(poi),
        ref_segment=preferred(poi(1));
    else
        poi=find(possible);
        if ~isempty(poi),
            ref_segment=possible(poi(1));
        else
            add_msg_board('ERROR: No reference segment can be identified.');
            return
        end;
    end;
    add_msg_board(sprintf('The reference segment is segment %i.',ref_segment));
else
    add_msg_board('ERROR: No model defined (BUNDLE statement missing).');
    return
end;

% extract table of residue ranges for helical segments

if isfield(constraints,'helices') || length(constraints.helices)<1,
    segment_ranges=zeros(length(constraints.helices),2);
    for k=1:length(constraints.helices),
        address=constraints.helices(k).adr;
        hyphen=strfind(address,'-');
        segment_ranges(k,1)=str2double(address(1:hyphen-1));
        segment_ranges(k,2)=str2double(address(hyphen+1:end));
    end;
else
    add_msg_board('ERROR: No helical segments defined');
    return
end;

% Generate the constraint site list, see Section 5.3

sites_per_protomer=0;

if ~isfield(constraints,'ref_chain') || isempty(constraints.ref_chain),
    constraints.ref_chain=model.current_chain;
end;

constraint_tags=':';

if isfield(constraints,'DEER'),
    if isempty(labels) || isempty(lindices),
        add_msg_board('Warning: All DEER constraints are ignored as the current structure does not contain spin labels.');
    else
        for k=1:length(constraints.DEER),
            [isref,isnew,residue]=check_site_address(constraints.DEER(k).adr1,stag,constraints.ref_chain,constraint_tags,'.NO');
            if isref && isnew,
                cindices=resolve_address(sprintf('[%s]%s%s',stag,constraints.ref_chain,residue));
                if isempty(cindices),
                    add_msg_board(sprintf('Warning: DEER constraint site %s does not exist in template. Constraints with this site will be ignored.',constraints.DEER(k).adr1));
                else
                   found=false;
                   for l=1:length(labels),
                       diff=cindices-lindices(l,:);
                       if sum(abs(diff))==0 && strcmpi(labels(l).name,constraints.DEER(k).label),
                           found=l;
                       end;
                   end;
                   if ~found,
                       add_msg_board(sprintf('Warning: DEER constraint site %s is not labeled or label type is wrong. Constraints with this site will be ignored.',constraints.DEER(k).adr1));
                   else
                       sites_per_protomer=sites_per_protomer+1;
                       constraint_sites{sites_per_protomer}.chain=constraints.ref_chain;
                       constraint_sites{sites_per_protomer}.ref=true;
                       constraint_sites{sites_per_protomer}.res=residue;
                       constraint_sites{sites_per_protomer}.segment=get_segment(residue,segment_ranges);
                       constraint_sites{sites_per_protomer}.coor=labels(found).xyz;
                       constraint_tags=[constraint_tags residue '.NO:'];
                   end;
                end;
            end;
            [isref,isnew,residue]=check_site_address(constraints.DEER(k).adr2,stag,constraints.ref_chain,constraint_tags,'.NO');
            if isref && isnew,
                cindices=resolve_address(sprintf('[%s]%s%s',stag,constraints.ref_chain,residue));
                if isempty(cindices),
                    add_msg_board(sprintf('Warning: DEER constraint site %s does not exist in template. Constraints with this site will be ignored.',constraints.DEER(k).adr2));
                else
                   found=false;
                   for l=1:length(labels),
                       diff=cindices-lindices(l,:);
                       if sum(abs(diff))==0 && strcmpi(labels(l).name,constraints.DEER(k).label),
                           found=l;
                       end;
                   end;
                   if ~found,
                       add_msg_board(sprintf('Warning: DEER constraint site %s is not labeled or label type is wrong. Constraints with this site will be ignored.',constraints.DEER(k).adr2));
                   else
                       sites_per_protomer=sites_per_protomer+1;
                       constraint_sites{sites_per_protomer}.chain=constraints.ref_chain;
                       constraint_sites{sites_per_protomer}.ref=true;
                       constraint_sites{sites_per_protomer}.res=residue;
                       constraint_sites{sites_per_protomer}.segment=get_segment(residue,segment_ranges);
                       constraint_sites{sites_per_protomer}.coor=labels(found).xyz;
                       constraint_tags=[constraint_tags residue '.NO:'];
                   end;
                end;
            end;
        end;
    end;
end;

if isfield(constraints,'direct'),
    for k=1:length(constraints.direct),
        [isref,isnew,residue]=check_site_address(constraints.direct(k).adr1,stag,constraints.ref_chain,constraint_tags,'.CA');
        if isref && isnew,
            cindices=resolve_address(sprintf('[%s]%s%s.CA',stag,constraints.ref_chain,residue));
            if isempty(cindices) || length(cindices)~=5,
                add_msg_board(sprintf('Warning: Direct constraint site %s does not exist in template. Constraints with this site will be ignored.',constraints.direct(k).adr1));
            else
               [message,xyz]=get_atom(cindices,'coor');
               if message.error,
                    add_msg_board(sprintf('Warning: Failed to find Calpha coordinates fo direct constraint site %s (Errror: %s). Constraints with this site will be ignored.',constraints.direct(k).adr1,message.text));
               else
                   sites_per_protomer=sites_per_protomer+1;
                   constraint_sites{sites_per_protomer}.chain=constraints.ref_chain;
                   constraint_sites{sites_per_protomer}.ref=true;
                   constraint_sites{sites_per_protomer}.res=residue;
                   constraint_sites{sites_per_protomer}.segment=get_segment(residue,segment_ranges);
                   constraint_sites{sites_per_protomer}.coor=xyz;
                   constraint_tags=[constraint_tags residue '.CA:'];
               end;
            end;
        end;
        [isref,isnew,residue]=check_site_address(constraints.direct(k).adr2,stag,constraints.ref_chain,constraint_tags,'.CA');
        if isref && isnew,
            cindices=resolve_address(sprintf('[%s]%s%s.CA',stag,constraints.ref_chain,residue));
            if isempty(cindices) || length(cindices)~=5,
                add_msg_board(sprintf('Warning: Direct constraint site %s does not exist in template. Constraints with this site will be ignored.',constraints.direct(k).adr2));
            else
               [message,xyz]=get_atom(cindices,'coor');
               if message.error,
                    add_msg_board(sprintf('Warning: Failed to find Calpha coordinates fo direct constraint site %s (Errror: %s). Constraints with this site will be ignored.',constraints.direct(k).adr2,message.text));
               else
                   sites_per_protomer=sites_per_protomer+1;
                   constraint_sites{sites_per_protomer}.chain=constraints.ref_chain;
                   constraint_sites{sites_per_protomer}.ref=true;
                   constraint_sites{sites_per_protomer}.res=residue;
                   constraint_sites{sites_per_protomer}.segment=get_segment(residue,segment_ranges);
                   constraint_sites{sites_per_protomer}.coor=xyz;
                   constraint_tags=[constraint_tags residue '.CA:'];
               end;
            end;
        end;
    end;
end;

% Generate the segment list, see Section 5.4

sites_per_protomer=0; % the original site list includes sites that are not within the bundle
site_addresses=':';

k=0;
for sp=1:length(constraints.helices),
    poi=find(~(constraints.bundle-sp),1);
    if ~isempty(poi),
        k=k+1;
        segment_list{k}.number=sp;
        segment_list{k}.chain=constraints.ref_chain;
        if sp==ref_segment,
            segment_list{k}.ref=true;
        else
            segment_list{k}.ref=false;
        end;
        segment_list{k}.range=segment_ranges(sp,:);
        [msg,coor_list]=get_object(sprintf('[%s]%s%s.CA',stag,constraints.ref_chain,constraints.helices(sp).adr),'coor');
        coor=zeros(length(coor_list),3);
        for kk=1:length(coor_list),
            coor(kk,:)=coor_list{kk};
        end;
        [p0,vec]=rmsd_line_3D(coor);
        if diagnostics,
            figure(1); clf;
            hold on;
            [m,n]=size(coor);
            for kc=1:m-1,            
                plot3([coor(kc,1),coor(kc+1,1)],[coor(kc,2),coor(kc+1,2)],[coor(kc,3),coor(kc+1,3)],'k');
            end;
            plot3([p0(1)-vec(1)/2,p0(1)+vec(1)/2],[p0(2)-vec(2)/2,p0(2)+vec(2)/2],[p0(3)-vec(3)/2,p0(3)+vec(3)/2],'r');
        end;
        vecp=vec/norm(vec);
        th=acos(vecp(3));
        if norm(vecp(1:2))>1e-6,
            vecp=vecp(1:2)/norm(vecp(1:2));
        else
            vecp=[0,0];
        end;
        phi=atan2(vecp(2),vecp(1));
        transmat1=affine('translation',-p0);
        transmat2=affine('Euler',[-phi,-th,0]);
        transmat=transmat2*transmat1;
        [mm,nn]=size(coor);
        coor=[coor;p0;vec];
        xyz=[coor ones(mm+2,1)];
        xyz(mm+2,4)=0;
        xyz=transmat*xyz';
        xyz=xyz';
        coor=xyz(1:mm,1:3);
        segment_list{k}.Calpha=coor;

        segment_list{k}.origin=xyz(mm+1,1:3);
        segment_list{k}.axis=xyz(mm+2,1:3);
        p0=segment_list{k}.origin;
        vec=segment_list{k}.axis;
        if diagnostics,
            figure(2); clf;
            hold on;
            [m,n]=size(coor);
            for kc=1:m-1,            
                plot3([coor(kc,1),coor(kc+1,1)],[coor(kc,2),coor(kc+1,2)],[coor(kc,3),coor(kc+1,3)],'k');
            end;
            plot3([p0(1)-vec(1)/2,p0(1)+vec(1)/2],[p0(2)-vec(2)/2,p0(2)+vec(2)/2],[p0(3)-vec(3)/2,p0(3)+vec(3)/2],'r');
            keyboard;
        end;
        sites=zeros(1,100);
        spoi=0;
        for kk=1:length(constraint_sites),
            if constraint_sites{kk}.segment==sp,
                tag=id2tag(kk,constraint_tags);
                sites_per_protomer=sites_per_protomer+1;
                site_addresses=[site_addresses tag ':']; 
                spoi=spoi+1;
                sites(spoi)=kk;
            end;
        end;
        if spoi==0,
            sites=[];
            add_msg_board(sprintf('Warning: Segment %i is unconstrained.',sp));
        else
            sites=sites(1:spoi);
        end;
        segment_list{k}.sites=sites;
        coor=zeros(spoi,3);
        for kk=1:spoi,
            coor(kk,:)=constraint_sites{sites(kk)}.coor;
        end;
        [mm,nn]=size(coor);
        xyz=[coor ones(mm,1)];
        xyz=transmat*xyz';
        xyz=xyz';
        coor=xyz(1:mm,1:3);
        segment_list{k}.sitecoor=coor;
        if diagnostics,
            figure(2);
            for kk=1:spoi,
                plot3(coor(kk,1),coor(kk,2),coor(kk,3),'bo');
            end;
            keyboard
        end;
        segment_list{k}.vary=ones(1,6);
        segment_list{k}.fixed=zeros(1,6);
        if sp==ref_segment, % the reference segment has no gamma Euler angle
            segment_list{k}.vary(3)=0;
            segment_list{k}.fixed(3)=0;
        end;
        if isfield(constraints,'perpendicular'),
            for kk=1:length(constraints.perpendicular),
                if constraints.perpendicular(kk)==sp,
                    segment_list{k}.vary(2)=0;
                    segment_list{k}.fixed(2)=pi/2;
                    segment_list{k}.vary(6)=0;
                    segment_list{k}.fixed(6)=0;
                end;
            end;
        end;
        if isfield(constraints,'parallel'),
            for kk=1:length(constraints.parallel),
                if constraints.parallel(kk)==sp,
                    segment_list{k}.vary(2)=0;
                    segment_list{k}.fixed(2)=0;
                    segment_list{k}.vary(3)=0;
                    segment_list{k}.fixed(3)=0;
                end;
            end;
        end;
    end;
end;

% Generating, indexing and processing the constraint list, see Section 5.5

% Generate chain tag string 

chain_tags=[':' constraints.ref_chain ':'];
if isfield(constraints,'chain_copies'),
    for k=1:length(constraints.chain_copies),
        chain_tags=[chain_tags constraints.chain_copies{k} ':'];
    end;
end;

con_num=0;

% scan DEER constraints

if isfield(constraints,'DEER'),
    for k=1:length(constraints.DEER),
        [structure,chain,modtag,residue] = dissect_address(constraints.DEER(k).adr1);
        s1adr=[residue '.NO'];
        siteid1=tag2id(s1adr,constraint_tags);
        sitenum1=tag2id(s1adr,site_addresses);
        cid1=tag2id(chain,chain_tags);
        [structure,chain,modtag,residue] = dissect_address(constraints.DEER(k).adr2);
        s2adr=[residue '.NO'];
        siteid2=tag2id(s2adr,constraint_tags);
        sitenum2=tag2id(s2adr,site_addresses);
        cid2=tag2id(chain,chain_tags);
        if isempty(sitenum1) || isempty(cid1) || isempty(sitenum2) || isempty(cid2),
            add_msg_board(sprintf('Warning: DEER constraint between %s and %s could not be assigned in model.',...
                constraints.DEER(k).adr1,constraints.DEER(k).adr2));
        else
            con_num=con_num+1;
            constraint_list{con_num}.type='NO';
            constraint_list{con_num}.adr1=constraints.DEER(k).adr1;
            constraint_list{con_num}.adr2=constraints.DEER(k).adr2;
            constraint_list{con_num}.sites=[siteid1 siteid2];
            coornum1=(cid1-1)*sites_per_protomer+sitenum1;
            coornum2=(cid2-1)*sites_per_protomer+sitenum2;
            constraint_list{con_num}.coor_pointer=[coornum1 coornum2];
            arg1=constraints.DEER(k).r;
            arg2=constraints.DEER(k).sigr;
            type=0;
            if arg1>0 && arg2>0, % mean distance/std. deviation
                constraint_list{con_num}.mode=0;
                constraint_list{con_num}.r=arg1;
                constraint_list{con_num}.sigr=arg2;
                constraint_list{con_num}.bounds=[arg1-tolerance*arg2,arg1+tolerance*arg2];                
            end;
            if arg1<0 && arg2>=0, % lower bound
                constraint_list{con_num}.mode=1;
                constraint_list{con_num}.bounds=[-arg1,very_long];                
            end;
            if arg1>=0 && arg2<0, % upper bound
                constraint_list{con_num}.mode=2;
                constraint_list{con_num}.bounds=[0,-arg2];                
            end;
            if arg1<0 && arg2<0, % lower bound/upper bound pair
                constraint_list{con_num}.mode=3;
                constraint_list{con_num}.bounds=[-arg1,-arg2];                
            end;                
        end;
    end;
end;

% scan direct constraints

if isfield(constraints,'direct'),
    for k=1:length(constraints.direct),
        [structure,chain,modtag,residue] = dissect_address(constraints.direct(k).adr1);
        s1adr=[residue '.CA'];
        siteid1=tag2id(s1adr,constraint_tags);
        sitenum1=tag2id(s1adr,site_addresses);
        cid1=tag2id(chain,chain_tags);
        [structure,chain,modtag,residue] = dissect_address(constraints.direct(k).adr2);
        s2adr=[residue '.CA'];
        siteid2=tag2id(s2adr,constraint_tags);
        sitenum2=tag2id(s2adr,site_addresses);
        cid2=tag2id(chain,chain_tags);
        if isempty(sitenum1) || isempty(cid1) || isempty(sitenum2) || isempty(cid2),
            add_msg_board(sprintf('Warning: Direct Calpha constraint between %s and %s could not be assigned in model.',...
                constraints.direct(k).adr1,constraints.direct(k).adr2));
        else
            con_num=con_num+1;
            constraint_list{con_num}.type='CA';
            constraint_list{con_num}.adr1=constraints.direct(k).adr1;
            constraint_list{con_num}.adr2=constraints.direct(k).adr2;
            constraint_list{con_num}.sites=[siteid1 siteid2];
            coornum1=(cid1-1)*sites_per_protomer+sitenum1;
            coornum2=(cid2-1)*sites_per_protomer+sitenum2;
            constraint_list{con_num}.coor_pointer=[coornum1 coornum2];
            arg1=constraints.direct(k).r;
            arg2=constraints.direct(k).sigr;
            type=0;
            if arg1>0 && arg2>0, % mean distance/std. deviation
                constraint_list{con_num}.mode=0;
                constraint_list{con_num}.r=arg1;
                constraint_list{con_num}.sigr=arg2;
                constraint_list{con_num}.bounds=[arg1-arg2,arg1+arg2];                
            end;
            if arg1<0 && arg2>=0, % lower bound
                constraint_list{con_num}.mode=1;
                constraint_list{con_num}.bounds=[-arg1,very_long];                
            end;
            if arg1>=0 && arg2<0, % upper bound
                constraint_list{con_num}.mode=2;
                constraint_list{con_num}.bounds=[0,-arg2];                
            end;
            if arg1<0 && arg2<0, % lower bound/upper bound pair
                constraint_list{con_num}.mode=3;
                constraint_list{con_num}.bounds=[-arg1,-arg2];                
            end;                
        end;
    end;
end;

% Monte Carlo trial loop, Section 5.6

% prepare statistics

if ~isempty(initial_model.segments),
    params=initial_model.info{1}.parameters;
    [mi,~]=size(params);
    models_used=zeros(1,mi);
else
    mi=1;
end;

fulfilled=zeros(1,length(constraint_list));
num_fulfilled=zeros(1,length(constraint_list));
helix_clashes=0;

% initialize data arrays
protomers=length(constraints.chain_copies)+1;
model_sites=protomers*sites_per_protomer;
% initialize site coordinate array
full_coor=zeros(model_sites,3);
% initialize parameter vector
parameter_vector=zeros(1,6*length(segment_list));

mck=0;
poi=1;
runtime=0;
max_fulfilled=0;
tic;

profile on

identity_affine=eye(4);

while mck<max_trials && poi<constraints.ensemble+1 && runtime<max_time, 
    mck=mck+1;
    runtime=toc;
    sel=round(0.5+mi*rand);
    if sel<1, sel=1; end;
    if sel>mi, sel=mi; end;
    % Generating a geometry (Section 5.6.1)
    for k=1:length(segment_list),
        bas=6*(k-1);
        ipoi=find(initial_model.segments==segment_list{k}.number);
        if isempty(ipoi),
            if segment_list{k}.vary(1), % alpha_i
                parameter_vector(bas+1)=2*pi*rand;
            else
                parameter_vector(bas+1)=segment_list{k}.fixed(1);
            end;
            if segment_list{k}.vary(2), % beta_i
                parameter_vector(bas+2)=acos(rand);
            else
                parameter_vector(bas+2)=segment_list{k}.fixed(2);
            end;
            if segment_list{k}.vary(3), % gamma_i
                parameter_vector(bas+3)=2*pi*rand;
            else
                parameter_vector(bas+3)=segment_list{k}.fixed(3);
            end;
            for kk=4:6, % translations x_i, y_i, z_i
                if segment_list{k}.vary(kk),
                    parameter_vector(bas+kk)=2*(rand-0.5)*max_trans;
                else
                    parameter_vector(bas+kk)=segment_list{k}.fixed(kk);
                end;
            end;
        else
            params=initial_model.info{ipoi}.parameters;
            parameter_vector(bas+1:bas+6)=params(sel,:);
        end;
    end;
    geometry(poi).transformation=parameter_vector;
    for k=1:length(segment_list),
        bas=6*(k-1);
        euler_list{k}=parameter_vector(bas+1:bas+3);
        translist{k}=parameter_vector(bas+4:bas+6);
    end;
    % transform coordinates for the reference chain
    parfor k=1:length(segment_list),
        euler=euler_list{k};
        transmat1=Euler_affine_fast(euler,identity_affine);
        translation=translist{k};
        transmat2=identity_affine;
        transmat2(1:3,4)=translation';
        coor=segment_list{k}.sitecoor;
        [mm,nn]=size(coor);
        transmat=transmat2*transmat1;
        xyz=[coor ones(mm,1)];
        xyz=transmat*xyz';
        xyz=xyz';
        xyz_list{k}=xyz(:,1:3);
    end;
    coorpoi=0;
    for k=1:length(segment_list),
        coor=xyz_list{k};
        [mm,~]=size(coor);
        full_coor(coorpoi+1:coorpoi+mm,:)=coor;
        coorpoi=coorpoi+mm;
    end;
    protomer_coor=full_coor(1:coorpoi,:);
    % make symmetry-related protomers
    symn=constraints.Cn;
    parfor j=2:constraints.Cn,
        phi=2*pi*(j-1)/symn;
        transmat=rotz_affine_fast(phi,identity_affine);
        [mm,~]=size(protomer_coor);
        xyz=[protomer_coor ones(mm,1)];
        xyz=transmat*xyz';
        xyz=xyz';
        xyz_list{j-1}=xyz(:,1:3);
    end;
    for j=2:constraints.Cn,
        full_coor((j-1)*coorpoi+1:j*coorpoi,:)=xyz_list{j-1};
    end;
    geometry(poi).sitecoor=full_coor;
    % Testing the model geometry
    rmsd=0;
    chisquare=0;
    failed=false;
    nsoft=0;
    success=0;
    failures=zeros(1,length(constraint_list));
    parfor k=1:length(constraint_list),
        cpoi=constraint_list{k}.coor_pointer;
        xyz1=full_coor(cpoi(1),:);
        xyz2=full_coor(cpoi(2),:);
        r=norm(xyz1-xyz2)/10; % conversion Angstrome -> nanometers
        if r<constraint_list{k}.bounds(1) || r>constraint_list{k}.bounds(2),
            failures(k)=1;
        else
            fulfilled(k)=fulfilled(k)+1;
            success=success+1;
        end;
        if constraint_list{k}.mode==0,
            d=r-constraint_list{k}.r;
            nsoft=nsoft+1;
            rmsd=rmsd+d^2;
            chisquare=chisquare+(d/constraint_list{k}.sigr)^2;
        end;
    end;
    failed=sum(failures);
    if success>0,
        num_fulfilled(success)=num_fulfilled(success)+1;
    end;
    if success>max_fulfilled,
        max_fulfilled=success;
        add_msg_board(sprintf('Best model so far fulfilled %i out of %i constraints.',success,length(constraint_list)));
    end;
    if failed,
        continue;
    end;
    if nsoft>0,
        geometry(poi).rmsd=sqrt(rmsd/nsoft);
        geometry(poi).chisquare=chisquare/nsoft;
    else
        geometry(poi).rmsd=0;
        geometry(poi).chisquare=0;
    end;
    % Test for helix-helix axis clashes
    % add_msg_board(sprintf('Now testing for helix-helix clashes in model %i',poi));
    failed=false;
    helices=constraints.Cn*length(segment_list);
    helix_ends=zeros(2*helices,3);
    for k=1:length(segment_list),
        bas=6*(k-1);
        coor=zeros(2,3);
        coor(1,:)=segment_list{k}.origin-segment_list{k}.axis/2;
        coor(2,:)=segment_list{k}.origin+segment_list{k}.axis/2;
        euler=parameter_vector(bas+1:bas+3);
        transmat1=Euler_affine_fast(euler,identity_affine);
        translation=parameter_vector(bas+4:bas+6);
        transmat2=identity_affine;
        transmat2(1:3,4)=translation';
        [mm,nn]=size(coor);
        transmat=transmat2*transmat1;
        xyz=[coor ones(mm,1)];
        xyz=transmat*xyz';
        xyz=xyz';
        bas2=2*(k-1);
        helix_ends(bas2+1:bas2+2,:)=xyz(:,1:3);
    end;
    protomer_hends=helix_ends(1:2*length(segment_list),:);
    % make symmetry-related helix axes
    for j=2:constraints.Cn,
        phi=2*pi*(j-1)/constraints.Cn;
        transmat=rotz_affine_fast(phi,identity_affine);
        [mm,nn]=size(protomer_hends);
        xyz=[protomer_hends ones(mm,1)];
        xyz=transmat*xyz';
        xyz=xyz';
        helix_ends((j-1)*2*length(segment_list)+1:2*j*length(segment_list),:)=xyz(:,1:3);
    end;
    geometry(poi).helix_ends=helix_ends;
    % test for clashes between any helix pairs
    for k=1:constraints.Cn*length(segment_list)-1,
        bas1=2*(k-1);
        for kk=k+1:constraints.Cn*length(segment_list),
            bas2=2*(kk-1);
            dist=vec_dist(helix_ends(bas1+1,:),helix_ends(bas1+2,:),helix_ends(bas2+1,:),helix_ends(bas2+2,:));
            % fprintf(1,'(%i,%i): %5.1f Å\n',k,kk,dist);
            if dist/10<constraints.helix_mindist,
                if ~failed,
                    helix_clashes=helix_clashes+1;
                end;
                failed=true;
            end;
        end;
    end;
    if ~failed,
        poi=poi+1;
        add_msg_board(sprintf('Geometry model %i was found.',poi-1)); 
        if ~isempty(initial_model.segments),
            models_used(sel)=models_used(sel)+1;
        end;
    end;
end;
poi=poi-1;
if poi>=1,
    geometry=geometry(1:poi);
else
    geometry=[];
end;
add_msg_board(sprintf('Ensemble with %i structures generated.',poi));
if mck>=max_trials,
    add_msg_board(sprintf('Maximum number of Monte Carlo trials (%i) reached.',max_trials));
end;
if runtime>=max_time,
    add_msg_board(sprintf('Maximum runtime (%i s) reached.',max_time));
end;

profile viewer

% save results in binary form
save(fname_mat,'constraints','geometry','constraint_list','segment_list','constraint_sites');

if ~isempty(geometry),
    [full_rms,max_trans]=plot_ensemble(constraints,constraint_list,segment_list,geometry);
end;

% save logfile
fid=fopen(fname_log,'wt');
fprintf(fid,'+++ MMM Helix Bundle Assembler logfile +++\n\n');
fprintf(fid,'Run limited to %i Monte Carlo trials with maximum runtime of %4.1f s.\n',max_trials,max_time);
fprintf(fid,'Run finished with %i Monte Carlo trials and runtime of %4.1f s.\n\n',mck,runtime);

fprintf(fid,'Reference chain: %s\n\n',constraints.ref_chain);
fprintf(fid,'There are %i active segments:\n',length(segment_list));
free_param=0;
for k=1:length(segment_list),
    rgb=color_grade(k,length(segment_list));
    fprintf(fid,'%i: residues %i-%i [rgb(%5.3f,%5.3f,%5.3f)]\n   varied: ',segment_list{k}.number,segment_list{k}.range(1),segment_list{k}.range(2),rgb(1),rgb(2),rgb(3));
    if segment_list{k}.vary(1), fprintf(fid,'alpha '); end;
    if segment_list{k}.vary(2), fprintf(fid,'beta '); end;
    if segment_list{k}.vary(3), fprintf(fid,'gamma '); end;
    if segment_list{k}.vary(4), fprintf(fid,'dx '); end;
    if segment_list{k}.vary(5), fprintf(fid,'dy '); end;
    if segment_list{k}.vary(6), fprintf(fid,'dz '); end;
    fprintf(fid,'\n');
    free_param=free_param+sum(segment_list{k}.vary);
    if sum(segment_list{k}.vary)<6,
        fprintf(fid,'   fixed: ');
        if ~segment_list{k}.vary(1), fprintf(fid,'alpha= %4.1f° ',180*segment_list{k}.fixed(1)/pi); end;
        if ~segment_list{k}.vary(2), fprintf(fid,'beta= %4.1f° ',180*segment_list{k}.fixed(2)/pi); end;
        if ~segment_list{k}.vary(3), fprintf(fid,'gamma= %4.1f° ',180*segment_list{k}.fixed(3)/pi); end;
        if ~segment_list{k}.vary(4), fprintf(fid,'dx= %3.1f nm ',segment_list{k}.fixed(4)); end;
        if ~segment_list{k}.vary(5), fprintf(fid,'dy= %3.1f nm ',segment_list{k}.fixed(5)); end;
        if ~segment_list{k}.vary(6), fprintf(fid,'dz= %3.1f nm ',segment_list{k}.fixed(6)); end;
    end;
    fprintf(fid,'\n');
end;
fprintf(fid,'\n');

fprintf(fid,'There are %i free parameters at DEER tolerance of %3.1f*sigr.\n\n',free_param,tolerance);
fprintf(fid,'There are %i active constraints:\n\n',length(constraint_list));
for k=1:length(constraint_list),
    fprintf(fid,'%i.%s: [%s,%s] with bounds (%4.2f,%4.2f) nm\n   mode: ',k,...
        constraint_list{k}.type,constraint_list{k}.adr1,constraint_list{k}.adr2,...
        constraint_list{k}.bounds(1),constraint_list{k}.bounds(2));
    switch constraint_list{k}.mode
        case 0
            fprintf(fid,'soft (r= %4.2f, sigr= %4.2f) nm',constraint_list{k}.r,constraint_list{k}.sigr);
        case 1
            fprintf(fid,'Lower bound');
        case 2
            fprintf(fid,'Upper bound');
        case 3
            fprintf(fid,'Upper and lower bound');
    end;
    fprintf(fid,'\n');
    fprintf(fid,'   Fulfilled in %5.2f%% of all trials.\n',100*fulfilled(k)/mck);
end;
fprintf(fid,'\n');

if ~isempty(initial_model.segments),
    fprintf(fid,'Assembly based on initial model for segments ');
    for k=1:length(initial_model.segments)-1,
        fprintf(fid,'%i, ',initial_model.segments(k));
    end;
    fprintf(fid,'%i.\n',initial_model.segments(end));
    really_used=sum(models_used>0);
    fprintf(fid,'Out of %i initial models %i were used.\n',length(models_used),really_used);
    if really_used<=10,
        fprintf(fid,'\n   Usage statistics:\n');
        total_used=sum(models_used);
        for k=1:length(models_used),
            if models_used(k)>0,
                fprintf(fid,'Initial model %i used for %5.1f%% of final models.\n',k,100*models_used(k)/total_used);
            end;
        end;
    end;
    fprintf(fid,'\n');
end;

fprintf(fid,'Success statistics:\n');
for k=1:length(num_fulfilled),
    fprintf(fid,'%2i constraints are fulfilled in %8.5f%% of all trials (%i).\n',k,100*num_fulfilled(k)/mck,num_fulfilled(k));
end;
fprintf(fid,'\n');

if num_fulfilled(end)>0,
    fprintf(fid,'Of geometries fulfilling all constraints, %5.2f%% have helix-helix clashes.\n',100*helix_clashes/num_fulfilled(end));
end;
if ~isempty(geometry),
    fprintf(fid,'Ensemble with %i models was created.\n',length(geometry));
    fprintf(fid,'Mean r.m.s.d. with respect to first model is %5.2f Å\n',full_rms);
    fprintf(fid,'Maximum segment translation in any direction is %5.2f Å\n',max_trans);
end;

fclose(fid);


function [isref,isnew,residue]=check_site_address(address,stag,ref_chain,constraint_tags,type)

isnew=false;

[structure,chain,modtag,residue] = dissect_address(address);
if ~length(structure)>2 && ~strcmpi(structure(2:end-1),stag),
    add_msg_board('Warning: Constraint site address %s not in current structure.',address);
end;
if ~isempty(modtag),
    add_msg_board('Warning: Chain model number ignored in constraint site address %s.',address);
end;

isref=strcmpi(chain,ref_chain);

if isref,
    found=tag2id([residue type],constraint_tags);
    if isempty(found),
        isnew=true;
    else
        isnew=false;
    end;
end;

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
        end;
    end;
end;

function segment=get_segment(residue,segment_ranges)

[m,n]=size(segment_ranges);
resnum=str2double(residue);
mindiff=10000;
for k=1:m,
    if resnum>= segment_ranges(k,1) && resnum<=segment_ranges(k,2),
        segment=k;
        return;
    end;
    diff=min(abs([segment_ranges(k,1)-resnum,resnum-segment_ranges(k,2)])); % minimum distance from this segment
    if diff<mindiff,
        mindiff=diff;
        segment=k;
    end;
end;

function [full_rms,max_trans]=plot_ensemble(constraints,constraint_list,segment_list,geometry)
% plots an ensemble of structures

diagnostics=false;

figure(13); clf;
if isfield(constraints,'Cn'),
    n=constraints.Cn;
    plot3([0,0],[0,0],[-30,30],'Linewidth',3,'Color','k'); % Cn symmetry axis
else
    n=1;
end;
hold on
full_rms=0;
max_trans=0;
for k=1:length(geometry),
    trafo=geometry(k).transformation;
    for kk=1:length(trafo),
        if mod(kk,6)>3,
            if abs(trafo(kk))>max_trans, max_trans=abs(trafo(kk)); end;
        end;
    end;
    if k==1,
        template=plot_geometry(segment_list,constraint_list,geometry(k).transformation,n);
        new_coor=template;
    else
        [new_coor,rms]=plot_geometry(segment_list,constraint_list,geometry(k).transformation,n,template);
        fprintf(1,'Geometry %i superimposes on geometry 1 with r.m.s.d. of %4.2f Å.\n',k,rms);  
        full_rms=full_rms+rms^2;
    end;
    reslist=[];
    for k=1:length(segment_list),
        reslist=[reslist segment_list{k}.range(1):segment_list{k}.range(2)];
    end;
    hcoor{1}=new_coor;
    for j=2:n,
        phi=2*pi*(j-1)/n;
        transmat=affine('rotz',phi);
        [mm,nn]=size(new_coor);
        xyz=[new_coor ones(mm,1)];
        xyz=transmat*xyz';
        xyz=xyz';
        hcoor{j}=xyz(:,1:3);
    end;

    for kk=1:length(constraint_list),
        [~,chain,~,residue]=dissect_address(constraint_list{kk}.adr1);
        if strcmpi(chain,'(A)'),
            cid1=1;
        else
            cid1=2;
        end;
        [~,poi1]=min(abs(reslist-str2double(residue)));
        [~,chain,~,residue]=dissect_address(constraint_list{kk}.adr2);
        if strcmpi(chain,'(A)'),
            cid2=1;
        else
            cid2=2;
        end;
        coor1=hcoor{cid1};
        coor2=hcoor{cid2};
        [~,poi2]=min(abs(reslist-str2double(residue)));
        if diagnostics,
            fprintf(1,'Constraint: %s.(%s,%s):',constraint_list{kk}.type,constraint_list{kk}.adr1,constraint_list{kk}.adr2);
            fprintf(1,' bounds: (%5.2f,%5.2f) nm\n',constraint_list{kk}.bounds(1),constraint_list{kk}.bounds(2));
            fprintf(1,'Geometry  : CA.(%i,%i):',reslist(poi1),reslist(poi2));
            fprintf(1,' found : %5.2f nm\n',norm(coor1(poi1,:)-coor2(poi2,:))/10);
        end;
    end;
end;
full_rms=sqrt(full_rms/(length(geometry)-1));
fprintf(1,'Mean r.m.s.d. with respect to structure 1: %5.2f Å.\n',full_rms);
fprintf(1,'Maximum translation: %5.2f Å.\n',max_trans);

function [new_coor,rms]=plot_geometry(segment_list,constraint_list,parameter_vector,n,template)
% plots a helix bundle geometry for a given parameter vector of 
% transformation parameters and C_n symmetry with multiplicity n with respect to the z axis 

segments=zeros(length(segment_list),2);
poi=0;
identity_affine=eye(4);
for k=1:length(segment_list),
    bas=6*(k-1);
    euler=parameter_vector(bas+1:bas+3);
    transmat1=Euler_affine_fast(euler,identity_affine);
    translation=parameter_vector(bas+4:bas+6);
    transmat2=translation_affine_fast(translation,identity_affine);
    coor=segment_list{k}.Calpha;
    [mm,nn]=size(coor);
    segments(k,:)=[poi+1,poi+mm];
    poi=poi+mm;
    transmat=transmat2*transmat1;
    xyz=[coor ones(mm,1)];
    xyz=transmat*xyz';
    xyz=xyz';
    Calpha=xyz(:,1:3);
    if k==1,
        new_coor=Calpha;
    else
        new_coor=[new_coor;Calpha];
    end;
end;
if nargin>4,
    [rms,new_coor]=rmsd_superimpose(template,new_coor);
end;
for k=1:length(segment_list),
    Calpha=new_coor(segments(k,1):segments(k,2),:);
    plot_segment(Calpha,k,length(segment_list));
    if n>1,
        % make and plot symmetry-related segments
        for j=2:n,
            phi=2*pi*(j-1)/n;
            transmat=affine('rotz',phi);
            [mm,nn]=size(Calpha);
            xyz=[Calpha ones(mm,1)];
            xyz=transmat*xyz';
            xyz=xyz';
            Calpha_j=xyz(:,1:3);
            plot_segment(Calpha_j,k,length(segment_list));
        end;
    end;
end;



function plot_segment(Cacoor,curr,num)
% plots a helical segment with color coding its number
% Cacoor Calpha coordinates
% curr   current segment number
% num    total number of segments

rgb=color_grade(curr,num);
[m,n]=size(Cacoor);
plot3([Cacoor(1,1),Cacoor(1,1)],[Cacoor(1,2),Cacoor(1,2)],[Cacoor(1,3),Cacoor(1,3)],'k>');
plot3([Cacoor(end,1),Cacoor(end,1)],[Cacoor(end,2),Cacoor(end,2)],[Cacoor(end,3),Cacoor(end,3)],'ko');
for kc=1:m-1,            
    plot3([Cacoor(kc,1),Cacoor(kc+1,1)],[Cacoor(kc,2),Cacoor(kc+1,2)],[Cacoor(kc,3),Cacoor(kc+1,3)],'LineWidth',2,'Color',rgb);
end;

function initial_model=get_initial_model(fname)

initial_model.segments=[];
try
    load(fname);
catch
    return
end;

if ~exist('segment_list','var') || ~exist('geometry','var') || isempty(geometry),
    return
end;

initial_model.segments=zeros(1,length(segment_list));
for k=1:length(segment_list),
    initial_model.segments(k)=segment_list{k}.number;
    initial_model.info{k}.parameters=zeros(length(geometry),6);
end;

initial_model.chisquare=zeros(1,length(geometry));
initial_model.rmsd=zeros(1,length(geometry));

for k=1:length(geometry),
    if isfield(geometry,'chisquare'),
        initial_model.chisquare(k)=geometry(k).chisquare;
    end;
    if isfield(geometry,'rmsd'),
        initial_model.rmsd(k)=geometry(k).rmsd;
    end;
    trafo=geometry(k).transformation;
    for kk=1:length(segment_list),
        bas=6*(kk-1);
        initial_model.info{kk}.parameters(k,:)=trafo(bas+1:bas+6);
    end;
end;


