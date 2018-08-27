function refine_ensemble(snum,restraints,idCode)
% Refines an ensemble of structures using a list of user supplied 
% restraints by calling Modeller
% each structure in the ensemble is refined individually
%
% snum          (optional) index of the structure to be refined, defaults 
%               to current structure
% restraints    restraint lists in the format supplied by rd_restraints, if
%               not present, user is prompted to load restraints
% idCode        identifier code for target structure (cut or filled to four 
%               characters)

global model
global help_files
global third_party
global general

my_path=pwd;
cd(general.restraint_files);

if nargin<1,
    snum=model.current_structure;
end;

if nargin<2,
    [fname,pname]=uigetfile('*.dat','Load restraints from file');
    if isequal(fname,0) || isequal(pname,0)
        add_msg_board('ERROR: Structure refinement cancelled.');
        add_msg_board('Restraint loading cancelled by user');
        return
    else
        reset_user_paths(pname);
        general.restraint_files=pname;
        restraints=rd_restraints(fullfile(pname,fname));
    end;
end;

if nargin<3,
    answer = inputdlg('Please provide identifier code for target structure (four characters)','Ensemble refinement');
    idCode=strtrim(char(answer));
    if length(idCode)~=4,
        add_msg_board('Warning: Identifier code must have four characters');
        add_msg_board('Automatic adjustment.');
    end;
end;

if length(idCode)>4,
    idCode=idCode(1:4);
end;
while length(idCode)<4,
    idCode=[idCode 'x'];
end;

add_msg_board(sprintf('Refined ensemble has identifier %s.',idCode));

if isempty(restraints),
    add_msg_board('ERROR: Ensemble refinement aborted.');
    add_msg_board('Empty restraint list.');
    return
end;

hfig=gcf;
set(hfig,'Pointer','watch');

entry=strcat(help_files,'third_party.html#Modeller');

dospath=which([third_party.modeller_version  '.exe']);
if isempty(dospath),
    message.error=2;
    message.text='Modeller software not found on Matlab path.';
    add_msg_board('This feature requires Modeller from the Sali lab');
    add_msg_board('ERROR: Modeller could not be found on the Matlab path.');
    add_msg_board('Please check whether Modeller is installed and the path set.');
    add_msg_board('(see also help browser)');
    set(hfig,'Pointer','arrow');
    webcall(entry,'-helpbrowser');
    return
end;
[modpath, modcmd] = fileparts(dospath);

clean_targ_seq=model.structures{snum}(1).sequence;

models=length(model.structures{snum}(1).xyz);
for mnum=1:models,
    add_msg_board(sprintf('Refining model %i out of %i models.',mnum,models));
    [DEER,cancelled]=process_DEER_restraints(snum,mnum,restraints);
    if cancelled,
        add_msg_board('ERROR: Ensemble refinement aborted.');
        add_msg_board('Processing of DEER constraints failed.');
        set(hfig,'Pointer','arrow');
        return
    end;
    [direct,cancelled]=process_direct_restraints(snum,mnum,restraints);
    if cancelled,
        add_msg_board('ERROR: Ensemble refinement aborted.');
        add_msg_board('Processing of direct constraints cancelled.');
        set(hfig,'Pointer','arrow');
        return
    end;
    if isempty(DEER) && isempty(direct),
       add_msg_board('ERROR: Ensemble refinement aborted.');
        add_msg_board('No experimental restraints.');
        set(hfig,'Pointer','arrow');
        return
    end;
    filename='temp.pdb';

    fname=fullfile(general.tmp_files, filename);
    msg=sprintf('Template structure saved as PDB file: %s',fname);
    add_msg_board(msg);
    [message,info]=wr_pdb_paradigm(fname,'temp',mnum);
    if message.error,
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    end;
    basname=sprintf('temp_to_%s',idCode); 
    [alignfile,message]=seq2pir(basname,'temp',idCode,info.seq,info.first,info.first_chain,info.last,info.last_chain,'',info.first);

    if message.error,
        add_message_board('ERROR: Ensemble refinement.');
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    end;

    add_msg_board(sprintf('Alignment file saved as: %s',alignfile));
    [algpath, algfile, algext] = fileparts(alignfile);
    alignfile=strcat(algfile,algext);

    runfilename=sprintf('model_temp_to_%s.py',idCode);
    runfile=fullfile(general.tmp_files, runfilename);
    message=mk_modeller_input(runfile,alignfile,idCode,snum,mnum,restraints,DEER,direct,info.first);

    if message.error,
        add_msg_board(message.text);
        add_message_board('ERROR: Modelling aborted.');
        set(hwin,'Pointer','arrow');
        return
    end;

    msg='Running Modeller. Please be patient...';

    [batcmd,message]=mk_modeller_bat(modpath,modcmd,runfile);
    if message.error,
        add_msg_board('ERROR: Modeller could not be initialized');
        add_msg_board(message.text);
        set(hfig,'Pointer','arrow');
        return
    end;

    add_msg_board(msg);
    drawnow;
    my_dir=pwd;
    cd(modpath);
    [s, w] = dos(batcmd);

    if s~=0,
        rem=w;
        while ~isempty(rem),
            [token,rem]=strtok(rem,char(10));
            if ~isempty(token),
                add_msg_board(token);
            end;
        end;
        message.error=2;
        message.text='Modeller error.';
        add_msg_board('ERROR: Modeller did not run successfully.');
        set(hfig,'Pointer','arrow');
        cd(my_dir);
        return
    end;
    cd(my_dir);

end;
[message,snum,stag,models]=rd_pdb_Modeller_ensemble(idCode,1,models,general.tmp_files);
if message.error,
    add_msg_board(message.text);
else
    add_msg_board(sprintf('Ensemble read into structure %i with tag [%s].',snum,stag));
    add_msg_board(sprintf('This ensemble has %i chain models',models));
end;
if isfield(restraints,'chains'),
    chain_tags=restraints.chains;
else
    chain_tags='';
end;
correct_Modeller_ensemble(snum,info.first,chain_tags);

make_report(snum,DEER);

id=[stag 'r'];
exists=tag2id(id,model.structure_tags);
while ~isempty(exists),
    id=[id 'a'];
    exists=tag2id(id,model.structure_tags);
end;
[snum1,message]=copy_structure(snum,id);
add_msg_board('Now removing spin labels...');
[repnum,msg]=replace(snum1,':R1A:IA1:',false,false);
add_msg_board(sprintf('%i spin labels were removed.',repnum));
add_msg_board('Now repacking side chains...');
clean_targ_seq=clean_targ_seq(info.first:end);
for modnum=1:models,
    message=repack(snum1,modnum,0,clean_targ_seq);
end;
msg='Refined ensemble was succesfully imported.';
add_msg_board(msg);
cd(my_path);

set(hfig,'Pointer','arrow');


function [DEER,cancelled]=process_DEER_restraints(snum,mnum,restraints)

global model

cancelled=false;

if ~isfield(restraints,'DEER'),
    DEER=[];
    return;
end;

restraints=label_all_sites(snum,mnum,restraints);

labels=label_information(model.sites);

for k=1:length(restraints.DEER),
    adr1=restraints.DEER(k).adr1;
    ind1=resolve_address(adr1);
    adr2=restraints.DEER(k).adr2;
    ind2=resolve_address(adr2);
    DEER(k).r=restraints.DEER(k).r;
    DEER(k).sigr=restraints.DEER(k).sigr;
    DEER(k).indices=[ind1;ind2];
    DEER(k).adr1=adr1;
    DEER(k).adr2=adr2;
    DEER(k).label=restraints.DEER(k).label;
    DEER(k).T=restraints.DEER(k).T;
    f1=false;
    f2=false;
    for l=1:length(labels),
        diff1=ind1-labels(l).indices;
        if sum(abs(diff1))==0,
            f1=true;
            DEER(k).xyz1=labels(l).xyz;
            DEER(k).rmsd1=labels(l).rmsd;
        end;
        diff2=ind2-labels(l).indices;
        if sum(abs(diff2))==0,
            f2=true;
            DEER(k).xyz2=labels(l).xyz;
            DEER(k).rmsd2=labels(l).rmsd;
        end;
    end;
    if ~f1 || ~f2,
        add_msg_board('ERROR: Automatic rotamer computation error.');
        add_msg_board('Please mail gunnar.jeschke@phys.chem.ethz.ch');
        cancelled=true;
        DEER=[];
        return;
    end;
end;

function restraints=label_all_sites(snum,mnum,restraints)
% Given a restraint array, all spin labels are computed that are required
% for specifying DEER restraints
% restraints for which the labeling sites do not exist are removed
% labels are attached
%
% restraints    array of existing restraints
%

global model
global hMain

if isempty(restraints.DEER),
    return
else
    DEER=restraints.DEER;
    poi=0;
end;

if isfield(model,'sites'),
    labels=label_information(model.sites);
else
    labels=[];
end;

if ~isempty(DEER),
    add_msg_board(sprintf('Generating labels for %i DEER restraints',length(DEER)));
end;
if ~isempty(labels),
    add_msg_board(sprintf('Structure already contains %i labels',length(labels)));
end;

allowed=0;
% check whether sites are already labelled and whether all restraint sites
% do exist
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
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        ind1=resolve_address(adr1);
        if isempty(ind1),
            add_msg_board(sprintf('Warning: Restraint %i removed as it has first label at site %s',k,adr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        end;
        ind1(3)=mnum;
        adr1=mk_address(ind1);
        found=false;
        for l=1:length(labels),
            diff=ind1-lindices(l,:);
            if sum(abs(diff))==0 && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,adr1));
            end;
        end;
        adr2=restraints.DEER(k).adr2;
        ind2=resolve_address(adr2);
        if isempty(ind2),
            add_msg_board(sprintf('Warning: Restraint %i removed as it has second label at site %s',k,adr2));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        end;
        ind2(3)=mnum;
        adr2=mk_address(ind2);
        found=false;
        for l=1:length(labels),
            diff=ind2-lindices(l,:);
            if sum(abs(diff))==0  && strcmpi(labels(l).name,restraints.DEER(k).label),
                found=true;
            end;
        end;
        if ~found,
            for l=1:length(to_do_list),
                if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr2;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label %s at site %s will be generated.',restraints.DEER(k).label,adr2));
            end;
        end;
        allowed=allowed+1;
        DEER(allowed)=restraints.DEER(k);
    end;
    if ~isempty(to_do_list),
        add_msg_board('Warning: New labels are computed in a structure that may already contain attached spin labels.');
        add_msg_board('This may lead to erroneous rotamer distributions.');
        add_msg_board('Consider starting from an unlabeled template or from a template where rotamers are computed but labels are not attached.');
    end;
else
    poi=0;
    to_do_list{1}=' ';
    label_list{1}=' ';
    for k=1:length(restraints.DEER),
        adr1=restraints.DEER(k).adr1;
        ind1=resolve_address(adr1);
        if isempty(ind1),
            add_msg_board(sprintf('Warning: Restraint %i removed as it has first label at site %s',k,adr1));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        else
            ind1(3)=mnum;
            adr1=mk_address(ind1);
            found=false;
            for l=1:length(to_do_list),
                if strcmp(adr1,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr1;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr1));
            end;
        end;
        adr2=restraints.DEER(k).adr2;
        ind2=resolve_address(adr2);
        if isempty(ind2),
            add_msg_board(sprintf('Warning: Restraint %i removed as it has second label at site %s',k,adr2));
            add_msg_board(sprintf('This site does not exist in current template %s',mk_address(snum)));
            continue
        else
            ind2(3)=mnum;
            adr2=mk_address(ind2);
            found=false;
            for l=1:length(to_do_list),
                if strcmp(adr2,to_do_list{l}) && strcmpi(label_list(l),restraints.DEER(k).label),
                    found=true;
                end;
            end;
            if ~found,
                poi=poi+1;
                to_do_list{poi}=adr2;
                label_list{poi}=restraints.DEER(k).label;
                T_list(poi)=restraints.DEER(k).T;
                add_msg_board(sprintf('Rotamers for label at site %s will be generated.',adr2));
            end;
        end;
        allowed=allowed+1;
        DEER(allowed)=restraints.DEER(k);
    end;
end;
if allowed<length(restraints.DEER),
    add_msg_board(sprintf('Warning: %i out of %i DEER restraints had to be removed.',...
        length(restraints.DEER)-allowed,length(restraints.DEER)));
else
    add_msg_board('All DEER restraints can be kept.');
end;

restraints.DEER=DEER(1:allowed);

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('rotamers %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
    end;
end;

for k=1:length(to_do_list),
    if ~strcmp(to_do_list{k},' '),
        command=sprintf('label %s %s %i',to_do_list{k},label_list{k},T_list(k));
        hMain.store_undo=false;
        cmd(hMain,command);
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
            labels(poi).rmsd=NOpos_rmsd(NOpos);
        end;
    end;
end;

function [direct,cancelled]=process_direct_restraints(snum,mnum,restraints)

global model

cancelled=false;
if ~isfield(restraints,'direct'),
    direct=[];
    return;
end;

if ~isempty(alignment),
end;

md=length(restraints.direct);
direct=zeros(md,5);

[network,masses,cindices]=coarse_residues(mk_address(snum));

[mn,nn]=size(cindices);
poi=0;
add_msg_board(sprintf('Checking and processing %i direct C_alpha-C_alpha restraints',md));
for k=1:md,
    possible=true;
    adr1=restraints.direct(k).adr1;
    ind1=resolve_address(adr1);
    if isempty(ind1),
        possible=false;
        add_msg_board(sprintf('Warning: Specified residue %s does not exist in template structure.',adr1));
        continue
    end;
    ind1(3)=mnum;
    net1=0;
    for kk=1:mn,
        match=sum(abs(ind1-cindices(kk,:)));
        if match==0,
           net1=kk;
           break;
        end;
    end;
    adr2=restraints.direct(k).adr2;
    ind2=resolve_address(adr2);
    if isempty(ind2),
        possible=false;
        add_msg_board(sprintf('Warning: Specified residue %s does not exist in template structure.',adr2));
        continue;
    end;
    ind2(3)=mnum;
    net2=0;
    for kk=1:mn,
        match=sum(abs(ind2-cindices(kk,:)));
        if match==0,
           net2=kk;
           break;
        end;
    end;
    if net1==0 || net2==0,
        possible=false;
        add_msg_board(sprintf('Warning: Constrained C_alpha-C_alpha distance %s to %s does not exist in template structure.',adr1,adr2));
    end;
    if possible,
        poi=poi+1;
        direct(poi,1)=net1;
        direct(poi,2)=net2;    
        direct(poi,4)=restraints.direct(k).r;
        direct(poi,5)=restraints.direct(k).sigr;
    end;
end;
if poi<md,
    add_msg_board(sprintf('Warning: %i out of %i direct restraints had to be removed.',md-poi,md));
else
    add_msg_board('All direct restraints can be kept.');
end;
direct=direct(1:poi,:);

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


function message=mk_modeller_input(runfile,alignfile,idCode,snum,mnum,restraints,DEER,direct,ares)

global model

tightness=0.085; % tightness of DEER distance restraint fitting, 0.085 is optimized
%tightness=tightness*sqrt(2);

message.error=0;
message.text='No error.';

if isfield(restraints,'helices'),
    helices=restraints.helices;
else
    helices=[];
end;
if isfield(restraints,'strands'),
    strands=restraints.strands;
else
    strands=[];
end;
if isfield(restraints,'sheets'),
    sheets=restraints.sheets;
else
    sheets=[];
end;
chain_tags=model.chain_tags{snum};
if length(model.structures{snum})==1,
    single_chain=true;
else
    single_chain=false;
end;



fid=fopen(runfile,'wt');
if fid==-1,
    message.error=2;
    message.text='Modeller input file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'from modeller import *\n');
fprintf(fid,'from modeller.automodel import *\n\n');
fprintf(fid,'log.verbose()\n');
fprintf(fid,'env = environ(rand_seed=-8123, restyp_lib_file=''$(LIB)/restyp_EPR.lib'', copy=None)\n\n');

fprintf(fid,'env.io.atom_files_directory = [''.'', ''../atomfiles'']\n\n');

fprintf(fid,'env.io.hetatm = True\n\n');

fprintf(fid,'class MyModel(automodel):\n');
fprintf(fid,'    def special_restraints(self, aln):\n');
fprintf(fid,'        rsr = self.restraints\n');
fprintf(fid,'        at  = self.atoms\n');

if ~isempty(DEER),
    % make list of labeled residues
    reslist=zeros(1,400);
    pairlist=zeros(length(DEER),2);
    poi=0;
    for k=1:length(DEER),
        [stag,ctag1,modelnum,resnum1]=mk_address_parts(DEER(k).indices(1,:));
        if DEER(k).indices(1,2)==1,
            resnum1=resnum1-ares+1;
        end;
        [stag,ctag2,modelnum,resnum2]=mk_address_parts(DEER(k).indices(2,:));
        if DEER(k).indices(1,2)==1,
            resnum2=resnum2-ares+1;
        end;
        pairlist(k,:)=[resnum1 resnum2];
        if isempty(find(reslist(1:poi)==resnum1,1)),
            poi=poi+1;
            reslist(poi)=resnum1;
            chain_assign{poi}=ctag1;
        end;
        if isempty(find(reslist(1:poi)==resnum2,1)),
            poi=poi+1;
            reslist(poi)=resnum2;
            chain_assign{poi}=ctag2;
        end;
    end;
    reslist=sort(reslist(1:poi));
    for k=1:length(reslist);
        if single_chain,
            fprintf(fid,'        NO%i = at[''N1:%i''], at[''O1:%i'']\n',reslist(k),reslist(k),reslist(k));
        else
            fprintf(fid,'        NO%i = at[''N1:%i:%s''], at[''O1:%i:%s'']\n',reslist(k),reslist(k),chain_assign{k},reslist(k),chain_assign{k});
        end;
        fprintf(fid,'        label%i = pseudo_atom.gravity_center(NO%i)\n',reslist(k),reslist(k));
        fprintf(fid,'        rsr.pseudo_atoms.append(label%i)\n',reslist(k));
    end;
    fprintf(fid,'\n');
    
    for k=1:length(reslist);
        if single_chain,
            fprintf(fid,'	r%i = rigid_body(self.residue_range(''%i'', ''%i''))\n',reslist(k),reslist(k),reslist(k));
        else
            fprintf(fid,'	r%i = rigid_body(self.residue_range(''%i:%s'', ''%i:%s''))\n',reslist(k),reslist(k),chain_assign{k},reslist(k),chain_assign{k});
        end;
        fprintf(fid,'	rsr.rigid_bodies.append(r%i)\n',reslist(k));
    end;
    fprintf(fid,'\n');

    for k=1:length(DEER);
        fprintf(fid,'	rsr.add(forms.gaussian(group=physical.xy_distance,\n');
        fprintf(fid,'                               feature=features.distance(label%i,\n',pairlist(k,1));
        fprintf(fid,'                                                         label%i),\n',pairlist(k,2));
        fprintf(fid,'                               mean=%6.3f, stdev=%5.3f))\n',10*DEER(k).r,tightness);
    end;
    fprintf(fid,'\n');
end;

if ~isempty(helices),
    for k=1:length(helices);
        [cid,resnum1,resnum2]=dissect_address(helices(k).adr,chain_tags);
        if isempty(cid), cid=1; end;
        if cid==1,
            resnum1=resnum1-ares+1;
            resnum2=resnum2-ares+1;
        end;
        ctag1=id2tag(cid,chain_tags);
        ctag2=ctag1;
        fprintf(fid,'	rsr.add(secondary_structure.alpha(self.residue_range(');
        if single_chain,
            fprintf(fid,'''%i:'',''%i:'')))\n',resnum1,resnum2);
        else
            fprintf(fid,'''%i:%s'',''%i:%s'')))\n',resnum1,ctag1,resnum2,ctag2);
        end;
    end;
    fprintf(fid,'\n');
end;

if ~isempty(strands),
    for k=1:length(strands);
        [cid,resnum1,resnum2]=dissect_address(strands(k).adr,chain_tags);
        if isempty(cid), cid=1; end;
        if cid==1,
            resnum1=resnum1-ares+1;
            resnum2=resnum2-ares+1;
        end;
        ctag1=id2tag(cid,chain_tags);
        ctag2=ctag1;
        fprintf(fid,'	rsr.add(secondary_structure.strand(self.residue_range(');
        if single_chain,
            fprintf(fid,'''%i:'',''%i:'')))\n',resnum1,resnum2);
        else
            fprintf(fid,'''%i:%s'',''%i:%s'')))\n',resnum1,ctag1,resnum2,ctag2);
        end;
    end;
    fprintf(fid,'\n');
end;

if ~isempty(sheets),
    for k=1:length(sheets);
        ind1=resolve_address(sheets(k).adr1);
        ind2=resolve_address(sheets(k).adr2);
        [stag,ctag1,modelnum,resnum1]=mk_address_parts(ind1);
        if ind1(2)==1,
            resnum1=resnum1-ares+1;
        end;
        [stag,ctag2,modelnum,resnum2]=mk_address_parts(ind2);
        if ind2(2)==1,
            resnum2=resnum2-ares+1;
        end;
        fprintf(fid,'	rsr.add(secondary_structure.sheet(at[,\n');
        if single_chain,
            fprintf(fid,'''N:%i''], at[''O:%i''],\n',resnum1,resnum2);
            fprintf(fid,'                               sheet_h_bonds=%i))\n',str2double(sheets(k).length));
        else
            fprintf(fid,'''N:%i:%s''], at[''O:%i:%s''],\n',resnum1,ctag1,resnum2,ctag2);
            fprintf(fid,'                               sheet_h_bonds=%i))\n',str2double(sheets(k).length));
        end;
    end;
    fprintf(fid,'\n');
end;

if ~isempty(direct),
    for k=1:length(direct);
        ind1=resolve_address(direct(k).adr1);
        ind2=resolve_address(direct(k).adr2);
        [stag,ctag1,modelnum,resnum1]=mk_address_parts(ind1);
         if ind1(2)==1,
            resnum1=resnum1-ares+1;
        end;
       [stag,ctag2,modelnum,resnum2]=mk_address_parts(ind2);
        if ind2(2)==1,
            resnum2=resnum2-ares+1;
        end;
        fprintf(fid,'	rsr.add(forms.gaussian(group=physical.xy_distance,\n');
        if single_chain,
            fprintf(fid,'                               feature=features.distance(at[''CA:%i''],\n',resnum1);
            fprintf(fid,'                                                         at[''CA:%i'']),\n',resnum2);
        else
            fprintf(fid,'                               feature=features.distance(at[''CA:%i:%s''],\n',resnum1,ctag1);
            fprintf(fid,'                                                         at[''CA:%i:%s'']),\n',resnum2,ctag2);
        end;
        fprintf(fid,'                               mean=%6.3f, stdev=0.020))\n',10*direct(k).r);
    end;
    fprintf(fid,'\n');
end;

fprintf(fid,'a = MyModel(env, alnfile = ''%s'',\n',alignfile);
fprintf(fid,'            knowns = ''%s'', sequence = ''%s'',\n','temp',idCode);
fprintf(fid,'            assess_methods=(assess.DOPE, assess.GA341))\n');

fprintf(fid,'a.starting_model= %i\n',mnum);
fprintf(fid,'a.ending_model  = %i\n\n',mnum);

fprintf(fid,'a.make()\n');

fclose(fid);

function [cid,resnum,resnum2]=dissect_address(adr,chain_tags)

resnum=[];
resnum2=[];
cid=[];
cprefix=strfind(adr,'(');
csuffix=strfind(adr,')');
hyphen=strfind(adr,'-');
if ~isempty(csuffix),
    if csuffix<length(adr),
        if isempty(hyphen),
            resnum=str2double(adr(csuffix+1:end));
        else
            resnum=str2double(adr(csuffix+1:hyphen-1));
            if hyphen<length(adr),
                resnum2=str2double(adr(hyphen+1:end));
            end;
        end;
    else
        return;
    end;
else
    if isempty(hyphen),
        resnum=str2double(adr);
    else
        resnum=str2double(adr(1:hyphen-1));
        if hyphen<length(adr),
            resnum2=str2double(adr(hyphen+1:end));
        end;
    end;
end;
if isnan(resnum),
    resnum=[];
    return;
end;
if isnan(resnum2),
    resnum2=[];
    return;
end;

if ~isempty(cprefix) && ~isempty(csuffix),
    ctag=adr(cprefix+1:csuffix-1);
    cid=tag2id(ctag,chain_tags);
end;

function make_report(snum,DEER)

global model 
global general
global hMain

modnum=length(model.structures{snum}(1).xyz);

filename=sprintf('%s_restraint_matching.txt',id2tag(snum,model.structure_tags));
fname=fullfile(general.tmp_files, filename);
fid=fopen(fname,'wt');
if fid==-1,
    add_msg_board('Restraint matching report file could not be written');
else
    fprintf(fid,'Restraint matching report for structure %s\n',id2tag(snum,model.structure_tags));
    fprintf(fid,'(all distances in nm)\n');
end;

full_rmsd=0;
full_num=0;
for km=1:modnum,
    if fid~=-1,
        fprintf(fid,'\n--- Restraint matching for model %i ---\n\n',km);
        fprintf(fid,'Resid. 1  Resid. 2   r(exp)  sr(exp) r(mod)  |dr/sr|  Matching\n');
    end;
    ma=0;
    rmsd=0;
    num=0;
    fom=0; % figure of merit
    if ~isempty(DEER),
        for k=1:length(DEER),
            lindices=DEER(k).indices;
            res1adr=mk_address([snum lindices(1,2) km lindices(1,4)]);
            [stag,ctag1,modelnum,resnum1]=mk_address_parts([snum lindices(1,2) km lindices(1,4)]);
            r1adrshort=sprintf('(%s)%i',ctag1,resnum1);
            [msg,N1coor]=get_object([res1adr '.N1'],'coor');
            [msg,O1coor]=get_object([res1adr '.O1'],'coor');
            DEER(k).xyz1=(N1coor+O1coor)/2;
            res2adr=mk_address([snum lindices(2,2) km lindices(2,4)]);
            [stag,ctag2,modelnum,resnum2]=mk_address_parts([snum lindices(2,2) km lindices(2,4)]);
            r2adrshort=sprintf('(%s)%i',ctag2,resnum2);
            [msg,N1coor]=get_object([res2adr '.N1'],'coor');
            [msg,O1coor]=get_object([res2adr '.O1'],'coor');
            DEER(k).xyz2=(N1coor+O1coor)/2;
            r0=norm(DEER(k).xyz1-DEER(k).xyz2)/10;
            rmsd=rmsd+(r0-DEER(k).r)^2;
            num=num+1;
            full_num=full_num+1;
            det=abs(r0-DEER(k).r)/DEER(k).sigr;
            fom=fom+det^2;
            if det>2,
                match_str='### not matched ###';
            elseif det>1,
                match_str='*** poorly matched ***';
            else
                match_str='well matched';
            end;
            if fid~=-1,
                fprintf(fid,'%10s%10s%8.2f%8.2f%8.2f %8.2f    %s\n',r1adrshort,r2adrshort,DEER(k).r,DEER(k).sigr,r0,det,match_str);
            end;
            if DEER(k).r+DEER(k).sigr>ma,
                ma=DEER(k).r+DEER(k).sigr;
            end;
            if r0>ma,
                ma=r0;
            end;
        end;
        if num==0, num=1; end;
        fprintf(fid,'\nMean r.m.s.d. of all restraints for model %i: %5.2f nm\n',km,sqrt(rmsd/num));
    end;
    full_rmsd=full_rmsd+rmsd;
end;
if fid~=-1,
    full_rmsd=sqrt(full_rmsd/full_num);
    fprintf(fid,'\nMean r.m.s.d. of all restraints for all models: %5.2f nm\n',full_rmsd);
    fclose(fid);
    hMain.report_file=fname;
    report_editor;
end;

function [batname,message]=mk_modeller_bat(modpath,modcmd,jobfile,is_bckg)

global general

if nargin<4,
    is_bckg=false;
end;

batname='run_modeller.bat';
envfile=[modpath filesep 'modenv.bat'];
batfile=[modpath filesep  batname];

message.error=0;
message.text='No error.';

fid=fopen(batfile,'wt');
if fid==-1,
    message.error=2;
    message.text='Modeller batch file could not be written.';
    add_msg_board('ERROR: File could not be opened for writing.');
    add_msg_board(message.text);
    return;
end;

rfid=fopen(envfile,'rt');
if rfid==-1,
    message.error=3;
    message.text='Modeller environment file could not be read.';
    add_msg_board('ERROR: Modeller environment file could not be read.');
    add_msg_board(message.text);
    return;
end;

while 1
    tline = fgetl(rfid);
    if ~ischar(tline), break, end
    fprintf(fid,'%s\n',tline);
end;

fprintf(fid,'echo modeller installed\n');
k=strfind(general.tmp_files,':');
if ~isempty(k),
    fprintf(fid,'%s\n',general.tmp_files(1:k));
end;

fprintf(fid,'cd %s\n',general.tmp_files);
if is_bckg,
    fprintf(fid,'start %s %s\n',modcmd,jobfile);
else
    fprintf(fid,'%s %s\n',modcmd,jobfile);
end;

fclose(fid);

