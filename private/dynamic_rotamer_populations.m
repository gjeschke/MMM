function [rotamer_stats,message]=dynamic_rotamer_populations(adr,backbone_stats,T,label,MD2PDB,calibration,radius)
% function rotamer_stats=dynamic_rotamer_populations(adr,T,label,MD2PDB,calibration,radius)
%    Updates populations in rotamer statistics based on SCWRL4 repacking
%    energies
%    SCWRL4 must be installed
%
% adr       residue address, must address a single residue, can also be a
%           1x4 index vector for a residue
% backbone_stats     rotamer coordinates and populations as provided by
%                   get_rotamers assuming only interactions with backbone
%                   atoms
% T                 labeling temperature
% label             type of label
% MD2PDB            translation table between MD and PDB atom order
% calibration       rotamer library calibration
% radius            (optional) radius around the labeling site C_alpha atom
%                   in which SCWRL4 is requested to repack sidechains (Å),
%                   defaults to infinity (whole structure is repacked)
%
% message   error message with fields message.error (integer, 0 no error)
%           and message.text (error description)
%
% G. Jeschke, 2010


global model
global general
global help_files
global third_party
global residue_defs

rotamer_stats.all_potentials=[]; % an empty structure should be returned when an error is encountered

scale_SCWRL4=0.2;

origin=backbone_stats.loc_frame_Ca; % Ca coordinates of labeling site

int_pop=calibration.pop;    % get non-normalized internal populations from the rotamers
int_pop=int_pop/sum(calibration.pop);   % normalize populations


conv_factor=(4.1868*1000); % conversion kcal/mol to J/mol to SI
gas_un=8.314472;    % universal gas constant in SI (J/(mol*K)       

message.error=0;
message.text='Dynamic rotamer populations are computed.';

if isa(adr,'char'),
    [indices,message]=resolve_address(adr);
    if message.error,
        message.text='ERROR: Residue address could not be resolved.';
        return;
    end;
else
    indices=adr;
end;

[m,n]=size(indices);
if m==1,
    indices=indices(indices>0);
    [m,n]=size(indices);
end;
if m~=1 || n~=4,
    message.error=1;
    message.text='Addressed object is not a single residue.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

entry=strcat(help_files,'third_party.html#SCWRL4');

dospath=which('scwrl4.exe');
if isempty(dospath),
    message.error=5;
    message.text='SCWRL4 software not found on Matlab path.';
    add_msg_board('Sidechain repacking requires SCWRL4 by Krivov, Shapovalov, and Dunbrack');
    add_msg_board('ERROR: SCWRL4 could not be found on the Matlab path.');
    add_msg_board('Please check whether SCWRL4 is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end;

snum=indices(1);
mnum=1;
snum0=model.current_structure;
snum=indices(1);
model.current_structure=snum;
cnum=indices(2);
rnum=indices(4);

chains=length(model.structures{snum}(:));
longseq='';
mutated=0;
for c=1:chains,
    residues=length(model.structures{snum}(c).residues{1}.info);
    newr=0;
    for r=1:residues,
        rtag=model.structures{snum}(c).residues{1}.info(r).name;
        atags=model.structures{snum}(c).residues{1}.info(r).atom_tags;
        backbone=true;
        id=tag2id('N',atags);
        if isempty(id), backbone=false; end;
        id=tag2id('O',atags);
        if isempty(id), backbone=false; end;
        id=tag2id('CA',atags);
        if isempty(id), 
            backbone=false;
            pointer=[];
        else
            pointer=model.structures{snum}(c).residues{1}.info(r).atom_numbers{id};
        end;
        id=tag2id('C',atags);
        if isempty(id), backbone=false; end;
        if backbone,
            newr=newr+1;
            cslc=tag2id(upper(rtag),upper(residue_defs.restags),residue_defs.single_letter_code);
            if isempty(cslc),
                cslc='G';
            end;
            if c==cnum && r==rnum,
                cslc='g'; % the residue to be mutated is lowercase glycine
                mutated=newr;
            else
                cslc=upper(cslc); % other residues are uppercase
            end;
            if nargin>=7, % repacking radius is defined
                if ~isempty(pointer)
                    xyz=model.structures{snum}(c).xyz{1}(pointer(1,1),:); % only first location of C_alpha is considered
                    if norm(xyz-origin)>radius, % residues outside packing radius are not repacked
                        cslc=lower(cslc);
                    end;
                    % Bfactor=model.structures{snum}(c).Bfactor{1}(pointer(1,1));
                else
                    cslc=lower(cslc); % if C_alpha is missing, residue is not repacked
                end;
            end;
            longseq=[longseq cslc];
        end;
    end;
end;

if mutated==0,
    message.error=6;
    message.text='Residue to be mutated not found or has incomplete backbone.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return
end;

seqfile=[general.tmp_files 'tmp.seq'];
fid=fopen(seqfile,'wt');
if fid==-1,
    message.error=7;
    message.text='Sequence file could not be written.';
    add_msg_board('ERROR: No mutation performed.');
    add_msg_board(message.text);
    return;
end;

fprintf(fid,'%s\n',longseq);
fclose(fid);

add_msg_board('Saving PDB input file...');
set(gcf,'Pointer','watch');
drawnow;

outfile=[general.tmp_files 'tmp.pdb'];
message=wr_pdb(outfile,'MMMM',mnum,':',true);

if message.error,
    add_msg_board(message.text);
    return
end;

framefile=[general.tmp_files 'frame.pdb'];
add_msg_board('Saving PDB frame file (heteroatoms)...');
message=wr_het_frame(framefile,'MHET',mnum,':',true);

rotamers=length(backbone_stats.all_rotamers);
frame = fileread(framefile);
% rotamers=10;
energies=zeros(1,rotamers);

add_msg_board('Computing sidechain interaction energy...');
comp_status=status_figure('Rotamers: Close to stop.');

for rotnum=1:rotamers,
    drawnow
    if mod(rotnum,round(rotamers/100))==0    % just a counter
        comp_status=status_figure((rotnum-1)/rotamers,comp_status);
        if ~comp_status, 
            message.error=9;
            message.text='### Warning ### Rotamer analysis stopped.';
            return 
        end;
    end;    
%     add_msg_board(sprintf('Computing dynamic interaction energy for rotamer %i...',rotnum));
    rotamer=sprintf_rotamer_pdb(backbone_stats,rotnum,label,MD2PDB,false);
    cframe=[general.tmp_files 'cframe.pdb'];
    cfid=fopen(cframe,'w');
    fprintf(cfid,'%s%s',frame,rotamer);
    fclose(cfid);
    infile=[general.tmp_files 'scwrl4_out.pdb'];

    comd=[dospath ' -i ' outfile ' -o ' infile];
    if frame,
        comd=[comd ' -f ' cframe];
    end;
    comd=[comd ' -s ' seqfile ' -h -t'];
    [s, w] = dos(comd);
    en=1e12;
    if s~=0,
        message.error=8;
        message.text='SCWRL4 error.';
        add_msg_board('ERROR: SCWRL4 did not run successfully.');
        set(gcf,'Pointer','arrow');
        return
    else
        rem=w;
        while ~isempty(rem),
            [token,rem]=strtok(rem,char(10));
            if ~isempty(token),
                % add_msg_board(token);
                if strfind(token,'Total minimal energy of the graph'),
                    [comment,energy]=strtok(token,'=');
                    en=str2num(energy(2:end));
                end;
            end;
        end;
    end;
    energies(rotnum)=en*conv_factor;
%     add_msg_board(sprintf('Energy is %10.1f kJ/mol\n',conv_factor*en/1000));

end;


energies=energies-min(energies);
ext_pop=exp(-scale_SCWRL4*energies/(gas_un*T));
ext_pop=ext_pop.*backbone_stats.all_potentials.ext_pop*backbone_stats.all_potentials.ext_net; % add backbone interaction energy
partition_function=sum(ext_pop.*int_pop)/sum(int_pop);
ext_net=sum(ext_pop);
ext_pop=ext_pop/sum(ext_pop);
pop_rot=ext_pop.*int_pop; % get populations for all rotamers
pop_rot=pop_rot/sum(pop_rot); % normalize full rotamer populations


% figure(7); clf;
% plot(backbone_stats.all_potentials.ext_pop,'k.');
% hold on
% plot(ext_pop,'ro');
% set(gca,'FontSize',14);
% axis([-5,rotamers+5,-0.005,max([max(ext_pop),max(backbone_stats.all_potentials.ext_pop)])+0.005]);

set(gcf,'Pointer','arrow');

rotamer_stats.all_potentials.ext_pop=ext_pop;
rotamer_stats.all_potentials.pop_rot=pop_rot;
rotamer_stats.all_potentials.partition_function=partition_function;
rotamer_stats.all_potentials.ext_net=ext_net;
rotamer_stats.ext_poten_type='dynamic_SCWRL4';

NOstats_all=backbone_stats.NOall;
NOstats_all(:,4)=pop_rot;
rotamer_stats.NOall=NOstats_all; % NO-centers and weights for all rotamers in a single matrix;
rotamer_stats.all_rotamers=backbone_stats.all_rotamers; % contains coordinates for all rotated rotamers for current position
rotamer_stats.labels_own=backbone_stats.labels_own; % rotamer coordinates in the local residue frame (Ca is always [0,0,0]);
rotamer_stats.loc_frame_Ca=backbone_stats.loc_frame_Ca; % global (protein) coordinates of the local frame origin (Ca alpha of the mutated residue);
        
drawnow;

% add the reference, if it does not yet exist
scwrl4_ref=true;
id=tag2id('Krivov:2009_scwrl4',third_party.tags,[],'|');
if ~isempty(id),
    if isfield(model,'auto_references'),
        if ~isempty(find(id==model.auto_references, 1)),
            scwrl4_ref=false;
        end;
    else
        model.auto_references=[];
    end;
    if scwrl4_ref,
        if ~isfield(model,'references'),
            model.references(1)=third_party.references(id);
        elseif isempty(model.references)
            model=rmfield(model,'references');
            model.references(1)=third_party.references(id);
        else
            model.references(end+1)=third_party.references(id);
        end;
        model.auto_references(end+1)=id;
    end;
end;

model.current_structure=snum0;

message.error=0;
message.text='Dynamic rotamer populations are computed.';

if comp_status, status_figure(1); end;

