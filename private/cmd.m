function handles=cmd(handles,command_line)
% function handles=cmd(handles,command_line)
%
% The command interpreter of MMM
% parses input commands and executes them
%
% the flag store_undo determines whether the command is stored in the undo
% list (default 1, store), used with store_undo=0 when cmd is called from
% the undo function

global hMain

if get(hMain.checkbox_log,'Value')
    fid=fopen(hMain.logfile,'a+');
    fprintf(fid,'cmd> %s\n',command_line);
    fclose(fid);
end

set(hMain.MMM,'Pointer','watch');
drawnow

commands=':atompair:attach:beacons:bckg:bilabel:blscan:camup:color:colorscheme:compact:copy:delete:detach:dihedrals:distance:domain:download:dssp:echo:helix:help:hide:inertiaframe:label:libcomp:libtest:lock:locrmsd:locorder:loop:mass:motion:mushroom:new:ortho:pdbload:persp:plot:radgyr:redo:refrmsd:remodel:repack:replace:report:rmsd:rotamers:SAS:scopy:select:sheet:show:statistics:symmetry:synonym:transparency:undo:undefine:unlock:unselect:view:wrseq:zoom:';

[cmd,args]=strtok(command_line); % separate command and arguments

test=tag2id(cmd,commands);

if isempty(test)
    % allow for abbreviated commands
    poi=strfind(commands,cmd);
    if length(poi)==1, % the abbreviated command is unique
        cmd=strtok(commands(poi:end),':');
    elseif length(poi)>1,
        poi0=poi;
        pp=0;
        for k=1:length(poi0),
            if poi0(k)>1 && strcmp(commands(poi0(k)-1),':'),
                pp=pp+1;
                poi(pp)=poi0(k);
            end;
        end;
        if pp>0,
            poi=poi(1:pp);
        else
            poi=[];
        end;
    end;
    if length(poi)>1,
        add_msg_board(sprintf('ERROR: Command is not unique. "%s" could be:',cmd));
        for k=1:length(poi),
            pcmd=strtok(commands(poi(k):end),':');
            add_msg_board(sprintf('%s',pcmd));
        end;
        add_msg_board('No action taken.');
    elseif isempty(poi),
        add_msg_board(sprintf('ERROR: Command "%s" is not a valid command',cmd));
    end;
end;
% List of all recognized commands, alternative spellings allowed,
% please keep alphabetic order of primary spelling
switch cmd,
    case 'atompair'
        handles=atompair(handles,args);
    case 'attach'
        handles=attach(handles);
    case {'background','bckg'}
        handles=set_background(handles,args);
    case 'beacons'
        handles=get_beacons(handles,args);
    case 'bilabel'
        handles=bilabel(handles,args);
    case 'blscan'
        handles=blscan(handles,args);
    case 'camup'
        handles=def_camup(handles,args);
    case 'color' 
        handles=color(handles,args);
    case 'colorscheme' 
        handles=colorscheme(handles,args);
    case 'compact'
        handles = ensemble_compact(handles,args);
    case 'copy' 
        handles=copy_3D(handles);
    case 'delete'
        handles=delete_cmd(handles,args);
    case 'detach'
        handles=detach(handles);
    case 'dihedrals'
        handles=dihedrals(handles,args);
    case 'distance'
        handles=distance(handles,args);
    case 'domain'
        handles=def_domain(handles,args);
    case 'download'
        handles=download(handles,args);
    case 'dssp'
        handles=dssp_assign(handles,args);
    case 'echo'
        handles=echo(handles,args);
    case 'helix'
        handles=def_helix(handles,args);
    case 'help'
        handles=help_cmd(handles,args,commands);
    case 'hide' 
        handles=hide(handles,args);
    case 'inertiaframe'
        inertiaframe;
    case 'label'
        handles=label(handles,args);
    case 'libcomp'
        handles=libcomp(handles,args);
    case 'libtest'
        handles=libtest(handles,args);
    case 'lock'
        handles=lock(handles,args);
    case 'locorder'
        handles=local_order(handles,args);
    case 'locrmsd'
        handles=local_rmsd(handles,args);
    case 'loop'
        handles=def_loop(handles,args);
    case 'mass' 
        handles=mass(handles,args);
    case 'motion' 
        handles=motion(handles,args);
    case 'mushroom' 
        handles=mushroom(handles,args);
    case 'new' 
        handles=new_model(handles,args);
    case {'orthographic','ortho'}
        handles=orthographic(handles);
    case 'pdbload'
        handles=pdbload(handles,args);
    case {'perspective','persp'}
        handles=perspective(handles);
    case 'plot'
        handles=plot_cmd(handles,args);
    case 'radgyr'
        handles=ensemble_rg(handles,args);
    case 'refrmsd'
        handles=ensemble_reference_rmsd(handles,args);
    case 'remodel'
        handles=remodel_section(handles,args);
    case 'repack'
        handles=repack_structure(handles,args);
    case 'report'
        handles=report(handles,args);
    case 'replace'
        handles=replace_restype(handles,args);
    case 'rmsd'
        handles=ensemble_rmsd(handles,args);
    case 'rotamers'
        handles=rotamers(handles,args);
    case 'SAS'
        handles=SAS(handles,args);
    case 'scopy'
        handles=scopy(handles,args);
    case 'select'
        handles=select(handles,args);
    case 'sheet'
        handles=def_sheet(handles,args);
    case 'show' 
        handles=show(handles,args);
    case 'statistics'
        handles = coil_statistics(handles,args);
    case 'symmetry'
        handles=symmetry(handles);
    case 'synonym'
        handles=mk_synonym(handles,args);
    case 'transparency' 
        handles=transparency(handles,args);
    case 'uncolor' 
        handles=uncolor(handles,args);
    case 'undefine' 
        handles=undefine(handles,args);
%     case 'unlabel'
%         handles=unlabel(handles,args);
    case 'unlock'
        handles=unlock(handles,args);
    case 'unselect'
        handles=unselect(handles,args);
    case 'untransparency' 
        handles=untransparency(handles,args);
%     case 'xyz' 
%         handles=xyz(handles,args);
    case 'view'
        handles=def_view(handles,args);
    case 'wrseq'
        handles=write_sequence(handles,args);
    case 'zoom'
        handles=zoom(handles,args);
    case 'null'
        add_msg_board('Null command for unnecessary undo');
end;

set(handles.edit_command_line,'String',command_line);

set(hMain.MMM,'Pointer','arrow');

% Command execution functions
% --------------------------
% these may also be defined in other files, all of them have to accept and
% return variable 'handles', they may (and usually will) modify handles
% optionally they may accept arguments 'args'
%
% Undo/Redo issue:
% Commands should generally store an inverse command in the undo list and
% record the command in the redo list, using function 'cmd_history' defined
% in a separate file, if no inverse command exists, they should call
% cmd_history with empty or no inverse command argument to produce a user warning 
% and execute only when the user accepts that, the command is rejected when
% the veto flag is set:
%
% [handles,veto]=cmd_history(handles,command,undo_cmd)
%
% 'undo' and 'redo' themselves, as well as file input and output are not
% executed by the command interpreter and are not stored in the history
% variables

function handles = atompair(handles,args)

global hMain
global rotamer_libraries

command=sprintf('atompair %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: atompair adr1 adr2');
    add_msg_board('Only atom addresses are allowed.');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) ~= 2
    add_msg_board('ERROR: Command requires two arguments.');
    return
end;


adr1=char(myargs{1}(1));
adr2=char(myargs{1}(2));

[msg1,xyz1] = get_object(adr1,'coor');

if msg1.error ~= 0 || isempty(xyz1)
    add_msg_board(sprintf('ERROR: First atom coordinate not found (%s)',msg1.text));
    return
end

[msg2,xyz2] = get_object(adr2,'coor');
if msg2.error ~= 0 || isempty(xyz2)
    add_msg_board(sprintf('ERROR: Second atom coordinate not found (%s)',msg2.text));
    return
end

dist = norm(xyz1-xyz2);

add_msg_board(sprintf('Distance %s - %s is %5.2f Å',adr1,adr2,dist));


function handles=set_background(handles,args)
% Sets background of model window
global hMain
global hModel

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: bckg black|grey|white');
    return
end;

command=sprintf('bckg %s',args);
switch hMain.color,
    case 'black'
        undo_cmd=sprintf('bckg black');
    case {'grey','gray'}
        undo_cmd=sprintf('bckg grey');
    case 'white'
        undo_cmd=sprintf('bckg white');
end;

bckg=strtok(args);

switch bckg,
    case {'black','k'}
        [handles,veto]=cmd_history(handles,command,undo_cmd);
        hMain.color='black';
        if hMain.detached,
            set(hModel.figure,'Color','k');
        else
            set(handles.panel_model,'BackgroundColor','k');
            set(handles.panel_model,'ForegroundColor',[216,41,0]/255);
            set(handles.panel_model,'HighlightColor',[255,177,100]/255);
            set(handles.panel_model,'ShadowColor',[222,125,0]/255);
        end;
    case {'grey','gray'}
        [handles,veto]=cmd_history(handles,command,undo_cmd);
        hMain.color='grey';
        if hMain.detached,
            set(hModel.figure,'Color',[0.941,0.941,0.941]);
        else
            set(handles.panel_model,'BackgroundColor',[0.941,0.941,0.941]);
            set(handles.panel_model,'ForegroundColor','k');
            set(handles.panel_model,'HighlightColor','w');
            set(handles.panel_model,'ShadowColor',[128,128,128]/255);
        end;
    case {'white','w'}
        [handles,veto]=cmd_history(handles,command,undo_cmd);
        hMain.color='white';
        if hMain.detached,
            set(hModel.figure,'Color','w');
        else
            set(handles.panel_model,'BackgroundColor','w');
            set(handles.panel_model,'ForegroundColor','k');
            set(handles.panel_model,'HighlightColor','w');
            set(handles.panel_model,'ShadowColor',[128,128,128]/255);
        end;
    otherwise,
        add_msg_board('Warning: Background color can be only black, grey, or white');
end;

function handles=color(handles,args)

global graph_settings

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: color address SVG_color_name');
    add_msg_board('or   : color address red green blue');
    add_msg_board('with red, green, blue between 0 and 1.');
    return
end;

command=sprintf('color %s',strtrim(args));
undo_cmd=sprintf('uncolor %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    colname=char(myargs{1}(2));
else
    add_msg_board('ERROR: color argument missing.');
    return
end;
if length(myargs{1})>=4,
    red=str2double(colname);
    green=str2double(char(myargs{1}(3)));
    blue=str2double(char(myargs{1}(4)));
    if ~isnan(red) && ~isnan(green) && ~isnan(blue),
        rgb=[red green blue];
        if max(rgb) > 1
            rgb = rgb/255;
        end
    else
        rgb=[];
    end;
else
    colind=tag2id(colname,graph_settings.color_tags);
    if ~isempty(colind),
        rgb=graph_settings.colors(colind,:);
    else
        rgb=[];
        add_msg_board('ERROR: Named color does not exist or less than three RGB values.');
        add_msg_board('--- SVG colors ---');
        msg='';
        fail=0;
        id=0;
        while ~fail,
            id=id+1;
            com=id2tag(id,graph_settings.color_tags);
            if isempty(com),
                fail=1;
                add_msg_board(msg);
            else
                if length(msg)+length(com)<63,
                    if isempty(msg), 
                        msg=com;
                    else
                        msg=sprintf('%s %s',msg,com);
                    end;
                else
                    add_msg_board(msg);
                    msg=com;
                end;
            end;
        end;
    end;
end;
if ~isempty(rgb),
    if sum(rgb>1)+sum(rgb<0)==0,
        if strcmp(address(1),'$'), % surface request
            msg=set_surface(address,'color',rgb);
            if msg.error,
                if iscell(msg.text)
                    for k=1:length(msg.text),
                        add_msg_board(msg.text{k});
                    end,
                else
                    add_msg_board(msg.text);
                end;
            end;
        else
            [message,argout]=set_object(address,'color',{rgb});
        end;
        % special handling for selected objects
        if strcmp(strtrim(address),'*')
            msg=sprintf('Color of selected objects was changed to [%5.3f,%5.3f,%5.3f].',rgb);
            add_msg_board(msg);
            highlight_selection
        end;
    else
        add_msg_board('ERROR: Values of the RGB vector must be between 0 and 1.');  
    end;
end;

function handles=download(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: download pdb_identifier');
    return
end;

command=sprintf('download %s',strtrim(args));
undo_cmd=' ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
pdbid=strtrim(char(myargs{1}(1)));
[msg,snum]=download_pdb(pdbid);
if isempty(snum),
    add_msg_board('No structure downloaded');
else
    add_msg_board(sprintf('Structure [%s] downloaded as structure number %i',pdbid,snum));
end;
if msg.error,
    add_msg_board(msg.text);
end;
drawnow;


function handles=pdbload(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: pdbload filename');
    return
end;

command=sprintf('pdbload %s',strtrim(args));
undo_cmd=' ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
fname=strtrim(char(myargs{1}(1)));
poi=strfind(fname,'.');
if isempty(poi),
    add_msg_board('File extension missing. Trying .pdb, .ent in this sequence.');
    tname=[fname '.pdb'];
    tname=which(tname);
    if isempty(tname),
        tname=[fname '.ent'];
        tname=which(tname);
        if isempty(tname),
            add_msg_board(sprintf('ERROR: Specified PDB file %s not found on current path.',tname));
            return
        end;
    end;
    fname=tname;
end;
[msg,snum]=load_pdb(fname);
if isempty(snum),
    add_msg_board('No structure loaded');
else
    add_msg_board(sprintf('Local file %s loaded as structure number %i',fname,snum));
end;
if msg.error,
    add_msg_board(msg.text);
end;
drawnow;

function handles=uncolor(handles,args)

global graph_settings

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: uncolor address');
    return
end;

command=sprintf('uncolor %s',strtrim(args));
undo_cmd='';
[handles,~]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
[message,argout]=set_object(address,'uncolor');
% special handling for selected objects
if strcmp(strtrim(address),'*')
    msg=sprintf('Color of selected objects was restored to previous one.');
    add_msg_board(msg);
    indices=resolve_address('*');
    [m,n]=size(indices);
    for ko=1:m, % loop over all objects
        idepth=length(find(indices(ko,:)>0)); % determine type of current object
        cindices=indices(ko,1:idepth);
        [msg,graphics]=get_object(cindices,'graphics');
        if graphics.mode==1,
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    set(graphics.objects(k),'Color',graph_settings.selected_color);
                end;
            end;
        elseif graphics.mode>1,
            for k=1:length(graphics.objects),
                if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                    set(graphics.objects(k),'FaceColor',graph_settings.selected_color);
                end;
            end;
        end;
    end;
end;

function handles=rotamers(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: rotamers address labelname temperature [no_context]');
    add_msg_board('where labelname is:');
    add_msg_board('IA-PROXYL|MTSL');
    add_msg_board('and temperature is:');
    add_msg_board('cryo|ambient or a number (in K)');
    return
end;

command=sprintf('rotamers %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);
if veto,
    add_msg_board('Rotamers command cancelled by user');
    return
end;

myargs=textscan(args,'%s');
address=char(myargs{1}(1));

if length(myargs{1})>1,
    label=char(myargs{1}(2));
else
    add_msg_board('Assuming default label MTSL');
    label='MTSL';
end;

if length(myargs{1})>2,
    temperature=char(myargs{1}(3));
else
    add_msg_board('Assuming default ambient temperature (298 K)');
    temperature='ambient';
end;

if length(myargs{1})>3,
    no_context=str2double(char(myargs{1}(4)));
else
    no_context= false;
end;

indices=resolve_address(address);

msg=compute_rotamers(indices,label,temperature,no_context);
if msg.error,
    add_msg_board(sprintf('ERROR in rotamer computation: %s',msg.text));
end;

function handles=dihedrals(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: dihedrals adr1 adr2 (location)');
    add_msg_board('where adr1 is:');
    add_msg_board('address of the cystein residue');
    add_msg_board('and adr2 is:');
    add_msg_board('address of the MNT side group');
    return
end;

command=sprintf('dihedrals %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,'noundo');

myargs=textscan(args,'%s');
adr1=char(myargs{1}(1));

if length(myargs{1})>1,
    adr2=char(myargs{1}(2));
else
    add_msg_board('Usage: dihedrals adr1 adr2 (loc)');
    add_msg_board('where adr1 is:');
    add_msg_board('address of the cystein residue');
    add_msg_board('and adr2 is:');
    add_msg_board('address of the MNT side group');
    return
end;

if length(myargs{1})>2,
    loc=char(myargs{1}(3));
else
    loc='';
end;

N=sprintf('%s.N',adr1);
N_coor=get_coor_loc(N,loc);

CA=sprintf('%s.CA',adr1);
CA_coor=get_coor_loc(CA,loc);

CB=sprintf('%s.CB',adr1);
CB_coor=get_coor_loc(CB,loc);

SG=sprintf('%s.SG',adr1);
SG_coor=get_coor_loc(SG,loc);

S1=sprintf('%s.S1',adr2);
S1_coor=get_coor_loc(S1,loc);

C4=sprintf('%s.C4',adr2);
C4_coor=get_coor_loc(C4,loc);

C3=sprintf('%s.C3',adr2);
C3_coor=get_coor_loc(C3,loc);

C5=sprintf('%s.C5',adr2);
C5_coor=get_coor_loc(C5,loc);

CASD=0;
chi1=180*dihedral(N_coor,CA_coor,CB_coor,SG_coor)/pi;
rotvec=chi1;
if ~isempty(S1_coor),
    CASD=norm(CA_coor-S1_coor);
    chi2=180*dihedral(CA_coor,CB_coor,SG_coor,S1_coor)/pi;
    rotvec=[rotvec chi2];
    if ~isempty(C4_coor),
        chi3=180*dihedral(CB_coor,SG_coor,S1_coor,C4_coor)/pi;
        rotvec=[rotvec chi3];
        if ~isempty(C3_coor),
            chi4=180*dihedral(SG_coor,S1_coor,C4_coor,C3_coor)/pi;
            rotvec=[rotvec chi4];
            if ~isempty(C5_coor),
                chi5=180*dihedral(S1_coor,C4_coor,C3_coor,C5_coor)/pi;
                rotvec=[rotvec chi5];
            end;
        end;
    end;
end;

character=rotamer_character(rotvec);
str1=sprintf('Rotamer: %s (%5.1f',character,chi1);
for k=2:length(rotvec),
    str1=sprintf('%s,%5.1f',str1,rotvec(k));
end;
str1=sprintf('%s)',str1);
add_msg_board(str1);
str2=sprintf('C_alpha S_delta distance: %4.1f Å.',CASD);
add_msg_board(str2);
% clipboard('copy',[str1 ';' str2]);
clipboard('copy',sprintf('%s\t%5.1f',character,CASD));
for k=1:length(rotvec),
    fprintf(1,'%5.1f ',rotvec(k));
end;
fprintf(1,'\n');

function coor=get_coor_loc(adr,loc)

if ~isempty(loc),
    adr1=sprintf('%s:%s',adr,loc);
else
    adr1=adr;
end;
[msg,coor]=get_object(adr1,'xyz');
if isempty(coor),
    [msg,coor]=get_object(adr,'xyz');
end;

function handles=label(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: label address labelname temperature');
    add_msg_board('where labelname is:');
    add_msg_board('IA-PROXYL|MTSL');
    add_msg_board('and temperature is:');
    add_msg_board('cryo|ambient or a number (in K)');
    return
end;

command=sprintf('label %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);
if veto
    add_msg_board('Label command cancelled by user');
    return
end

myargs=textscan(args,'%s');
address=char(myargs{1}(1));

if length(myargs{1})>1
    label=char(myargs{1}(2));
else
    add_msg_board('Assuming default label MTSL');
    label='MTSL';
end

if length(myargs{1})>2
    temperature=char(myargs{1}(3));
else
    add_msg_board('Assuming default cryogenic temperature (175 K)');
    temperature='cryo';
end

indices=resolve_address(address);

msg=attach_spin_label(indices,label,temperature);
if msg.error
    add_msg_board(sprintf('ERROR in spin labeling: %s',msg.text));
end

function handles = bilabel(handles,args)

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: bilabel address1 address2 labelname [torsionpot]');
    add_msg_board('where labelname is:');
    add_msg_board('NTA|IDA|IDS');
    return
end

command=sprintf('bilabel %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,'noundo');
if veto
    add_msg_board('Bilabel command cancelled by user');
    return
end

myargs=textscan(args,'%s');
if length(myargs{1}) < 2
    add_msg_board('ERROR: bilabel command requires two residue addresses.');
    return
end

adr1 = char(myargs{1}(1));
adr2 = char(myargs{1}(2));

if length(myargs{1})>2
    label=char(myargs{1}(3));
else
    add_msg_board('Assuming default label IDA');
    label='IDA';
end

if length(myargs{1})>3
    torsionpot = 1000*str2double(char(myargs{1}(4)));
else
    torsionpot = [];
end

[ind1,msg] = resolve_address(adr1);

if msg.error
    add_msg_board(sprintf('ERROR: First address %s could not be resolved (%s)',adr1,msg.text));
    return
end
if length(ind1) ~= 4
    add_msg_board(sprintf('ERROR: First address %s does not specify a residue',adr1));
    return
end

[ind2,msg] = resolve_address(adr2);

if msg.error
    add_msg_board(sprintf('ERROR: Second address %s could not be resolved (%s)',adr2,msg.text));
    return
end
if length(ind2) ~= 4
    add_msg_board(sprintf('ERROR: Second address %s does not specify a residue',adr2));
    return
end

if ind1(1) ~= ind2(1)
    add_msg_board(sprintf('ERROR: Addresses %s and %s belong to different structures',adr1,adr2));
    return
end

if ind1(2) ~= ind2(2)
    add_msg_board(sprintf('ERROR: The two residues (%s and %s) must belong to the same chain',adr1,adr2));
    return
end

[msg,coor1] = get_object(ind1(1:3),'xyz_heavy');
if msg.error
    add_msg_board(sprintf('ERROR: Coordinates of the labelled chain could not be obtained (%s)',msg.text));
    return
end
[msg,elements] = get_object(ind1(1:3),'elements_heavy');
if msg.error
    add_msg_board(sprintf('ERROR: Element specifiers of the labelled chain could not be obtained (%s)',msg.text));
    return
end
ecoor = [elements,coor1];

two_pronged_label(adr1,adr2,label,ecoor,torsionpot);

function handles = blscan(handles,args)

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: blscan chain_address labelname outputfile [torsionpot]');
    add_msg_board('where labelname is:');
    add_msg_board('NTA|IDA');
    return
end;

command=sprintf('blscan %s',strtrim(args));
handles = cmd_history(handles,command,'noundo');

myargs=textscan(args,'%s');
if length(myargs{1}) < 1
    add_msg_board('ERROR: blscan command requires chain (model) address.');
    return
end

adr1 = char(myargs{1}(1));

if length(myargs{1})>1
    label=char(myargs{1}(2));
else
    add_msg_board('Assuming default label IDA');
    label='IDA';
end

[ind1,msg] = resolve_address(adr1);

if msg.error
    add_msg_board(sprintf('ERROR: Chain address %s could not be resolved (%s)',adr1,msg.text));
    return
end
if length(ind1) < 2 || length(ind1) > 3
    add_msg_board(sprintf('ERROR: Address %s does not specify a chain or chain model',adr1));
    return
end

if length(myargs{1})>2
    fname = char(myargs{1}(3));
else
    add_msg_board('Assuming default output file name bilabel_site_scan.txt');
    fname = 'bilabel_site_scan';
end

if length(myargs{1})>3
    if strcmpi(char(myargs{1}(4)),'native') || strcmpi(char(myargs{1}(4)),'true')
        native = 1;
    elseif str2double(char(myargs{1}(4))) == 1
        native = 1;
    else
        native = 0;
    end
else
    native = 0;
end

if length(myargs{1})>4
    torsionpot = 1000*str2double(char(myargs{1}(5)));
else
    torsionpot = [];
end

msg = bilabel_site_scan(ind1,label,fname,native,torsionpot);

if msg.error
    add_msg_board(sprintf('ERROR: Bilabel site scan failed (%s)',msg.text)); 
end


function handles=libtest(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: libtest testdata outfile library [no_context]');
    return
end;

command=sprintf('libtest %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);
if veto,
    add_msg_board('Libtest command cancelled by user');
    return
end;

myargs=textscan(args,'%s');
fname=char(myargs{1}(1));
if length(myargs{1})>1,
    fname2=char(myargs{1}(2));
else
    [path,name,ext] = fileparts(fname);
    name=strcat(name,'_prc',ext);
    fname2=fullfile(path,name);
end;

if length(myargs{1})>2,
    library=char(myargs{1}(3));
else
    library='';
end;

if length(myargs{1})>3,
    no_context=str2double(char(myargs{1}(4)));
else
    no_context = false;
end;

testset = rd_rotlib_test_data(fname);

if isempty(testset),
    add_msg_board('ERROR: Empty test data set.');
else
    testset=update_libtest(testset,fname2,handles,library,no_context);
end;
if ~isempty(testset),
    wr_rotlib_test_data(fname2,testset);
end;

function handles = libcomp(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: libcomp address lib1 lib2 >outfile');
    add_msg_board('   or: libcomp <infile >outfile');
    add_msg_board('If outfile is missing, output goes to message board.');
    return
end;

command=sprintf('libcomp %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);
if veto,
    add_msg_board('Libcomp command cancelled by user');
    return
end;

myargs=textscan(args,'%s');
outfile = '';
lib1 = '';
lib2 = '';
adr = '';
arg1 = char(myargs{1}(1));
if arg1(1) == '<',
    infile = arg1(2:end);
else
    adr = arg1;
    infile = '';
end;
if isempty(infile),
    if length(myargs{1})<3,
        add_msg_board('ERROR: Not enough arguments.');
        add_msg_board('Usage: libcomp address lib1 lib2 >outfile');
        add_msg_board('   or: libcomp <infile >outfile');
        add_msg_board('If outfile is missing, output goes to message board.');
        return
    end;   
    lib1 = char(myargs{1}(2));
    lib2 = char(myargs{1}(3));
    if length(myargs{1}) > 3,
        arg4 = char(myargs{1}(4));
        if arg4(1) == '>',
            outfile = arg4(2:end);
        else
            add_msg_board('ERROR: Fourth argument must begin with > or must not be present.');
            add_msg_board('Usage: libcomp address lib1 lib2 >outfile');
            add_msg_board('   or: libcomp <infile >outfile');
            add_msg_board('If outfile is missing, output goes to message board.');
            return
        end;
    end;
else
    if length(myargs{1}) > 1,
        arg2 = char(myargs{1}(2));
        if arg2(1) == '>',
            outfile = arg2(2:end);
        else
            add_msg_board('ERROR: Second argument must begin with > or must not be present.');
            add_msg_board('Usage: libcomp address lib1 lib2 >outfile');
            add_msg_board('   or: libcomp <infile >outfile');
            add_msg_board('If outfile is missing, output goes to message board.');
            return
        end;
    end;
end;

[err,msg] = do_libcomp(adr,lib1,lib2,infile,outfile);

if err,
    add_msg_board(sprintf('ERROR: %s',msg));
end;

function handles=colorscheme(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: colorscheme address scheme_type (address2) (sensitivity)');
    add_msg_board('where scheme_type is:');
    add_msg_board('secondary|chain|sequence|Bfactor|Bfactor_tight|charge|ensmble|difference');
    add_msg_board('and address2 is required only for scheme_type difference.');
    return
end;

command=sprintf('colorscheme %s',strtrim(args));
undo_cmd=sprintf('uncolor %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);
scheme='';

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    scheme=char(myargs{1}(2));
end;
if length(myargs{1})>2,
    arg3=char(myargs{1}(3));
    args={scheme arg3};
else
    args={scheme};
end;
if length(myargs{1})>3,
    arg4=char(myargs{1}(4));
    args={scheme arg3 arg4};
end;

if ~isempty(scheme),
    if strcmp(scheme,'difference') && length(myargs{1})<3,
        add_msg_board('ERROR: Color scheme "difference" requires a second address argument.');
    else
        [message,argout]=set_object(address,'colorscheme',args);
        if message.error,
            add_msg_board(sprintf('WARNING: %s',message.text));
        end;
        % special handling for selected objects
        if strcmp(strtrim(address),'*')
            msg=sprintf('Color of selected objects was changed to scheme "%s".',scheme);
            add_msg_board(msg);
            highlight_selection
        end;
    end;
else
    add_msg_board('ERROR: No color scheme specified.');
end;

function handles=repack_structure(handles,args)

global model

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: repack adr');
    add_msg_board('where adr is the address of the structure whose sidechains are to be repacked.');
    return
end;

command=sprintf('repack %s',strtrim(args));
undo_cmd='  ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=strtrim(char(myargs{1}(1)));

option='';
if length(myargs{1})>1,
    option=char(myargs{1}(2));
end;

crystal=false;
switch option
    case 'crystal'
        crystal=true;
end;

[indices,msg]=resolve_address(address);

if msg.error,
    add_msg_board('ERROR: Structure address could not be resolved.');
    add_msg_board(msg.text);
    return
end;

if length(indices)>1,
    add_msg_board('ERROR: Lower level object addressed. No sidechain repacking.');
end;

snum=indices(1);
models=length(model.structures{snum}(1).residues);
for modnum=1:models, % repack all models of current structure
    snum1=repack(snum,modnum,crystal);
    if isempty(snum1),
        add_msg_board(sprintf('ERROR: Model %i in structure %s could not be repacked.',modnum,address));
    end;
end;

function handles=replace_restype(handles,args)

global model

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: replace adr restype');
    add_msg_board('where adr is the adress of the structure, restype is the PDB code or pseudocode of non-native residues to be replaced by their native equivalents.');
    return
end;

command=sprintf('replace %s',strtrim(args));
undo_cmd='  ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
adr=strtrim(char(myargs{1}(1)));
replacements=[':' strtrim(char(myargs{1}(2))) ':'];

[indices,msg]=resolve_address(adr);

if msg.error,
    add_msg_board('ERROR: Structure address could not be resolved.');
    add_msg_board(msg.text);
    return
end;

selected=false;
if length(indices)>1,
    add_msg_board('Warning: Lower level object addressed. Only selected residues are replaced.');
    selected=true;
end;

[repnum,msg]=replace(indices,replacements,selected);

if ~msg.error,
    add_msg_board(sprintf('%i residues were replaced in structure %s.',repnum,adr));
else
    add_msg_board('ERROR: No residue replacement.');
    add_msg_board(msg.text);
end;

function handles=echo(handles,args)

command=sprintf('echo');
undo_cmd=sprintf('echo');
[handles,veto]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: echo text');
    return
else
    add_msg_board(strtrim(args));
end;


function handles=copy_3D(handles)

global hMain
global hModel

if ~hMain.detached,
    add_msg_board('ERROR: Visualization can be copied to clipboard only when model window is detached.');
    add_msg_board('Use command "detach" first.');
    return
end;

add_msg_board('Copying. Please be patient...');
if hMain.atom_graphics_auto,
    switch_it=true;
    adjust_atom_graphics(false);
else
    switch_it=false;
end;
print(hModel.figure,'-dbitmap');
if switch_it,
    adjust_atom_graphics(true);
end;
add_msg_board('Copying finished. Display will update shortly.');

function handles=lock(handles,args)

global hMain
global model

command=sprintf('lock');
undo_cmd=sprintf('unlock');
[handles,veto]=cmd_history(handles,command,undo_cmd);

set(handles.uitoggletool_lock,'State','on');
model.locked=1;
if isfield(hMain,'hierarchy_window') && ishandle(hMain.hierarchy_window),
    hhandles=guidata(hMain.hierarchy_window);
    set(hhandles.pushbutton_loop,'Enable','off');
    set(hhandles.pushbutton_helix,'Enable','off');
    set(hhandles.pushbutton_strand,'Enable','off');
end;
if isfield(hMain,'hierarchy_window_large') && ishandle(hMain.hierarchy_window_large),
    hhandles=guidata(hMain.hierarchy_window_large);
    set(hhandles.pushbutton_loop,'Enable','off');
    set(hhandles.pushbutton_helix,'Enable','off');
    set(hhandles.pushbutton_strand,'Enable','off');
end;

function handles=unlock(handles,args)

global hMain
global model

command=sprintf('unlock');
undo_cmd=sprintf('lock');
[handles,veto]=cmd_history(handles,command,undo_cmd);

set(handles.uitoggletool_lock,'State','off');
model.locked=0;
if isfield(hMain,'hierarchy_window') && ishandle(hMain.hierarchy_window),
    hhandles=guidata(hMain.hierarchy_window);
    set(hhandles.pushbutton_loop,'Enable','on');
    set(hhandles.pushbutton_helix,'Enable','on');
    set(hhandles.pushbutton_strand,'Enable','on');
end;
if isfield(hMain,'hierarchy_window_large') && ishandle(hMain.hierarchy_window_large),
    hhandles=guidata(hMain.hierarchy_window_large);
    set(hhandles.pushbutton_loop,'Enable','on');
    set(hhandles.pushbutton_helix,'Enable','on');
    set(hhandles.pushbutton_strand,'Enable','on');
end;

function handles=dssp_assign(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: dssp address');
    return
end;

command=sprintf('dssp %s',strtrim(args));
undo_cmd='  ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
indices=resolve_address(address);
reassign_secondary(indices);

function handles=transparency(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: transparency address falpha');
    add_msg_board('where falpha is between 0 (transparent) and 1 (opaque).');
    return
end;

command=sprintf('transparency %s',strtrim(args));
undo_cmd=sprintf('untransparency %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

if veto
    return
end;

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    alpha=str2double(char(myargs{1}(2)));
else
    add_msg_board('Transparency value missing.');
    return
end;
if isnan(alpha) || alpha>1
    add_msg_board('Transparency value must be a number not greater than 1.');
    return
end;
if strcmp(address(1),'$'), % surface request
    msg=set_surface(address,'transparency',alpha);
    if msg.error,
        if iscell(msg.text)
            for k=1:length(msg.text),
                add_msg_board(msg.text{k});
            end,
        else
            add_msg_board(msg.text);
        end;
    end;
else
    [message,argout]=set_object(address,'transparency',{alpha});
end;

function handles=untransparency(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: untransparency address');
    return
end;

command=sprintf('untransparency %s',strtrim(args));
undo_cmd='';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
[message,argout]=set_object(address,'untransparency');

function handles=SAS(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: SAS address tag radius');
    add_msg_board('where address must address a single structure, chain, coordinate set, or residue');
    add_msg_board('and tag is a tag for future access of the surface');
    add_msg_board('radius is an optional argument for the probe radius, it defaults to 1.5 Å');
    return
end;

command=sprintf('SAS %s',strtrim(args));
undo_cmd='  ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=strtrim(char(myargs{1}(1)));
if length(myargs{1})>1,
    tag=strtrim(char(myargs{1}(2)));
else
    add_msg_board('Usage: SAS address tag radius');
    add_msg_board('where address must address a single structure, chain, coordinate set, or residue');
    add_msg_board('and tag is a tag for future access of the surface');
    add_msg_board('radius is an optional argument for the probe radius, it defaults to 1.5 Å');
    return
end;

if length(myargs{1})>2,
    radius=str2double(char(myargs{1}(3)));
    if isnan(radius),
        add_msg_board('Warning: Wrong probe radius argument. Reverting to default radius 1.5 Å');
    end;
else
    radius=1.5;
end;

msg=mk_SAS(address,tag,radius);

if msg.error,
    add_msg_board(msg.text);
end;

function handles=scopy(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: scopy adr newtag');
    add_msg_board('where adr is the address of the structure to be copied');
    add_msg_board('and newtag is the new structure tag');
    return
end;

command=sprintf('scopy %s',strtrim(args));
undo_cmd='  ';
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=strtrim(char(myargs{1}(1)));
if length(myargs{1})>1,
    idCode=strtrim(char(myargs{1}(2)));
else
    add_msg_board('Usage: scopy adr newtag');
    add_msg_board('where adr is the address of the structure to be copied');
    add_msg_board('and newtag is the new structure tag');
    return
end;

[indices,msg]=resolve_address(address);

if msg.error,
    add_msg_board('ERROR: Structure address could not be resolved.');
    add_msg_board(msg.text);
    return
end;

if length(indices)>1,
    add_msg_board('Warning: Lower level object addressed, but whole structure will be copied.');
end;

[snum1,msg]=copy_structure(indices(1),idCode);

if msg.error,
    add_msg_board(msg.text);
end;

function handles=detach(handles)
global hMain
global hModel

if ~hMain.detached;
    handles.oldpos=get(handles.axes_model,'Position');
    hMain.oldpos=handles.oldpos;
    hMain.panel_model=handles.panel_model;
    hMain.axes_model=handles.axes_model;
    model_window;
    newpos=get(hModel.figure,'Position');
    switch hMain.color
        case 'black'
            set(hModel.figure,'Color','k');
        case 'grey'
            set(hModel.figure,'Color',[0.941,0.941,0.941]);
        case 'white'
            set(hModel.figure,'Color','w');
    end;
    rot_state=get(handles.uitoggletool_rotate,'State');
    if strcmp(rot_state,'on'),
        view3D(hMain.figure,'off');
        view3D(hModel.figure,'rot');
    end;
    zoom_state=get(handles.uitoggletool_zoom,'State');
    if strcmp(zoom_state,'on'),
        view3D(hMain.figure,'off');
        view3D(hModel.figure,'zoom');
    end;
    pan_state=get(handles.uitoggletool_pan,'State');
    if strcmp(pan_state,'on'),
        view3D(hMain.figure,'off');
        view3D(hModel.figure,'pan');
    end;
    newpos(1)=0.01;
    newpos(2)=0.01;
    newpos(3)=0.98;
    newpos(4)=0.98;
    set(handles.axes_model,'Parent',hModel.figure);
    set(handles.axes_model,'Position',newpos);
    set(handles.panel_model,'Title','Model Detached');
    hMain.detached=1;
    command=sprintf('detach');
    undo_cmd=sprintf('attach');
    [handles,veto]=cmd_history(handles,command,undo_cmd);
end;

function handles=hide(handles,args)

global model

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: hide address');
    return
end;

mode='';

command=sprintf('hide %s',strtrim(args));
undo_cmd=sprintf('show %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
axes(handles.axes_model);

if ~isempty(strfind(address,'|')) % bilabel case
    sep = strfind(address,'|');
    adr1 = address(1:sep-1);
    adr2 = address(sep+1:end);
    [ind1,msg] = resolve_address(adr1);
    if msg.error
        add_msg_board(sprintf('ERROR: First residue address %s could not be resolved (%s). Aborting.',adr1,msg.error));
        return
    end
    if length(ind1) ~= 4
        add_msg_board(sprintf('ERROR: First address %s does not specify a residue. Aborting.',adr1));
        return
    end
    [ind2,msg] = resolve_address(adr2);
    if msg.error
        add_msg_board(sprintf('ERROR: Second residue address %s could not be resolved (%s). Aborting.',adr2,msg.error));
        return
    end
    if length(ind2) ~= 4
        add_msg_board(sprintf('ERROR: Second address %s does not specify a residue. Aborting.',adr2));
        return
    end
    if ~isfield(model,'bisites')
        add_msg_board('ERROR: Current model does not contain labels attached to two residues. Aborting.');
        return
    else
        check1 = [ind1;ind2];
        check2 = [ind2;ind1];
        for scan = 1:length(model.bisites)
            for site = 1:length(model.bisites{scan}.sites)
                if sum(sum(abs(check1-model.bisites{scan}.sites{site}.indices))) == 0 ...
                        || sum(sum(abs(check2-model.bisites{scan}.sites{site}.indices))) == 0
                    objects = model.bisites{scan}.sites{site}.gobjects;
                    for k = 1:length(objects)
                        delete(objects(k));
                    end
                end
            end
        end
    end
end

if strcmp(address(1),'$'), % surface request
    msg=set_surface(address,'visible','off');
    if msg.error,
        if iscell(msg.text)
            for k=1:length(msg.text),
                add_msg_board(msg.text{k});
            end,
        else
            add_msg_board(msg.text);
        end;
    end;
else
    [message,argout]=set_object(address,'hide');
end;
lighting gouraud

function handles=delete_cmd(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: delete address');
    return
end;

mode='';

command=sprintf('delete %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
axes(handles.axes_model);
if strcmp(address(1),'$'), % surface request
    msg=set_surface(address,'delete');
    if msg.error,
        if iscell(msg.text)
            for k=1:length(msg.text),
                add_msg_board(msg.text{k});
            end,
        else
            add_msg_board(msg.text);
        end;
    end;
else
    add_msg_board('Warning: Current version can delete only surfaces');
end;
lighting gouraud

function handles=attach(handles)
global hMain
global hModel

if hMain.detached,
    set(handles.axes_model,'Position',hMain.oldpos);
    set(handles.axes_model,'Parent',hMain.panel_model);
    set(handles.panel_model,'Title','Model');
    view3D(hModel.figure,'off');
    close(hModel.figure);
    rot_state=get(handles.uitoggletool_rotate,'State');
    if strcmp(rot_state,'on'),
        view3D(hMain.figure,'rot');
    end;
    zoom_state=get(handles.uitoggletool_zoom,'State');
    if strcmp(zoom_state,'on'),
        view3D(hMain.figure,'zoom');
    end;
    pan_state=get(handles.uitoggletool_pan,'State');
    if strcmp(pan_state,'on'),
        view3D(hMain.figure,'pan');
    end;
    hMain.detached=0;
    command=sprintf('attach');
    undo_cmd=sprintf('detach');
    [handles,veto]=cmd_history(handles,command,undo_cmd);
end;

function handles=orthographic(handles)

set(handles.axes_model,'Projection','orthographic');
command=sprintf('ortho');
undo_cmd=sprintf('persp');
[handles,veto]=cmd_history(handles,command,undo_cmd);

function handles=perspective(handles)

set(handles.axes_model,'Projection','perspective');
command=sprintf('persp');
undo_cmd=sprintf('ortho');
[handles,veto]=cmd_history(handles,command,undo_cmd);

function handles=EPRdyn(handles,args)
% serves as an example of a non-undoable command
command=sprintf('EPRdyn %s',args);
[handles,veto]=cmd_history(handles,command);
if veto,
    msg=sprintf('Command was vetoed by user.');
else
    msg=sprintf('Command should be executed.');
end;
add_msg_board(msg);

function handles=plot_cmd(handles,args)

global graph_settings

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: plot adr1 adr2 [[width] color]');
    add_msg_board('Only atom or location addresses are allowed.');
    return
end

command = sprintf('plot %s',args);
[handles,veto]=cmd_history(handles,command);
if veto
    msg=sprintf('Command was vetoed by user.');
else
    msg=sprintf('Command should be executed.');
end

myargs=textscan(args,'%s');
adr1=char(myargs{1}(1));
adr2=char(myargs{1}(2));
[msg1,coor1] = get_object(adr1,'coor');
[msg2,coor2] = get_object(adr2,'coor');

if msg1.error
    add_msg_board(msg1.text);
    return;
end
if msg2.error
    add_msg_board(msg1.text);
    return;
end

if length(myargs{1})>2
    lwidth=str2double(char(myargs{1}(3)));
else
    add_msg_board('Linewidth argument missing. Default is 6.');
    lwidth = 6;
end

if length(myargs{1})>3
    if length(myargs{1})>4
        rgb = zeros(1,3);
        for k = 1:3
            rgb(k) = str2double(char(myargs{1}(3+k)));
        end
    else
        colname=char(myargs{1}(4));
        colind=tag2id(colname,graph_settings.color_tags);
        if ~isempty(colind)
            rgb=graph_settings.colors(colind,:);
        else
            add_msg_board('Warning: This SVG color does not exist. Defaulting to darkgreen.');
            rgb=[0,0.39,0];
        end
    end
else
    add_msg_board('Color argument missing. Default is darkgreen.');
    rgb = [0,0.39,0];
end

linestyle ='-';
if length(myargs{1})>6
    linestyle = char(myargs{1}(7));
end

axes(handles.axes_model);
plot3([coor1(1),coor2(1)],[coor1(2),coor2(2)],[coor1(3),coor2(3)],'LineStyle',linestyle,'LineWidth',lwidth,'Color',rgb);


function handles=show(handles,args)

global model
global graph_settings

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: show address mode');
    add_msg_board('Type "help show" for on-line help on modes.');
    return
end;

mode='';

command=sprintf('show %s',strtrim(args));
undo_cmd=sprintf('hide %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    mode=char(myargs{1}(2));
else
    add_msg_board('Warning: No mode specified.');
    add_msg_board('Default mode may apply or nothing will be displayed');
end;
axes(handles.axes_model);

if exist('mode','var') && ~isempty(strfind(address,'|')) && strcmp(mode,'label') % bilabel case
    sep = strfind(address,'|');
    adr1 = address(1:sep-1);
    adr2 = address(sep+1:end);
    [ind1,msg] = resolve_address(adr1);
    if msg.error
        add_msg_board(sprintf('ERROR: First residue address %s could not be resolved (%s). Aborting.',adr1,msg.error));
        return
    end
    if length(ind1) ~= 4
        add_msg_board(sprintf('ERROR: First address %s does not specify a residue. Aborting.',adr1));
        return
    end
    [ind2,msg] = resolve_address(adr2);
    if msg.error
        add_msg_board(sprintf('ERROR: Second residue address %s could not be resolved (%s). Aborting.',adr2,msg.error));
        return
    end
    if length(ind2) ~= 4
        add_msg_board(sprintf('ERROR: Second address %s does not specify a residue. Aborting.',adr2));
        return
    end
    rgb = [];
    if length(myargs{1})>2
        colstr = char(myargs{1}(3));
        if length(myargs{1}) >= 5
            rgb = [str2double(myargs{1}(3)),str2double(myargs{1}(4)),str2double(myargs{1}(5))];
            if max(rgb) > 1
                rgb = rgb/255;
            end
        else
            colind=tag2id(colstr,graph_settings.color_tags);
            if ~isempty(colind)
                rgb = graph_settings.colors(colind,:);
            else
                add_msg_board('Warning: Invalid color specification (using default color)');
                rgb = [];
            end
        end
    end
    if ~isfield(model,'bisites')
        add_msg_board('ERROR: Current model does not contain labels attached to two residues. Aborting.');
        return
    else
        check1 = [ind1;ind2];
        check2 = [ind2;ind1];
        for scan = 1:length(model.bisites)
            for site = 1:length(model.bisites{scan}.sites)
                if sum(sum(abs(check1-model.bisites{scan}.sites{site}.indices))) == 0 ...
                        || sum(sum(abs(check2-model.bisites{scan}.sites{site}.indices))) == 0
                    model.bisites{scan}.sites{site}.gobjects = plot_bilabel(model.bisites{scan}.sites{site}.coor(:,1:4),ind1(1:2),rgb);
                end
            end
        end
    end
    lighting gouraud
    return
end

if strcmp(address(1),'$'), % surface request
    msg=set_surface(address,'visible','on');
    if msg.error,
        if iscell(msg.text)
            for k=1:length(msg.text),
                add_msg_board(msg.text{k});
            end,
        else
            add_msg_board(msg.text);
        end;
    end;
elseif length(myargs{1}) <= 2
    [message,argout]=set_object(address,'show',{mode});
else
    argin{1} = mode;
    argin{2} = char(myargs{1}(3));
    if length(myargs{1}) > 3
        argin{3} = char(myargs{1}(4));
    end;
    [message,argout]=set_object(address,'show',argin);
end;

lighting gouraud

if strcmp(strtrim(address),'*')
    highlight_selection
end;

function handles=motion(handles,args)

global graph_settings
global model
global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: motion move_adr target_adr');
    add_msg_board('Type "help motion" for on-line help.');
    return
end;

command=sprintf('motion %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);

if isempty(args),
    add_msg_board('ERROR: Address arguments missing.');
end;
myargs=textscan(args,'%s');
adr1=char(myargs{1}(1));
if length(myargs{1})>1,
    adr2=char(myargs{1}(2));
else
    add_msg_board('ERROR: Target structure address missing.');
end;
axes(handles.axes_model);
[message,argout]=set_object('[:]','hide');
[transmat,target,moving,alpha]=opt_superposition_ensemble(adr1,adr2);
if ~isempty(target) && ~isempty(moving) && ~isempty(alpha),
    if isfield(model,'motion'),
        k=length(model.motion)+1;
    else
        k=1;
    end;
    [m,n]=size(target);
    id_target=resolve_address(adr1);
    id_moving=resolve_address(adr2);
    if isempty(id_target),
        add_msg_board('ERROR: Address of target structure is wrong.');
        return;
    else
        id_target=id_target(1);
    end;
    if isempty(id_moving),
        add_msg_board('ERROR: Address of moving structure is wrong.');
        return;
    else
        id_moving=id_moving(1);
    end;
    tag1=id2tag(id_target,model.structure_tags);
    tag2=id2tag(id_moving,model.structure_tags);
    tag=['motion:' tag1 '_to_' tag2];
    model.motion(k).tag=tag;
    model.motion(k).color=graph_settings.motion_arrow;
    model.motion(k).transparency=alpha;
    model.motion(k).radius=0.4;
    model.motion(k).template=moving;
    model.motion(k).target=target;
    model.motion(k).active=1;
    model.motion(k).snum=id_target;
    model.motion(k).transform=0;
    model.motion(k).gobjects=[];
    
    transform_structure(id_moving,transmat);
    [message,argout]=set_object(adr1,'show',{'string'});
    [message,argout]=set_object(adr1,'color',{graph_settings.decent_string});
     msg=set_surface(['$' tag],'show');
else
    add_msg_board('ERROR: Structure superposition failed.');
end;

camlookat(hMain.axes_model);
lighting gouraud

function handles = new_model(handles,args)

global model
global hMain
global MMM_info

old_undo = hMain.store_undo;
if ~isempty(args)
    myargs=textscan(args,'%s');
    if ~isempty(myargs{1}) && strcmp(char(myargs{1}),'!')
        hMain.store_undo = false;
    end
end

command=sprintf('new %s',strtrim(args));
[handles,veto] = cmd_history(handles,command);
hMain.store_undo = old_undo;

if veto
    add_msg_board('Command "new" cancelled by user.');
    return
end

if hMain.detached
    set(handles.axes_model,'Position',hMain.oldpos);
    set(handles.axes_model,'Parent',hMain.panel_model);
    set(handles.panel_model,'Title','Model');
    view3D(hModel.figure,'off');
    close(hModel.figure);
    rot_state=get(handles.uitoggletool_rotate,'State');
    if strcmp(rot_state,'on')
        view3D(hMain.figure,'rot');
    end
    zoom_state=get(handles.uitoggletool_zoom,'State');
    if strcmp(zoom_state,'on')
        view3D(hMain.figure,'zoom');
    end
    pan_state=get(handles.uitoggletool_pan,'State');
    if strcmp(pan_state,'on')
        view3D(hMain.figure,'pan');
    end
    hMain.detached=0;
end


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
hMain.virgin=0;

handles.MMM.Name = MMM_info.title;

clear global model
model = [];

function handles=mk_synonym(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: synonym address new_name');
    add_msg_board('where address must specify a structure or chain.');
    return
end;

newtag='';

command=sprintf('synonym %s',strtrim(args));
undo_cmd=sprintf('unsynonym %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
oldtag=char(myargs{1}(1));
if length(myargs{1})>1,
    newtag=char(myargs{1}(2));
else
    add_msg_board('ERROR: New name is missing.');
    return
end;
synonym(oldtag,newtag);

function handles=def_helix(handles,args)

global model
global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: helix chain_address name start_residue end_residue [TM]');
    return
end;

if model.locked,
    add_msg_board('ERROR: Must unlock secondary structure before redefinition');
    return;
end;

newtag='';

command=sprintf('helix %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,'null');

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    newtag=char(myargs{1}(2));
else
    add_msg_board('ERROR: Name of helix is missing');  
    return
end;
if length(myargs{1})>2,
    sstart=str2double(char(myargs{1}(3)));
else
    add_msg_board('ERROR: Start of residue range is missing');  
    return
end;
if length(myargs{1})>3,
    send=str2double(char(myargs{1}(4)));
else
    add_msg_board('ERROR: End of residue range is missing');  
    return
end;
type='';
if length(myargs{1})>4,
    type=strtrim(char(myargs{1}(5)));
end;
if strcmpi(type,'TM') || strcmpi(type,'transmembrane'),
    TM=1;
else
    TM=0;
end;
if hMain.hierarchy_display,
%     if hMain.large,
%         close(hMain.hierarchy_window_large);
%     else
%         close(hMain.hierarchy_window);
%     end;
    close(hMain.hierarchy_window);
    hMain.hierarchy_display=0;
end;    
secondary(address,'helix',newtag,[sstart,send],TM);

function handles=def_sheet(handles,args)

global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: sheet chain_address name start_residue end_residue [TM]');
    return
end;

newtag='';

command=sprintf('sheet %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,'null');

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    newtag=char(myargs{1}(2));
else
    add_msg_board('ERROR: Name of sheet is missing');  
    return
end;
if length(myargs{1})>2,
    sstart=str2double(char(myargs{1}(3)));
else
    add_msg_board('ERROR: Start of residue range is missing');  
    return
end;
if length(myargs{1})>3,
    send=str2double(char(myargs{1}(4)));
else
    add_msg_board('ERROR: End of residue range is missing');  
    return
end;
type='';
if length(myargs{1})>4,
    type=char(myargs{1}(5));
end;
if strcmpi(type,'TM') || strcmpi(type,'transmembrane'),
    TM=1;
else
    TM=0;
end;
if hMain.hierarchy_display,
%     if hMain.large,
%         close(hMain.hierarchy_window_large);
%     else
%         close(hMain.hierarchy_window);
%     end;
    close(hMain.hierarchy_window);
    hMain.hierarchy_display=0;
end;    
secondary(address,'sheet',newtag,[sstart,send],TM);

function handles=def_loop(handles,args)

global model
global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: loop chain_address name start_residue end_residue');
    return
end;

newtag='';

if model.locked,
    add_msg_board('ERROR: Must unlock secondary structure before redefinition');
    return;
end;

command=sprintf('loop %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,'null');

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    newtag=char(myargs{1}(2));
else
    add_msg_board('ERROR: Name of loop is missing');  
    return
end;
if length(myargs{1})>2,
    sstart=str2double(char(myargs{1}(3)));
else
    add_msg_board('ERROR: Start of residue range is missing');  
    return
end;
if length(myargs{1})>3,
    send=str2double(char(myargs{1}(4)));
else
    add_msg_board('ERROR: End of residue range is missing');  
    return
end;
if hMain.hierarchy_display,
%     if hMain.large,
%         close(hMain.hierarchy_window_large);
%     else
%         close(hMain.hierarchy_window);
%     end;
    close(hMain.hierarchy_window);
    hMain.hierarchy_display=0;
end;    
secondary(address,'loop',newtag,[sstart,send]);

function handles=symmetry(handles)

global hMain

silent=~hMain.store_undo;

command=sprintf('symmetry');
[handles,veto]=cmd_history(handles,command);
if veto,
    add_msg_board('Transformation to symmetry frame cancelled by user.');
else
    [msg,snum]=symmetry_transform;
    if msg.error,
        add_msg_board(msg.text);
    end;
    if ~msg.error || msg.error>100,
        if ~silent && ~isempty(snum),
            [message,argout]=set_object(snum,'hide');
        end;
    end;
end;

function handles=def_domain(handles,args)

global hMain

newtag='';

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: domain (+) name address_of_residues_to_add');
    return
end;

command=sprintf('domain %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,'null');

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
extend=0;
if length(myargs{1})>1,
    if strcmp(myargs{1}(1),'+')
        address=char(myargs{1}(2));
        extend=1;
    end;
else
    if strcmp(myargs{1}(1),'+')
        add_msg_board('ERROR: Residue address missing.');
    else
        add_msg_board('ERROR: Domain name is missing.');
    end;
    return
end;
tag_pos=2+extend;
if length(myargs{1})>=tag_pos,
    newtag=char(myargs{1}(tag_pos));
else
    add_msg_board('ERROR: Domain name is missing');  
    return
end;
create=2;
msg=domain(newtag,create,address);
create=1;
if extend,
    create=0;
    if msg.error~=8,
        add_msg_board('Warning: Domain does not yet exist and will be created.');
    end;
elseif hMain.store_undo,
    if msg.error==8,
        add_msg_board(msg.text);
        answer=questdlg('Do you want to redefine domain?','Domain already exists','Redefine','Cancel','Cancel');
        switch answer
            case 'Redefine'
                create=1;
            case 'Cancel'
                add_msg_board('Use syntax "domain + address name" to add to an existing domain.');
                return
        end;
    end;
end;
msg=domain(newtag,create,address);
if msg.error,
    add_msg_board(sprintf('ERROR: %s',msg.text));
end;

function handles=select(handles,args)
    
global model
global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: select address');
    return
end;

command=sprintf('select %s',strtrim(args));
undo_cmd=sprintf('unselect %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
indices=resolve_address(address);

if ~isfield(model,'selections') || ~isfield(model,'selected');
    model.selections=0;
    model.selected={};
end;

[m,n]=size(indices);
for k=1:m,
    idepth=length(find(indices(k,:)>0)); % determine type of current object
    cindices=indices(k,1:idepth);
    sec_indices=resolve_secondary_address(address,indices(k,:));
    if ~isempty(sec_indices),
        snum=indices(k,1);
        cnum=indices(k,2);
        switch sec_indices(1)
            case 0
                model.structures{snum}(cnum).loop_defs{sec_indices(2)}.selected=1;
            case 1
                model.structures{snum}(cnum).helix_defs{sec_indices(2)}.selected=1;
            case 2
                model.structures{snum}(cnum).sheet_defs{sec_indices(2)}.selected=1;
        end;
    end;
    if model.selections==0,
        model.selections=1;
        model.selected{1}=cindices;
    else
        store=1;
        for kk=1:model.selections,
            sindices=model.selected{kk};
            sdepth=length(find(sindices>0)); % determine type of current object
            sindices=sindices(1:sdepth);
            if sdepth==idepth
                if cindices==sindices
                    store=0;
                    add_msg_board('Double selection of same object ignored');
                end;
            end;
        end;
        if store
            model.selections=model.selections+1;
            model.selected{model.selections}=cindices;
        end;
    end;
end;
highlight_selection;
% Display of current selection

if model.selections>=1,
    add_msg_board('--- Selected ---');
    for k=1:model.selections,
        % [message,prev_info]=get_object(model.selected{k},'info');
        % add_msg_board(prev_info{1});
        adr=mk_address(model.selected{k},1);
        add_msg_board(adr);
    end;
else
    add_msg_board('--- Nothing selected ---');
end;

if length(myargs{1})<=1,
    if hMain.hierarchy_display,
%         if hMain.large,
%             hwhandle=hMain.hierarchy_window_large;
%             hhandles=guidata(hwhandle);
%             hierarchy_window_large('sequence_display',hhandles);
%         else
%             hwhandle=hMain.hierarchy_window;
%             hhandles=guidata(hwhandle);
%             hierarchy_window('sequence_display',hhandles);
%         end;
        hwhandle=hMain.hierarchy_window;
        hhandles=guidata(hwhandle);
        hierarchy_window('sequence_display',hhandles);
    end;
end;

function handles=unselect(handles,args)
    
global model
global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: unselect address');
    return
end;

command=sprintf('unselect %s',strtrim(args));
undo_cmd=sprintf('select %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));
indices=resolve_address(address);

if strcmp(strtrim(address),'*') % unselect all
    unselect_all_secondary;
end;

if ~isfield(model,'selections') || ~isfield(model,'selected');
    add_msg_board('Unselect failed as nothing is selected.');
    model.selections=0;
    model.selected={};
end;

if model.selections==0 || isempty(model.selected);
    add_msg_board('Unselect failed as nothing is selected.');
end;

unlight_selection;
[m,n]=size(indices);

for k=1:m,
    sec_indices=resolve_secondary_address(address,indices(k,:));
    if ~isempty(sec_indices),
        snum=indices(k,1);
        cnum=indices(k,2);
        switch sec_indices(1)
            case 0
                model.structures{snum}(cnum).loop_defs{sec_indices(2)}.selected=0;
            case 1
                model.structures{snum}(cnum).helix_defs{sec_indices(2)}.selected=0;
            case 2
                model.structures{snum}(cnum).sheet_defs{sec_indices(2)}.selected=0;
        end;
    end;
end;

new_selections=0;
new_selected={};
for kk=1:model.selections,
    sindices=model.selected{kk};
    sdepth=length(find(sindices>0)); % determine type of current object
    sindices=sindices(1:sdepth);
    store=1;
    for k=1:m,
        idepth=length(find(indices(k,:)>0)); % determine type of current object
        cindices=indices(k,1:idepth);
        if sdepth==idepth
            if cindices==sindices
                store=0;
            end;
        end;
    end;
    if store
        new_selections=new_selections+1;
        new_selected{new_selections}=sindices;
    end;
end;
model.selections=new_selections;
model.selected=new_selected;
highlight_selection;
% Display of current selection

if model.selections>=1,
    add_msg_board('--- Selected ---');
    for k=1:model.selections,
        [message,prev_info]=get_object(model.selected{k},'info');
        add_msg_board(prev_info{1});
    end;
else
    add_msg_board('--- Nothing selected ---');
end;

if hMain.hierarchy_display,
%     if hMain.large,
%         hwhandle=hMain.hierarchy_window_large;
%         hhandles=guidata(hwhandle);
%         hierarchy_window_large('sequence_display',hhandles);
%     else
%         hwhandle=hMain.hierarchy_window;
%         hhandles=guidata(hwhandle);
%         hierarchy_window('sequence_display',hhandles);
%     end;        
    hwhandle=hMain.hierarchy_window;
    hhandles=guidata(hwhandle);
    hierarchy_window('sequence_display',hhandles);
end;


function handles=zoom(handles,args)

global hMain;
global model

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: zoom in|out');
    return
end;

myargs=textscan(args,'%s');
mode=char(myargs{1}(1));

all_graphics={};

switch lower(mode)
    case 'in'
        command=sprintf('zoom in');
        undo_cmd=sprintf('zoom out');
        [handles,veto]=cmd_history(handles,command,undo_cmd);
        indices=resolve_address('*');
        [m,n]=size(indices);
        poi=0;
        for ko=1:m, % loop over all objects
            idepth=length(find(indices(ko,:)>0)); % determine type of current object
            cindices=indices(ko,1:idepth);
            [message,allindices]=get_object(cindices,'descendants');
            [ma,na]=size(allindices);
            for ka=0:ma, % ka=0 takes care of the selected object, k>0 of its descendants
                if ka>0,
                    cindices=allindices(ka,:);
                end;
                cindices=cindices(cindices>0);
                if length(cindices)==3, % selection of secondary structure elements
                    if isfield(model.structures{cindices(1)}(cindices(2)).residues{cindices(3)},'secondary_graphics')
                        secgraphics=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.secondary_graphics;
                        if ~isempty(secgraphics),
                            for k=1:length(secgraphics),
                                graphics=secgraphics(k);
                                if ~isempty(graphics),
                                    poi=poi+1;
                                    all_graphics{poi}=graphics;
                                end;
                            end;
                        end;
                    end;
                end;
                if ~isempty(cindices)
                    [msg,graphics]=get_object(cindices,'graphics');
                    if ~isempty(graphics),
                        poi=poi+1;
                        all_graphics{poi}=graphics;
                    end;
                end;
            end;
        end;
        m=length(all_graphics);
        if m<1,
            add_msg_board('### Warning ###: Selection has no displayed graphics.');
            return;
        end;
        ghandles=zeros(1,10000);
        poi=0;
        fail=1;
        for k=1:m,
            if isfield(all_graphics{k},'objects'),
                if ~isempty(all_graphics{k}.objects),
                    poi2=poi+length(all_graphics{k}.objects);
                    ghandles(poi+1:poi2)=all_graphics{k}.objects;
                    poi=poi2;
                    fail=0;
                end;
            end;
        end;
        if fail,
            add_msg_board('### Warning ### Selection has no displayed graphics.');
        else
            ghandles=ghandles(1:poi2);
%             fprintf(1,'Starting camlookat\n');
%             profile on
%             tic,
            camlookat(ghandles);
%             toc,
%             fprintf(1,'Back from camlookat\n');
%             profile viewer
        end;
    case 'out'
        command=sprintf('zoom out');
        undo_cmd=sprintf('zoom in');
        [handles,veto]=cmd_history(handles,command,undo_cmd);
        camlookat(hMain.axes_model);
    otherwise
        add_msg_board('ERROR: Unrecognized zoom mode.');
        add_msg_board('type "help zoom" to get on-line help.');
end;

function handles=undefine(handles,args)

global model
global hMain

if model.locked,
    add_msg_board('ERROR: Must unlock secondary structure before undefining');
    return;
end;

command=sprintf('undefine %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);
if veto,
    add_msg_board('Secondary structure undefinition cancelled by user.');
    return
end;

if ~isempty(args) && ~isempty(strtrim(args)),
    myargs=textscan(args,'%s');
    address=strtrim(char(myargs{1}(1)));
    indices=resolve_address(address);
    if length(indices)==1,
        [msg,indices]=get_object(indices,'children');
    end;
else
    indices=zeros(1,2);
    indices(1)=model.current_structure;
    cid=tag2id(model.current_chain,model.chain_tags{model.current_structure});
    indices(2)=cid;
end;
[m,n]=size(indices);
if n~=2,
    add_msg_board('ERROR: Supplied address must be for a structure or chain.');
    return;
end;
if ~isempty(indices),
    if hMain.hierarchy_display,
%         if hMain.large,
%             close(hMain.hierarchy_window_large);
%             hMain.hierarchy_display=0;
%         else
%             close(hMain.hierarchy_window);
%             hMain.hierarchy_display=0;
%         end;
        close(hMain.hierarchy_window);
        hMain.hierarchy_display=0;
    end;    
    for k=1:m,
        set_chain(indices(k,:),'undefine');
    end;
else
    add_msg_board('ERROR: No chain addressed for secondary structure undefine.');
end;

% function handles=xyz(handles,args)
% 
% myargs=textscan(args,'%s');
% address=char(myargs{1}(1));
% [message,xyz]=get_object(address,'xyz');
% disp(xyz);

function handles=mass(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: mass address');
    return
end;


myargs=textscan(args,'%s');
address=char(myargs{1}(1));
if length(myargs{1})>1,
    option=myargs{1}(2);
    if ~strcmp(option,'-water'),
        option='-nowater';
    end;
else
    option='-nowater';
end;
[message,mass]=get_object(address,'mass',option);
if message.error,
    add_msg_board(message.text);
    return
end;
indices=resolve_address(address); % to care about symbolic address * and abbreviated addresses
address=mk_address(indices);

if iscell(mass),
    masses=cat(1,mass{:});
    mass=sum(masses);
end,
add_msg_board(sprintf('Mass of %s is %6.2f g/mol',address,mass));


function handles=help_cmd(handles,args,commands)


global help_files

if ~isempty(args) && ~isempty(strtrim(args)),
    myargs=textscan(args,'%s');
    command=strtrim(char(myargs{1}(1)));
    id=tag2id(command,commands);
    if ~isempty(id),
        entry=strcat(help_files,'commands.html#',command);
        webcall(entry,'-helpbrowser');
    else
        add_msg_board(sprintf('ERROR: Command "%s" does not exist.',command));
    end;
else
    entry=strcat(help_files,'commands.html');
    webcall(entry,'-helpbrowser');
    add_msg_board('--- List of commands ---');
    msg='';
    fail=0;
    id=0;
    while ~fail,
        id=id+1;
        com=id2tag(id,commands);
        if isempty(com),
            fail=1;
            add_msg_board(msg);
        else
            if length(msg)+length(com)<60,
                if isempty(msg),
                    msg=com;
                else
                    msg=sprintf('%s %s',msg,com);
                end;
            else
                add_msg_board(msg);
                msg=com;
            end;
        end;
    end;
    add_msg_board('Use "help cmd" to get online help on command "cmd"');
end;

function handles=def_view(handles,args)

global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: view x|-x|y|-y|z|-z|memory');
    add_msg_board('or   : view xval yval zval');
    add_msg_board('or   : view xval yval zval xup yup zup');
    return
end;


undo_cmd=sprintf('noundo');


myargs=textscan(args,'%s');

if length(myargs{1})<3,
	mode=char(myargs{1}(1));
    switch mode
        case {'x','+x'}
            command=sprintf('view x');
            handles=cmd_history(handles,command,undo_cmd);
            set_view([1,0,0]);
        case '-x'
            command=sprintf('view -x');
            handles=cmd_history(handles,command,undo_cmd);
            set_view([-1,0,0]);
        case {'y','+y'}
            command=sprintf('view y');
            handles=cmd_history(handles,command,undo_cmd);
            set_view([0,1,0]);
        case '-y'
            command=sprintf('view -y');
            handles=cmd_history(handles,command,undo_cmd);
            set_view([0,-1,0]);
        case {'z','+z'}
            command=sprintf('view z');
            handles=cmd_history(handles,command,undo_cmd);
            set_view([0,0,1]);
        case '-z'
            command=sprintf('view -z');
            handles=cmd_history(handles,command,undo_cmd);
            set_view([0,0,-1]);
        case 'memory'
            set(hMain.axes_model,'CameraUpVector',hMain.view_memory_camup);
            set(hMain.axes_model,'CameraTarget',hMain.view_memory_camtar);
            set(hMain.axes_model,'CameraViewAngle',hMain.view_memory_camview);
            set(hMain.axes_model,'CameraPosition',hMain.view_memory_campos);
            set(hMain.axes_model,'Projection',hMain.view_memory_camproj);
            set_view;
        otherwise
            add_msg_board(sprintf('ERROR: View direction "%s" is not defined.',mode));
    end;
    if ~strcmp(mode,'memory'),
        camlookat(hMain.axes_model);
    end;
else
    x=str2double(char(myargs{1}(1)));
    y=str2double(char(myargs{1}(2)));
    z=str2double(char(myargs{1}(3)));
    if ~isnan(x) && ~isnan(y) && ~isnan(z),
        vec=[x,y,z];
        nv=norm(vec);
        command=sprintf('view %5.3f %5.3f %5.3f',x/nv,y/nv,z/nv);
        if length(myargs{1})>=6,
            xup=str2double(char(myargs{1}(4)));
            yup=str2double(char(myargs{1}(5)));
            zup=str2double(char(myargs{1}(6)));
            if ~isnan(xup) && ~isnan(yup) && ~isnan(zup),
                vec=[x,y,z,xup,yup,zup];
                nv2=norm([xup,yup,zup]);
                command=sprintf('view %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f',x/nv,y/nv,z/nv,xup/nv2,yup/nv2,zup/nv2);
                handles=cmd_history(handles,command,undo_cmd);
            else
                add_msg_board('Warning: Camera up arguments must be three numbers.');
                add_msg_board('Camera up vector is unchanged.');
            end;
        else
            handles=cmd_history(handles,command,undo_cmd);
        end;
        set_view(vec);
    else
        add_msg_board('ERROR: Vector argument must contain three or six numbers.');
    end;
end;

function handles=def_camup(handles,args)

global hMain

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: camup x|-x|y|-y|z|-z');
    add_msg_board('or   : camup xup yup zup');
    return
end;


myargs=textscan(args,'%s');

if length(myargs{1})<3,
	mode=char(myargs{1}(1));
    switch mode
        case {'x','+x'}
            set_camup([1,0,0]);
        case '-x'
            set_camup([-1,0,0]);
        case {'y','+y'}
            set_camup([0,1,0]);
        case '-y'
            set_camup([0,-1,0]);
        case {'z','+z'}
            set_camup([0,0,1]);
        case '-z'
            set_camup([0,0,-1]);
        otherwise
            add_msg_board(sprintf('ERROR: Camera up direction "%s" is not defined.',mode));
    end;
    camlookat(hMain.axes_model);
else
    xup=str2double(char(myargs{1}(1)));
    yup=str2double(char(myargs{1}(2)));
    zup=str2double(char(myargs{1}(3)));
    if ~isnan(xup) && ~isnan(yup) && ~isnan(zup),
        vec=[xup,yup,zup];
        set_camup(vec);
    else
        add_msg_board('ERROR: Vector argument must contain three numbers.');
    end;
end;

function handles=unlabel(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: unlabel address');
    add_msg_board('where address must resolve to labeled residues.');
    return
end;

command=sprintf('unlabel %s',strtrim(args));
undo_cmd=sprintf('label %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));

[indices,msg]=resolve_address(address);
if isempty(indices),
    add_msg_board('No objects selected by this address.');
else
    residues=0;
    [m,n]=size(indices);
    for k=1:m,
        cindices=indices(k,:);
        cindices=cindices(cindices>0);
        if length(cindices)==4,
            residues=residues+1;
            set_residue(cindices,'unlabel');
        end;
    end;
    if residues>0,
        add_msg_board(sprintf('%i residue(s) unlabeled.',residues));
    else
        add_msg_board('No residues selected by this address. Nothing unlabeled.');
    end;
end;

function handles=write_sequence(handles,args)

global residue_defs
global hMain
global general
global model

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: wrseq address');
    return
end;

command=sprintf('wrseq %s',strtrim(args));
undo_cmd=sprintf('wrseq %s',strtrim(args));
[handles,veto]=cmd_history(handles,command,undo_cmd);

myargs=textscan(args,'%s');
address=char(myargs{1}(1));

[indices,msg]=resolve_address(address);
if isempty(indices),
    add_msg_board('ERROR: No objects selected by this address.');
else
    [m,n]=size(indices);
    if n~=4,
        add_msg_board('ERROR: Selected objects are not residues.');
    else
        seq='';
        for k=1:m,
            tlc=model.structures{indices(k,1)}(indices(k,2)).residues{indices(k,3)}.info(indices(k,4)).name;
            aa=tag2id(tlc,upper(residue_defs.restags));
            slc=residue_defs.single_letter_code(aa);
            seq=[seq slc];
        end;
        add_msg_board('Sequence for selection:');
        add_msg_board(seq);
        outfile=[general.tmp_files 'seq.txt'];
        fid=fopen(outfile,'wt');
        fprintf(fid,'%s\n',seq);
        fclose(fid);
        hMain.report_file=outfile;
        report_editor;
    end;        
end;

function handles=report(handles,args)

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: report filename information [...]');
    add_msg_board('Usage: report filename label address [brief/full]');
    add_msg_board('Usage: report filename distance adr1 adr2 [brief/full]');
    add_msg_board('Only residue addresses are allowed.');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) < 2,
    add_msg_board('ERROR: Command requires at least two arguments.');
    return
end;

fname = char(myargs{1}(1));
information = char(myargs{1}(2));

error_type = 0;
switch information
    case 'label'
        if length(myargs{1}) < 3,
            add_msg_board('ERROR: Residue address missing.');
            return;
        else
            adr1=char(myargs{1}(3));
        end;
        [msg,label] = get_object(adr1,'label');
        if msg.error,
            add_msg_board('ERROR: No label information.');
            if isfield(msg,'txt'),
                add_msg_board(msg.txt);
            end;
            return
        end;
        if length(myargs{1}) > 3,
            mode = char(myargs{1}(4));
        else
            mode = 'brief';
        end;
        fid = fopen(fname,'at');
        if fid == -1,
            add_msg_board('ERROR: Output file could not be opened.');
            return;
        end;
        switch mode
            case 'brief'
                for k = 1:length(label),
                    fprintf(fid,'%5.3f\t%5.3f\t%5.3f\t%5.3f %% x,y,z,rmsd [nm] for %s (%s) at %i K\n',label(k).xyz,label(k).rmsd,adr1,label(k).name,label(k).T);
                end;
            case 'full'
                for k = 1:length(label),
                    fprintf(fid,'%% Residue %s labeled with %s at %i K:\n',adr1,label(k).name,label(k).T);
                    fprintf(fid,'%% xNO/nm\tyNO/nm\tzNO/nm\tpopulation\n');
                    [m,n] = size(label(k).NOpos);
                    for kk = 1:m,
                        fprintf(fid,'  %6.3f\t%6.3f\t%6.3f\t%6.3f\n',label(k).NOpos(kk));
                    end;
                end;
            otherwise
                add_msg_board(sprintf('ERROR: Unknown output mode: %s',mode));            
        end;
        fclose(fid);
    case 'distance'
        if length(myargs{1}) < 4,
            add_msg_board('ERROR: Distance report requires two residue addresses.');
            return;
        else
            adr1=char(myargs{1}(3));
            adr2=char(myargs{1}(4));
        end;
        [msg,label1] = get_object(adr1,'label');
        if msg.error,
            add_msg_board(sprintf('ERROR: No label information for %s.',adr1));
            add_msg_board(msg.txt);
            return
        end;
        [msg,label2] = get_object(adr2,'label');
        if msg.error,
            add_msg_board(sprintf('ERROR: No label information for %s.',adr2));
            add_msg_board(msg.txt);
            return
        end;
        if length(myargs{1}) > 4,
            mode = char(myargs{1}(5));
        else
            mode = 'brief';
        end;
        [rax,distr]=get_distribution(label1(1).NOpos,label2(1).NOpos,0.1);
        fid = fopen(fname,'at');
        if fid == -1,
            add_msg_board('ERROR: Output file could not be opened.');
            return;
        end;
        switch mode
            case 'brief'
                for k = 1:length(rax),
                    fprintf(fid,'%6.2f%8.4f\n',rax(k),distr(k));
                end;
            case 'full'
                fprintf(fid,'%% Residue %s labeled with %s at %i K:\n',adr1,label1(1).name,label1(1).T);
                fprintf(fid,'%% Residue %s labeled with %s at %i K:\n',adr2,label2(1).name,label2(1).T);
                fprintf(fid,'%% r/nm  P(r)\n');
                for k = 1:length(rax),
                    fprintf(fid,'%6.2f%8.4f\n',rax(k),distr(k));
                end;
            otherwise
                add_msg_board(sprintf('ERROR: Unknown output mode: %s',mode));            
        end;
        fclose(fid);
    otherwise
    add_msg_board(sprintf('ERROR: Unknown information request: %s',information));
end;

adr1=char(myargs{1}(1));
adr2=char(myargs{1}(2));
[msg1,coor1] = get_object(adr1,'coor');
[msg2,coor2] = get_object(adr2,'coor');

function handles=distance(handles,args)

global hMain
global rotamer_libraries

command=sprintf('distance %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: distance adr1 adr2 label [filename]');
    add_msg_board('Only residue addresses are allowed.');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) < 3,
    add_msg_board('ERROR: Command requires at least three arguments.');
    return
end;


adr1=char(myargs{1}(1));
adr2=char(myargs{1}(2));
label=char(myargs{1}(3));

rindices1 = resolve_address(adr1);
[m1,n1] = size(rindices1);
if m1 ~= 1 || n1 ~=4,
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr1));
    return;
end;
rindices2 = resolve_address(adr2);
[m2,n2] = size(rindices2);
if m2 ~= 1 || n2 ~=4,
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr2));
    return
end;

mode = '';
rot_lib_name = '';

for k = 1:length(rotamer_libraries),
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label),
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T),
            if rotamer_libraries(k).T(kk) == 298,
                id = kk;
            end;
        end;
        if id >0,
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
            mode = 'single';
        end;
        if isfield(rotamer_libraries(k),'MC') && ~isempty(rotamer_libraries(k).MC),
            rot_lib_name = rotamer_libraries(k).MC;
            mode = 'MC';
        end;
    end;
end;

if length(myargs{1}) > 3,
    fname = char(myargs{1}(4));
else
    fname = '';
end;

set(hMain.MMM,'Pointer','watch');
drawnow

switch mode
    case 'single'
        load(rot_lib_name);
        label1 = rotamer_populations(rindices1,rot_lib);
        label2 = rotamer_populations(rindices2,rot_lib);
        [rax,distr]=get_distribution(label1(1).NOpos,label2(1).NOpos,0.1);
        distr = distr/sum(distr);
        rmean = 10*sum(rax.*distr);
        add_msg_board(sprintf('Mean distance : %5.1f Å',rmean));
        rdev = 10*rax - rmean;
        mom2 = sum(rdev.^2.*distr);
        add_msg_board(sprintf('Std. deviation: %5.1f Å',sqrt(mom2)));
        if ~isempty(fname),
            fid = fopen(fname,'at');
            if fid == -1,
                add_msg_board('ERROR: Output file could not be opened.');
                return;
            end;
            for k = 1:length(rax),
                fprintf(fid,'%6.2f%8.4f\n',rax(k),distr(k));
            end;
            fclose(fid);
        end;
    case 'MC'
        poi = strfind(rot_lib_name,'#');
        for k = 1:5,
            clib = sprintf('%s%i%s',rot_lib_name(1:poi-1),k,rot_lib_name(poi+1:end));
            load(clib);
            label1 = rotamer_populations(rindices1,rot_lib);
            label2 = rotamer_populations(rindices2,rot_lib);
            if isempty(label1) || isempty(label2),
                add_msg_board('ERROR: Site is too tight to be labelled with all libraries.');
                return
            end;
            NOpos1{k} = label1.NOpos;
            NOpos2{k} = label2.NOpos;
        end;
        [rax,distr]=get_distribution(NOpos1{1},NOpos2{1},0.1);
        distributions = zeros(25,length(distr));
        rmean = zeros(1,25);
        stddev = zeros(1,25);
        mdistr = zeros(size(distr));
        for k1 = 1:5,
            for k2 = 1:5,
                no = 5*(k1-1)+k2;
                [~,distr]=get_distribution(NOpos1{k1},NOpos2{k2},0.1);
                distr = distr/sum(distr);
                mdistr = mdistr + distr;
                distributions(no,:) = distr;
                rmean(no) = 10*sum(rax.*distr);
                % fprintf(1,'%i: <r> = %8.4f\n',no,rmean(no));
                rdev = 10*rax - rmean(no);
                stddev(no) = sqrt(sum(rdev.^2.*distr));
                % fprintf(1,'%i: sr = %8.4f\n',no,stddev(no));
            end;
        end;
        mdistr = mdistr/sum(mdistr);
        mrmean = 10*sum(rax.*mdistr);
        mdiff = std(rmean);
        mrdev = 10*rax - mrmean;
        mstddev = sqrt(sum(mrdev.^2.*mdistr));
        msdiff = std(stddev);
        add_msg_board(sprintf('Mean distance : %5.1f Å predicted with uncertainty %5.1f Å',mrmean,mdiff));
        add_msg_board(sprintf('Std. deviation: %5.1f Å predicted with uncertainty %5.1f Å',mstddev,msdiff));
        if ~isempty(fname),
            fid = fopen(fname,'at');
            if fid == -1,
                add_msg_board('ERROR: Output file could not be opened.');
                return;
            end;
            for k = 1:length(rax),
                fprintf(fid,'%6.2f',rax(k));
                for kk = 1:25,
                    fprintf(fid,'%8.4f',distributions(kk,k));
                end;
                fprintf(fid,'\n');
            end;
            fclose(fid);
        end;
    otherwise
        add_msg_board(sprintf('ERROR: No rotamer library for label %s.',label));
end;
set(hMain.MMM,'Pointer','arrow');
drawnow

function handles = coil_statistics(handles,args)

global hMain
global rotamer_libraries

command=sprintf('statistics %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: statistics adr label [filename]');
    add_msg_board('Only chain addresses are allowed, the chain must be continuous.');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) < 2,
    add_msg_board('ERROR: Command requires at least two arguments.');
    return
end;


adr = char(myargs{1}(1));
label=char(myargs{1}(2));
adr1 = strcat(adr,'{:}.CA'); % find all residues that have Calpha atoms

rindices1 = resolve_address(adr1);
[m1,n1] = size(rindices1);
if m1 < 1 || n1 ~= 5,
    add_msg_board(sprintf('ERROR: Address %s could not be resolved to a set of residues.',adr));
    return;
end;

rot_lib_name = '';

for k = 1:length(rotamer_libraries),
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label),
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T),
            if rotamer_libraries(k).T(kk) == 298,
                id = kk;
            end;
        end;
        if id > 0,
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
        end;
    end;
end;

if isempty(rot_lib_name),
    add_msg_board(sprintf('ERROR: No rotamer library at 298 K could be found for label %s.',label));
    return;
end;

if length(myargs{1}) > 2,
    fname = char(myargs{1}(3));
else
    fname = '';
end;

set(hMain.MMM,'Pointer','watch');
drawnow

load(rot_lib_name);
label1 = rotamer_populations(rindices1(:,1:4),rot_lib);

rindices2 = zeros(length(label1),4);
NOpos = zeros(length(label1),3);
for k = 1:length(label1),
    rindices2(k,:) = label1(k).indices;
    pop = label1(k).NOpos(:,4);
    xyz = label1(k).NOpos(:,1:3);
    pop=pop/sum(pop);
    xyz=pop'*xyz;
    NOpos(k,:) = xyz;
end;

max_res_dist = 0;
for k1 = 1:length(label1)-1,
    for k2 = k1+1:length(label1),
        if sum(abs(rindices2(k1,1:3)-rindices2(k2,1:3))) == 0, % same structure, same chain, same model
            res_dist = abs(rindices2(k1,4)-rindices2(k2,4));
            if res_dist > max_res_dist,
                max_res_dist = res_dist;
            end;
        end;
    end;
end;

ra = 0;
re = 150;
n = 150;
rax = linspace(ra,re,n+1);
histo = zeros(max_res_dist,n+1);

for k1 = 1:length(label1)-1,
    for k2 = k1+1:length(label1),
        if sum(abs(rindices2(k1,1:3)-rindices2(k2,1:3))) == 0, % same structure, same chain, same model
            res_dist = abs(rindices2(k1,4)-rindices2(k2,4));
            r = norm(NOpos(k1,:)-NOpos(k2,:));
            poi = 1 + round(n*(r-ra)/(re-ra));
            if poi > 0 && poi <= n+1,
                histo(res_dist,poi) = histo(res_dist,poi) + 1;
            end;
        end;
    end;
end;

mask = (rax >= 20).*(rax <= 60);
rdax = 1:max_res_dist;
idr = zeros(1,max_res_dist);
for k = 1:max_res_dist,
    hm = histo(k,:).*mask;
    idr(k) = sum(hm)/sum(histo(k,:));
end;

if ~isempty(fname),
    save(fname,'histo','rax','idr','rdax');
end;
set(hMain.MMM,'Pointer','arrow');
drawnow


function handles = get_beacons(handles,args)

global hMain
global rotamer_libraries

command=sprintf('beacons %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: beacons anchor label num [rmin] [rmax] [n]');
    add_msg_board('Only an address referring to one or two residues is allowed for anchor.');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) < 3,
    add_msg_board('ERROR: Command requires at least three arguments.');
    add_msg_board('Usage: beacons anchor label num [rmin] [rmax] [n]');
    add_msg_board('Only an address referring to one or two residues is allowed for anchor.');
    return
end;


anchor = char(myargs{1}(1));
label = char(myargs{1}(2));
nb = str2double(char(myargs{1}(3))); % number of beacons

if nb < 4 || nb > 6,
    add_msg_board('ERROR: Beacon search implemented only for 4 to 6 beacons.');
    return
end;

if length(myargs{1}) > 3, % minimum distance (nm) for a beacon to be considered
    rmin = str2double(char(myargs{1}(4)));
else
    rmin = 2.5; % defaults to 2.5 nm
end;
if rmin <= 0,
    rmin = 0.3;
end;

if length(myargs{1}) > 4, % maximum distance (nm) for a beacon to be considered
    rmax = str2double(char(myargs{1}(5)));
else
    rmax = 3.5; % defaults to 3.5 nm
end;

if length(myargs{1}) > 5, % number of proposals to be computed
    nmax = str2double(char(myargs{1}(6)));
else
    nmax = 20; % defaults to 20 sets
end;

aindices = resolve_address(anchor);
[ma,n1] = size(aindices);
if ma < 1 || n1 ~= 4,
    add_msg_board(sprintf('ERROR: Address %s could not be resolved to a residue address.',anchor));
    return;
end;


adr = mk_address(aindices(1,1:3));
fname = sprintf('beacons_%s_%s_%i.dat',adr,label,nb);

rindices1 = zeros(10000,4);
chains_done = zeros(ma,3);
poi = 0;
for k = 1:ma, % find all amino acid residue indices for all chains that contain anchor residues
    done = false;
    for kd = 1:k-1,
        if sum(abs(aindices(k,1:3)-chains_done(kd,:)))==0,
            done = true;
        end;
    end;
    if ~done,
        adr = mk_address(aindices(k,1:3));
        chains_done(k,:) = aindices(k,1:3);
        adr1 = strcat(adr,'.CA'); % find all residues that have Calpha atoms

        rindices0 = resolve_address(adr1);
        [m1,n1] = size(rindices0);
        if m1 < 1 || n1 ~= 5,
            add_msg_board(sprintf('ERROR: Address %s could not be resolved to a set of residues.',adr));
            return;
        end;
        rindices1(poi+1:poi+m1,:) = rindices0(:,1:4);
        poi = poi + m1;
    end;
end;
rindices1 = rindices1(1:poi,:);

% create virtual labels

rot_lib_name = '';

for k = 1:length(rotamer_libraries),
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label),
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T),
            if rotamer_libraries(k).T(kk) == 298,
                id = kk;
            end;
        end;
        if id > 0,
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
        end;
    end;
end;

if isempty(rot_lib_name),
    add_msg_board(sprintf('ERROR: No rotamer library at 298 K could be found for label %s.',label));
    return;
end;

set(hMain.MMM,'Pointer','watch');
drawnow

load(rot_lib_name);
label1 = rotamer_populations(rindices1(:,1:4),rot_lib);

% find the indices of all residues that could be labelled and the
% corresponding label positions
rindices2 = zeros(length(label1),4);
NOpos = zeros(length(label1),3);
NOanchor = zeros(1,ma);
apoi = 0;
for k = 1:length(label1),
    rindices2(k,:) = label1(k).indices;
    for ka = 1:ma,
        if sum(abs(label1(k).indices - aindices(ka,:))) == 0, % was an anchor residue found?
            apoi = apoi+1;
            NOanchor(apoi) = k;
        end;
    end;
    pop = label1(k).NOpos(:,4);
    xyz = label1(k).NOpos(:,1:3);
    pop=pop/sum(pop);
    xyz=pop'*xyz;
    NOpos(k,:) = xyz;
end;
NOanchor = NOanchor(1:apoi);

if length(NOanchor) < ma,
    add_msg_board(sprintf('ERROR: Not all anchor residues could be labeled for address %s.',anchor));
    return;
end;

xyz_anchor = NOpos(NOanchor,:);

% Prune label list to sites within accepted distance range and with at
% least two rotamers and a partition function of at least 0.1
poi = 0;
NOsites = zeros(length(label1),3);
sindices = zeros(length(label1),4);
r = zeros(1,apoi);
for k = 1:length(label1),
    for ka = 1:apoi,
        r(ka) = norm(NOpos(k,:)-xyz_anchor(ka,:))/10; % conversion Å to nm
    end;
    m = length(label1(k).rotamers); % number of rotamers
    z = label1(k).partition_function;
    if min(r) >= rmin && max(r) <= rmax && m >= 2 && z >= 0.1,
        poi = poi + 1;
        NOsites(poi,:) = NOpos(k,:);
        sindices(poi,:) = rindices2(k,:);
    end;
end;

NOsites = NOsites(1:poi,:);
sindices = sindices(1:poi,:);

add_msg_board(sprintf('%i possible beacon sites found for %i anchors. Optimizing PDOP...',poi,apoi));

% create the a vectors
avec = ones(poi,4);
for a = 1:apoi,
    avec_set{a} = avec;
end;
rsite = zeros(apoi,poi);
for a = 1:apoi,
    avec = avec_set{a};
    for j = 1:poi,
        rj = norm(NOsites(j,:)-xyz_anchor(a,:));
        rsite(a,j) = rj;
        avec(j,1:3) = (NOsites(j,:)-xyz_anchor(a,:))/rj;
    end;
    avec_set{a} = avec;
end;

beacons = zeros(nmax,4*nb); % beacon index array
nom = 1e6*ones(1,nmax); % number of merit
maxnom = 1e6;

H = zeros(nb,4);
for ka = 1:apoi,
    Hset{ka} = H;
end;
config = zeros(1,4*nb);
tic,
switch nb
    case 4
        for k1 = 1:poi-3,
            for ka = 1:apoi,
                avec = avec_set{a};
                H = Hset{ka};
                H(1,:) = avec(k1,:);
                Hset{ka} = H;
            end;
            config(1:4) = sindices(k1,:);
            for k2 = k1+1:poi-2,
                for ka = 1:apoi,
                    avec = avec_set{a};
                    H = Hset{ka};
                    H(2,:) = avec(k2,:);
                    Hset{ka} = H;
                end;
                config(5:8) = sindices(k2,:);
                for k3 = k2+1:poi-1,
                    for ka = 1:apoi,
                        avec = avec_set{a};
                        H = Hset{ka};
                        H(3,:) = avec(k3,:);
                        Hset{ka} = H;
                    end;
                    config(9:12) = sindices(k3,:);
                    for k4 = k3+1:poi,
                        for ka = 1:apoi,
                            avec = avec_set{a};
                            H = Hset{ka};
                            H(4,:) = avec(k4,:);
                            Hset{ka} = H;
                        end;
                        config(13:16) = sindices(k4,:);
                        cnom = 0;
                        for ka = 1:apoi,
                            H = Hset{ka};
                            cnom = cnom + trace(inv(H'*H));
                        end;
                        cnom = sqrt(cnom/apoi);
                        if cnom < maxnom,
                            nom(nmax) = cnom;
                            beacons(nmax,:) = config;
                            [~,sind] = sort(nom);
                            nom = nom(sind);
                            beacons = beacons(sind,:);
                            maxnom = nom(nmax);
                        end;
                    end;
                end;
            end;
        end;
    case 5
        for k1 = 1:poi-4,
            for ka = 1:apoi,
                avec = avec_set{a};
                H = Hset{ka};
                H(1,:) = avec(k1,:);
                Hset{ka} = H;
            end;
            config(1:4) = sindices(k1,:);
            for k2 = k1+1:poi-3,
                for ka = 1:apoi,
                    avec = avec_set{a};
                    H = Hset{ka};
                    H(2,:) = avec(k2,:);
                    Hset{ka} = H;
                end;
                config(5:8) = sindices(k2,:);
                for k3 = k2+1:poi-2,
                    for ka = 1:apoi,
                        avec = avec_set{a};
                        H = Hset{ka};
                        H(3,:) = avec(k3,:);
                        Hset{ka} = H;
                    end;
                    config(9:12) = sindices(k3,:);
                    for k4 = k3+1:poi-1,
                        for ka = 1:apoi,
                            avec = avec_set{a};
                            H = Hset{ka};
                            H(4,:) = avec(k4,:);
                            Hset{ka} = H;
                        end;
                        config(13:16) = sindices(k4,:);
                        for k5 = k4+1:poi,
                            for ka = 1:apoi,
                                avec = avec_set{a};
                                H = Hset{ka};
                                H(5,:) = avec(k5,:);
                                Hset{ka} = H;
                            end;
                            config(17:20) = sindices(k5,:);
                            cnom = 0;
                            for ka = 1:apoi,
                                H = Hset{ka};
                                cnom = cnom + trace(inv(H'*H));
                            end;
                            cnom = sqrt(cnom/apoi);
                            if cnom < maxnom,
                                nom(nmax) = cnom;
                                beacons(nmax,:) = config;
                                [~,sind] = sort(nom);
                                nom = nom(sind);
                                beacons = beacons(sind,:);
                                maxnom = nom(nmax);
                            end;
                        end;
                    end;
                end;
            end;
        end;
    case 6
        for k1 = 1:poi-5,
            for ka = 1:apoi,
                avec = avec_set{a};
                H = Hset{ka};
                H(1,:) = avec(k1,:);
                Hset{ka} = H;
            end;
            config(1:4) = sindices(k1,:);
            for k2 = k1+1:poi-4,
                for ka = 1:apoi,
                    avec = avec_set{a};
                    H = Hset{ka};
                    H(2,:) = avec(k2,:);
                    Hset{ka} = H;
                end;
                config(5:8) = sindices(k2,:);
                for k3 = k2+1:poi-3,
                    for ka = 1:apoi,
                        avec = avec_set{a};
                        H = Hset{ka};
                        H(3,:) = avec(k3,:);
                        Hset{ka} = H;
                    end;
                    config(9:12) = sindices(k3,:);
                    for k4 = k3+1:poi-2,
                        for ka = 1:apoi,
                            avec = avec_set{a};
                            H = Hset{ka};
                            H(4,:) = avec(k4,:);
                            Hset{ka} = H;
                        end;
                        config(13:16) = sindices(k4,:);
                        for k5 = k4+1:poi-1,
                            for ka = 1:apoi,
                                avec = avec_set{a};
                                H = Hset{ka};
                                H(5,:) = avec(k5,:);
                                Hset{ka} = H;
                            end;
                            config(17:20) = sindices(k5,:);
                            for k6 = k5+1:poi,
                                for ka = 1:apoi,
                                    avec = avec_set{a};
                                    H = Hset{ka};
                                    H(6,:) = avec(k6,:);
                                    Hset{ka} = H;
                                end;
                                config(21:24) = sindices(k6,:);
                                cnom = 0;
                                for ka = 1:apoi,
                                    H = Hset{ka};
                                    cnom = cnom + trace(inv(H'*H));
                                end;
                                cnom = sqrt(cnom/apoi);
                                if cnom < maxnom,
                                    nom(nmax) = cnom;
                                    beacons(nmax,:) = config;
                                    [~,sind] = sort(nom);
                                    nom = nom(sind);
                                    beacons = beacons(sind,:);
                                    maxnom = nom(nmax);
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
end;
toc,

fid = fopen(fname,'wt');
fprintf(fid,'%% Optimum beacon sites for localization of residue %s\n',anchor);
fprintf(fid,'%% Predicted distances between %4.1f and %4.1f nm\n',rmin,rmax);
fprintf(fid,'%% ');
for k = 1:nb,
    fprintf(fid,'site %i address  ',k);
end;
fprintf(fid,'Figure of merit\n');
for j = 1:nmax,
    for k = 1:nb,
        bas = 1 + 4*(k-1);
        adr = mk_address(beacons(j,bas:bas+3));
        fprintf(fid,'%16s',adr);
    end;
    fprintf(fid,'%15.3f\n',nom(j));
end;
fclose(fid);

set(hMain.MMM,'Pointer','arrow');
drawnow

hMain.report_file = fname;
report_editor;


function handles = ensemble_rmsd(handles,args)

global model

command=sprintf('rmsd %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: rmsd adr [type]');
    add_msg_board('type can be ''all'' or ''backbone'' and defaults to ''all''');
    return
end

myargs=textscan(args,'%s');
if length(myargs{1}) < 2
    type = 'all';
else
    type=char(myargs{1}(2));
    if ~strcmp(type,'all') && ~strcmp(type,'backbone')
        add_msg_board(sprintf('ERROR: Unknown r.m.s.d. type %s.',type));
        return
    end
end

adr = char(myargs{1}(1));

rindices1 = resolve_address(adr);
if isempty(rindices1)
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return
end
[m,n] = size(rindices1);
if n == 3
    add_msg_board('ERROR: Chain model selection is not allowed for ensemble r.m.s.d. computation.');
    return
elseif n < 3
    adr = strcat(adr,'{1}'); % select chain models 1 for structure or chain
    rindices1 = resolve_address(adr);
else
    rindices1(:,3) = ones(m,1); % select residues, atoms, or locations in chain models 1
end
if strcmp(type,'backbone')
    if n > 4
        add_msg_board('ERROR: Backbone r.m.s.d. cannot be requested for atoms or locations.');
        return
    else
        adr = strcat(adr,'.N,CA,C');
        rindices1 = resolve_address(adr);
    end
end

[m1,~] = size(rindices1);
if m1 < 1
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return;
end

mods = length(model.structures{rindices1(1,1)}(rindices1(1,2)).residues);
[msg,xyzcell] = get_object(rindices1,'xyz');
if msg.error
    add_msg_board(sprintf('ERROR: Coordinate retrieval failed (%s)',msg.text));
    return
end
xd = 0;
for k = 1:m1
    if iscell(xyzcell)
        xt = xyzcell{k};
    else
        xt = xyzcell;
    end
    [mt,nt] = size(xt);
    xd = xd + mt;
end
xyz = zeros(xd,nt);
xd = 0;
for k = 1:m1
    if iscell(xyzcell)
        xt = xyzcell{k};
    else
        xt = xyzcell;
    end
    [mt,~] = size(xt);
    xyz(xd+1:xd+mt,:) = xt;
    xd = xd + mt;
end
[mc,~] = size(xyz);

for k = 2:mods
    rindices1(:,3) = k*ones(m1,1);
    [~,xyzcell] = get_object(rindices1,'xyz');
    xd = 0;
    for kk = 1:m1
        if iscell(xyzcell)
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end
        [mt,nt] = size(xt);
        xd = xd + mt;
    end
    xyzc = zeros(xd,nt);
    xd = 0;
    for kk = 1:m1
        if iscell(xyzcell)
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end
        [mt,~] = size(xt);
        xyzc(xd+1:xd+mt,:) = xt;
        xd = xd + mt;
    end
    [mcc,~] = size(xyzc);
    if mcc~=mc
        add_msg_board('ERROR: Inconsistent coordinate sets in different models');
        return;
    end
    xyz = xyz + xyzc;
end
xyz = xyz/mods; % mean coordinates
rmsd = 0;
for k = 1:mods
    rindices1(:,3) = k*ones(m1,1);
    [~,xyzcell] = get_object(rindices1,'xyz');
    xd = 0;
    for kk = 1:m1
        if iscell(xyzcell)
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end
        [mt,nt] = size(xt);
        xd = xd + mt;
    end
    xyzc = zeros(xd,nt);
    xd = 0;
    for kk = 1:m1
        if iscell(xyzcell)
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end
        [mt,~] = size(xt);
        xyzc(xd+1:xd+mt,:) = xt;
        xd = xd + mt;
    end
    diff = xyz - xyzc;
    rmsd = rmsd + sum(sum(diff.^2));
end
rmsd = sqrt(rmsd/(m1*mods));
add_msg_board(sprintf('Ensemble r.m.s.d. for selection %s is %4.2f Å',adr,rmsd)); 

function handles = ensemble_reference_rmsd(handles,args)

global model

command=sprintf('refrmsd %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: refrmsd adr refadr [type]');
    add_msg_board('type can be ''all'' or ''backbone'' and defaults to ''all''');
    return
end

myargs=textscan(args,'%s');

if length(myargs{1}) < 2
    add_msg_board('ERROR: At least two address arguments are required.');
    return
end

if length(myargs{1}) < 3
    type = 'all';
else
    type=char(myargs{1}(3));
    if ~strcmp(type,'all') && ~strcmp(type,'backbone')
        add_msg_board(sprintf('ERROR: Unknown r.m.s.d. type %s.',type));
        return
    end
end

adr = char(myargs{1}(1));
refadr = char(myargs{1}(2));

rindices1 = resolve_address(adr);
if isempty(rindices1),
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return
end;
[m,n] = size(rindices1);
if n == 3,
    add_msg_board('ERROR: Chain model selection is not allowed for ensemble r.m.s.d. computation.');
    return
elseif n < 3
    adr = strcat(adr,'{1}'); % select chain models 1 for structure or chain
    rindices1 = resolve_address(adr);
else
    rindices1(:,3) = ones(m,1); % select residues, atoms, or locations in chain models 1
end;
if strcmp(type,'backbone'),
    if n > 4,
        add_msg_board('ERROR: Backbone r.m.s.d. cannot be requested for atoms or locations.');
        return
    else
        adr = strcat(adr,'.N,CA,C');
        rindices1 = resolve_address(adr);
    end;
end;

[m1,~] = size(rindices1);
if m1 < 1,
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return;
end;

rindices2 = resolve_address(refadr);
if isempty(rindices2),
    add_msg_board(sprintf('ERROR: Reference address %s could not be resolved.',refadr));
    return
end;
[~,nr] = size(rindices2);
if strcmp(type,'backbone'),
    if nr > 4,
        add_msg_board('ERROR: Backbone r.m.s.d. cannot be requested for atoms or locations.');
        return
    else
        adr = strcat(refadr,'.N,CA,C');
        rindices2 = resolve_address(adr);
    end;
end;
[mr,~] = size(rindices2);


mods = length(model.structures{rindices1(1,1)}(rindices1(1,2)).residues);
add_msg_board(sprintf('%i models in the ensemble.',mods));
[msg,xyzcell] = get_object(rindices2,'xyz');
if msg.error,
    add_msg_board(sprintf('ERROR: Coordinate retrieval of reference structure failed (%s)',msg.text));
    return
end;
xd = 0;
for k = 1:mr,
    if iscell(xyzcell),
        xt = xyzcell{k};
    else
        xt = xyzcell;
    end;
    [mt,nt] = size(xt);
    xd = xd + mt;
end;
xyz = zeros(xd,nt);
xd = 0;
for k = 1:m1,
    if iscell(xyzcell),
        xt = xyzcell{k};
    else
        xt = xyzcell;
    end;
    [mt,~] = size(xt);
    xyz(xd+1:xd+mt,:) = xt;
    xd = xd + mt;
end;
[mc,~] = size(xyz);

msq = 0;
for k = 1:mods,
    rindices1(:,3) = k*ones(m1,1);
    [~,xyzcell] = get_object(rindices1,'xyz');
    xd = 0;
    for kk = 1:m1,
        if iscell(xyzcell),
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end;
        [mt,nt] = size(xt);
        xd = xd + mt;
    end;
    xyzc = zeros(xd,nt);
    xd = 0;
    for kk = 1:m1,
        if iscell(xyzcell),
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end;
        [mt,~] = size(xt);
        xyzc(xd+1:xd+mt,:) = xt;
        xd = xd + mt;
    end;
    diff = xyz - xyzc;
    rmsd = sqrt(sum(sum(diff.^2))/mc);
%     [rmsd,coor2b,transmat]=rmsd_superimpose(xyz,xyzc);
%     fprintf(1,'RMSD of both structures is: %4.2f Å\n',rmsd);
    msq = msq + rmsd^2;
end;
rmsd = sqrt(msq/mods);
add_msg_board(sprintf('Ensemble r.m.s.d. for selection %s with respect to %s is %4.2f Å',adr,refadr,rmsd)); 
add_msg_board(sprintf('Ensemble rmsd_100 for selection %s with respect to %s is %4.2f Å',adr,refadr,rmsd/(1+log(sqrt(mc/100))))); 

function handles = local_rmsd(handles,args)

global model

command=sprintf('locrmsd %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: locrmsd adr');
    add_msg_board('''adr'' must be a peptide chain address');
    return
end

myargs=textscan(args,'%s');

adr = char(myargs{1}(1));

indices = resolve_address(adr);
if isempty(indices)
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return
end
[~,n] = size(indices);
if n ~= 2
    add_msg_board('ERROR: For local ensemble r.m.s.d. computation, a chain address must be given.');
    add_msg_board('Usage: locrmsd adr');
    add_msg_board('''adr'' must be a peptide chain address');
    return
end

snum = indices(1);
cnum = indices(2);

residues=length(model.structures{snum}(cnum).residues{1}.info);
resnum = zeros(1,residues);
for k = 1:residues
    resnum(k) = model.structures{snum}(cnum).residues{1}.info(k).number;
end

figure(1); clf;
set(gca,'FontSize',12);
plot(resnum,model.structures{snum}(cnum).CA_rmsd,'k.');
xlabel('Residue number');
ylabel('rmsd [Å]');
title(sprintf('Calpha ensemble rmsd for chain %s',adr));
fname = sprintf('ensemble_CA_rmsd_%s.mat',adr);
data = [resnum' model.structures{snum}(cnum).CA_rmsd'];
save(fname,'data');

function handles = local_order(handles,args)

global model

command=sprintf('locorder %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args))
    add_msg_board('Usage: locorder adr');
    add_msg_board('''adr'' must be a peptide or nucleic acid chain address');
    return
end

myargs=textscan(args,'%s');

adr = char(myargs{1}(1));

indices = resolve_address(adr);
if isempty(indices)
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return
end
[~,n] = size(indices);
if n ~= 2
    add_msg_board('ERROR: For local order analysis, a chain address must be given.');
    add_msg_board('Usage: locorder adr');
    add_msg_board('''adr'' must be a peptide or nucleic acid chain address');
    return
end

snum = indices(1);
cnum = indices(2);

if model.structures{snum}(cnum).seqtype == 2
    nucac = true;
elseif model.structures{snum}(cnum).seqtype == 1
    nucac = false;
else
    add_msg_board('ERROR: Chain is neither a peptide nor a nucleic acid chain. Aborting.');
    return
end

if length(model.structures{snum}(cnum).residues) < 2
    add_msg_board('ERROR: Not an ensemble structure. Aborting.');
    return
end

if nucac
    [resaxis,s,npN,npC,meanct,uncert,msg] = analyze_local_order_RNA(indices);
else
    [resaxis,s,npN,npC,meanct,uncert,msg] = analyze_local_order(indices);
end

if msg.error
    add_msg_board(sprintf('ERROR (locorder): %s',msg.text)); 
end

figure(1); clf;
set(gca,'FontSize',12);
errorbar(resaxis,s,uncert(1,:),uncert(2,:),'k.','MarkerSize',14);
xlabel('Residue number');
ylabel('Order parameter s');
title(sprintf('Local order parameter for chain %s',adr));

figure(2); clf;
set(gca,'FontSize',12);
hold on
plot(resaxis,npN,'.','Color',[0,0,0.75],'MarkerSize',14);
plot(resaxis,npC,'.','Color',[0.75,0,0],'MarkerSize',14);
xlabel('Residue number');
ylabel('n_p');
if nucac
    title(sprintf('Shape persistence numbers n_{p,C5''} and n_{p,C3''} for chain %s',adr));
else
    title(sprintf('Shape persistence numbers n_{p,N} and n_{p,C} for chain %s',adr));
end

fname = sprintf('ensemble_local_order_%s.mat',adr);
save(fname,'resaxis','s','npN','npC','meanct','uncert');

function inertiaframe

global model

num_models = length(model.structures{model.current_structure}(1).xyz);
for k = 1:num_models,
    xyz = model.structures{model.current_structure}(1).xyz{k};
    inertia=inertia_tensor(xyz);
    [v,d] = eig(inertia);
    % diag_inert = v'*inertia*v;
    transmat=eye(4,4);
    transmat(1:3,1:3)=v';
    matrices{k} = transmat;
end;

transform_structure_ensemble(model.current_structure,matrices);

% for k = 1:num_models,
%     xyz = model.structures{model.current_structure}(1).xyz{k};
%     inertia=inertia_tensor(xyz);
% end;

function handles = ensemble_rg(handles,args)
% radius of gyration for an ensemble

global model

command=sprintf('radgyr %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: radgyr adr [type]');
    add_msg_board('type can be ''all'' or ''backbone'' and defaults to ''all''');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) < 2,
    type = 'all';
else
    type=char(myargs{1}(2));
    if ~strcmp(type,'all') && ~strcmp(type,'backbone'),
        add_msg_board(sprintf('ERROR: Unknown atom selection type %s.',type));
        return
    end;
end;

adr = char(myargs{1}(1));

rindices1 = resolve_address(adr);
if isempty(rindices1),
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return
end;
[m,n] = size(rindices1);
if n == 3,
    add_msg_board('ERROR: Chain model selection is not allowed for ensemble r.m.s.d. computation.');
    return
elseif n < 3
    adr = strcat(adr,'{1}'); % select chain models 1 for structure or chain
    rindices1 = resolve_address(adr);
else
    rindices1(:,3) = ones(m,1); % select residues, atoms, or locations in chain models 1
end;
if n >= 4,
    add_msg_board('ERROR: Radius of gyration can be computed only for chains or complexes of chains.');
    return
elseif strcmp(type,'backbone'),
    adr = strcat(adr,'.N,CA,C');
    rindices1 = resolve_address(adr);
end;

[m1,~] = size(rindices1);
if m1 < 1,
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return;
end;

mods = length(model.structures{rindices1(1,1)}(rindices1(1,2)).residues);
[msg,xyzcell] = get_object(rindices1,'xyz');
if msg.error,
    add_msg_board(sprintf('ERROR: Coordinate retrieval failed (%s)',msg.text));
    return
end;

rg2 = 0;
natoms = 0;
rgvec = zeros(1,mods);
for k = 1:mods,
    rindices1(:,3) = k*ones(m1,1);
    [~,xyzcell] = get_object(rindices1,'xyz');
    xd = 0;
    for kk = 1:m1,
        if iscell(xyzcell),
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end;
        [mt,nt] = size(xt);
        xd = xd + mt;
    end;
    xyzc = zeros(xd,nt);
    xd = 0;
    for kk = 1:m1,
        if iscell(xyzcell),
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end;
        [mt,~] = size(xt);
        xyzc(xd+1:xd+mt,:) = xt;
        xd = xd + mt;
    end;
    [mcc,~] = size(xyzc);
    if k == 1,
        mc = mcc;
    end;
    if mcc~=mc,
        add_msg_board('ERROR: Inconsistent coordinate sets in different models');
        return;
    end;
    meanc=mean(xyzc,1);
    centered=xyzc-repmat(meanc,mcc,1);
    rgvec(k) = sqrt(sum(sum(centered.^2))/mcc);
    rg2 = rg2 + sum(sum(centered.^2));
    natoms = natoms + mcc;
end;
radgyr = sqrt(rg2/natoms);
add_msg_board(sprintf('Radius of gyration for selection %s is %4.2f Å',adr,radgyr)); 

save radgyr rgvec

function handles = ensemble_compact(handles,args)

global model

command=sprintf('compact %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: compact adr [type]');
    add_msg_board('type can be ''all'' or ''backbone'' and defaults to ''all''');
    return
end;

myargs=textscan(args,'%s');
if length(myargs{1}) < 2,
    type = 'all';
else
    type=char(myargs{1}(2));
    if ~strcmp(type,'all') && ~strcmp(type,'backbone'),
        add_msg_board(sprintf('ERROR: Unknown type %s for compacting.',type));
        return
    end;
end;

adr = char(myargs{1}(1));

rindices1 = resolve_address(adr);
if isempty(rindices1),
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return
end;
[m,n] = size(rindices1);
if n >= 3,
    add_msg_board('ERROR: Chain model, residue or atom selection is not allowed for ensemble compaction.');
    return
elseif n < 3
    adr = strcat(adr,'{1}'); % select chain models 1 for structure or chain
    rindices1 = resolve_address(adr);
else
    rindices1(:,3) = ones(m,1); % select residues, atoms, or locations in chain models 1
end;
if strcmp(type,'backbone'),
    adr = strcat(adr,'.N,CA,C');
    rindices1 = resolve_address(adr);
end;

[m1,~] = size(rindices1);
if m1 < 1,
    add_msg_board(sprintf('ERROR: Address %s could not be resolved.',adr));
    return;
end;

mods = length(model.structures{rindices1(1,1)}(rindices1(1,2)).residues);
[msg,xyzcell] = get_object(rindices1,'xyz');
if msg.error,
    add_msg_board(sprintf('ERROR: Coordinate retrieval failed (%s)',msg.text));
    return
end;
xd = 0;
for k = 1:m1,
    if iscell(xyzcell),
        xt = xyzcell{k};
    else
        xt = xyzcell;
    end;
    [mt,nt] = size(xt);
    xd = xd + mt;
end;
xyz = zeros(xd,nt);
xd = 0;
coor = cell(1,mods);
for k = 1:m1,
    if iscell(xyzcell),
        xt = xyzcell{k};
    else
        xt = xyzcell;
    end;
    [mt,~] = size(xt);
    xyz(xd+1:xd+mt,:) = xt;
    xd = xd + mt;
end;
[mc,~] = size(xyz);
coor{1} = xyz;

for k = 2:mods,
    rindices1(:,3) = k*ones(m1,1);
    [~,xyzcell] = get_object(rindices1,'xyz');
    xd = 0;
    for kk = 1:m1,
        if iscell(xyzcell),
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end;
        [mt,nt] = size(xt);
        xd = xd + mt;
    end;
    xyzc = zeros(xd,nt);
    xd = 0;
    for kk = 1:m1,
        if iscell(xyzcell),
            xt = xyzcell{kk};
        else
            xt = xyzcell;
        end;
        [mt,~] = size(xt);
        xyzc(xd+1:xd+mt,:) = xt;
        xd = xd + mt;
    end;
    [mcc,~] = size(xyzc);
    if mcc~=mc,
        add_msg_board('ERROR: Inconsistent coordinate sets in different models');
        return;
    end;
    coor{k} = xyzc;
end;

msqmat = zeros(mods);
for k=1:mods-1,
    coor1 = coor{k};
    for k2 = k+1:mods,
        coor2 = coor{k2};
        rms=rmsd_superimpose(coor1,coor2);
        msqmat(k,k2)=rms^2;
        msqmat(k2,k)=rms^2;
    end;
end;
rms = 1e6;
maxr = 0;
best = 0;
worst = 0;
for k=1:mods,
    crms=sqrt(sum(msqmat(k,:))/(mods-1));
    if crms < rms,
        rms = crms;
        best = k;
    end;
    if crms > maxr,
        maxr = crms;
        worst = k;
    end;
end;
if best == 0,
    add_msg_board('Compaction requires an ensemble of at least two conformations');
    return
end;
add_msg_board(sprintf('Representative model is no. %i with mean r.m.s.d. of %4.2f Å to all other models.',best,rms));
add_msg_board(sprintf('Most remote model is no. %i with mean r.m.s.d. of %4.2f Å to all other models.',worst,maxr));

xyz = model.structures{model.current_structure}(rindices1(1)).xyz{best};

inertia=inertia_tensor(xyz);
[v,d] = eig(inertia);
transmat=eye(4,4);
transmat(1:3,1:3)=v';

xyz = coor{best};
template = affine_trafo_coor(xyz,transmat);

for k = 1:mods,
    xyz = model.structures{model.current_structure}(rindices1(1)).xyz{k};
    [~,~,transmat]=rmsd_superimpose(template,coor{k});
    model.structures{model.current_structure}(rindices1(1)).xyz{k} = affine_trafo_coor(xyz,transmat);
end;

function handles=mushroom(handles,args)

global rotamer_libraries

command=sprintf('mushroom %s',strtrim(args));
undo_cmd='noundo';
[handles]=cmd_history(handles,command,undo_cmd);

if isempty(args) || isempty(strtrim(args)),
    add_msg_board('Usage: mushroom adr1 adr2 label');
    add_msg_board('or   : mushroom adr label');
    add_msg_board('Only residue addresses are allowed.');
    return
end;

myargs=textscan(args,'%s');
switch length(myargs{1})
    case 2
        adr1=char(myargs{1}(1));
        adr2='';
        label=char(myargs{1}(2));
    case 3
        adr1=char(myargs{1}(1));
        adr2=char(myargs{1}(2));
        label=char(myargs{1}(3));
    otherwise
        add_msg_board('ERROR: Command requires two or three arguments.');
end;



rindices1 = resolve_address(adr1);
[m1,n1] = size(rindices1);
if m1 ~= 1 || n1 ~=4,
    add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr1));
    return;
end;
if ~isempty(adr2)
    rindices2 = resolve_address(adr2);
    [m2,n2] = size(rindices2);
    if m2 ~= 1 || n2 ~=4,
        add_msg_board(sprintf('ERROR: Address %s is not a unique residue specification.',adr2));
        return
    end;
else
    rindices2 = [];
end;

rot_lib_name = '';

for k = 1:length(rotamer_libraries),
    if strcmpi(label,rotamer_libraries(k).tc) || strcmpi(label,rotamer_libraries(k).label),
        id = 0;
        for kk = 1:length(rotamer_libraries(k).T),
            if rotamer_libraries(k).T(kk) == 298,
                id = kk;
            end;
        end;
        if id >0,
            rot_lib_name = id2tag(id,rotamer_libraries(k).files);
        end;
    end;
end;
if isempty(rot_lib_name)
    add_msg_board(sprintf('ERROR: Unknown label %s.',label));
    return
end

load(rot_lib_name);
label1 = rotamer_populations(rindices1,rot_lib);
[nr1,~] = size(label1.NOpos);
[pst0,center0] = point_spread_tensor(label1.NOpos(:,1:3),rot_lib.calibration.pop');
[dircos0,diag0] = eig(pst0);

[pst,center] = point_spread_tensor(label1.NOpos(:,1:3),label1.NOpos(:,4));
[dircos,diag] = eig(pst);
[x,y,z] = ellipsoid(0,0,0,sqrt(2*diag(1,1)),sqrt(2*diag(2,2)),sqrt(2*diag(3,3)),20);
for k1 = 1:21,
    for k2 = 1:21,
        coor = [x(k1,k2) y(k1,k2) z(k1,k2)];
        coor = coor*dircos';
        x(k1,k2) = coor(1)+center(1);
        y(k1,k2) = coor(2)+center(2);
        z(k1,k2) = coor(3)+center(3);
    end;
end;
surf(x,y,z);


if ~isempty(rindices2)
    label2 = rotamer_populations(rindices2,rot_lib);
end

function handles = remodel_section(handles,args)

if isempty(args) || isempty(strtrim(args)) 
    add_msg_board('Usage: remodel adr1i adr1f adr2i adr2f');
    add_msg_board('where adr1i is the address of the initial residue of the target');
    add_msg_board('where adr1f is the address of the final residue of the target');
    add_msg_board('where adr2i is the address of the initial residue of the template');
    add_msg_board('where adr2f is the address of the final residue of the template');
    return
end

command=sprintf('remodel %s',strtrim(args));
[handles,veto]=cmd_history(handles,command);
if veto
    add_msg_board('remodel command cancelled by user');
    return
end

myargs=textscan(args,'%s');

if length(myargs{1}) < 4
    add_msg_board('Usage: remodel adr1i adr1f adr2i adr2f');
    add_msg_board('where adr1i is the address of the initial residue of the target');
    add_msg_board('where adr1f is the address of the final residue of the target');
    add_msg_board('where adr2i is the address of the initial residue of the template');
    add_msg_board('where adr2f is the address of the final residue of the template');
    return
end

adr1i = char(myargs{1}(1));
adr1f = char(myargs{1}(2));
adr2i = char(myargs{1}(3));
adr2f = char(myargs{1}(4));

replace_section(adr1i,adr1f,adr2i,adr2f);

