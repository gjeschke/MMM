function get_SAS_chain(chain,modnum)
% hObject    handle to menu_build_SAS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global graph_settings
global model
global chemistry
global general
global hMain
global third_party

if ~exist('modnum','var') || isempty(modnum)
    modnum = 1;
end

dospath=which('msms.exe');
if isempty(dospath)
    add_msg_board('ERROR: MSMS could not be found on the Matlab path.');
    add_msg_board('Please check whether MSMS is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    return
end

coor = chain.xyz{modnum};

indices=resolve_address('*');
indices=indices(indices>0);
[m,n]=size(indices);
if m==1,
    if n<=4,
        if n<2,
            adr=[mk_address(indices(1)) '(:){1}'];
            add_msg_board('Considering first coordinate set of each chain');
        elseif n<3,
            adr=[mk_address(indices(1:2)) '{1}'];
            add_msg_board('Considering first coordinate set of the chain');
        else
            adr=mk_address(indices);
        end;
        indices=resolve_address(adr);
    else
        add_msg_board('Selection addresses an atom or location');
        add_msg_board('Computing surface for current structure');
        adr=[mk_address(model.current_structure) '(:){1}'];
        indices=resolve_address(adr);
    end;
else
    add_msg_board('Selection does not address a single object');
    add_msg_board('Computing surface for current structure');
    adr=[mk_address(model.current_structure) '(:){1}'];
    indices=resolve_address(adr);
end;
[msg,coor]=get_object(indices,'xyz');
if iscell(coor),
    coor=cat(1,coor{:});
end;

% A piece of code that may be useful for coarse-grained visualization
% [k,v]=convhulln(coor);
% add_msg_board(sprintf('Volume of convex hull: %6.1f',v));
% h=trisurf(k,coor(:,1),coor(:,2),coor(:,3));
% set(h, 'FaceColor', [0,0,1], 'EdgeColor', 'none', 'FaceAlpha',0.75,'FaceLighting','gouraud','Clipping','off');
% set(h, 'CDataMapping','direct','AlphaDataMapping','none');


[msg,elements]=get_object(indices,'elements');
if iscell(elements),
    elements=cat(1,elements{:});
end;
vdW=zeros(size(elements));
for k=1:length(elements),
    if elements(k)>0 && elements(k)<length(chemistry.pse),
        vdW(k)=chemistry.pse(elements(k)).vdW;
    end;
end;

density=round(5000/length(elements));
if density>3, density=3; end;
if density<1, density=1; end;
dstring=sprintf(' -density %i',density);

if isfield(model,'SAS'),
    poi=length(model.SAS)+1;
else
    poi=1;
end;

answer = inputdlg('Please provide a tag for this surface',sprintf('Solvent accessible surface for %s',adr),1,{sprintf('SAS_%i',poi)});
stag=strtrim(char(answer));
found=0;
if poi>1,
    for k=1:poi-1,
        if strcmp(stag,model.SAS(k).tag),
            found=1;
        end;
    end;
end;
if found,
    add_msg_board('ERROR: Surface tag already existed.');
    add_msg_board('Please provide a new (unique) tag.');
    return
end;
if isempty(answer),
    add_msg_board('Surface computation canceled (no tag).');
    return
end;

set(hMain.MMM,'Pointer','watch');
drawnow;
outfile=[general.tmp_files stag '.xyzr'];
fid=fopen(outfile,'w');
if fid==-1,
    add_msg_board('ERROR: Coordinate file could not be opened for writing.');
    return
end;
for k=1:length(vdW),
    if vdW(k)>0,
        fprintf(fid,'%12.3f%12.3f%12.3f%6.2f',coor(k,1),coor(k,2),coor(k,3),vdW(k));
        if k<length(vdW),
            fprintf(fid,'\n');
        end;
    end;
end;
fclose(fid);

[pathstr, name, ext] = fileparts(outfile);
outfile=fullfile(pathstr,[name ext]);

msmsfile=[general.tmp_files stag];
[pathstr, name, ext] = fileparts(msmsfile);
msmsfile=fullfile(pathstr,[name ext]);
cmd=[dospath ' -if ' outfile ' -of ' msmsfile dstring];
[s, w] = dos(cmd);
if s~=0,
    add_msg_board('ERROR: MSMS did not run successfully.');
    set(hMain.MMM,'Pointer','arrow');
    return
else
    comments=textscan(w,'%s','Delimiter','\n');
    lines=comments{1};
    for k=1:length(lines),
        msg=char(lines(k));
        add_msg_board(msg);
    end;
end;

% add the reference, if it does not yet exist
msms_ref=true;
id=tag2id('Sanner:1996_msms',third_party.tags,[],'|');
if isfield(model,'auto_references'),
    if ~isempty(find(id==model.auto_references, 1)),
        msms_ref=false;
    end;
else
    model.auto_references=[];
end;
if msms_ref,
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

axes(handles.axes_model);

[tri,x,y,z,info]=rd_msms(msmsfile);
obj = trisurf(tri,x,y,z);
set(obj, 'FaceColor', graph_settings.SAS_color, 'EdgeColor', 'none', 'FaceAlpha',graph_settings.SAS_transparency,'FaceLighting','gouraud','Clipping','off');
set(obj, 'CDataMapping','direct','AlphaDataMapping','none');

dg.gobjects=obj;
dg.tag=['SAS:' stag];
dg.color=graph_settings.SAS_color;
dg.level=info.probe_r;
dg.transparency=graph_settings.SAS_transparency;
dg.active=true;

if isfield(model,'surfaces') && ~isempty(model.surfaces),
    model.surfaces(end+1)=dg;
else
    if isfield(model,'surfaces'),
        model=rmfield(model,'surfaces');
    end;
    model.surfaces(1)=dg;
end;

if ~isfield(model,'SAS')
    model.SAS(1).tri=tri;
    model.SAS(1).x=x;
    model.SAS(1).y=y;
    model.SAS(1).z=z;
    model.SAS(1).info=info;
    model.SAS(1).tag=stag;
    model.SAS(1).snum=indices(1);
else
    poi=length(model.SAS)+1;
    model.SAS(poi).tri=tri;
    model.SAS(poi).x=x;
    model.SAS(poi).y=y;
    model.SAS(poi).z=z;
    model.SAS(poi).info=info;
    model.SAS(poi).tag=stag;
    model.SAS(poi).snum=indices(1);
end;

set(hMain.MMM,'Pointer','arrow');


guidata(hObject,handles);
