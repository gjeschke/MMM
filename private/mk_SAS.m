function msg=mk_SAS(address,tag,radius)
% function msg=mk_SAS(address,tag,radius)
%
% create and display solvent accessible surface for the object addressed by
% address (must be a single object)
% 
% tag       tag by which the surface can be addressed later
% radius    (optional) probe radius, defdaults to 1.5 Å
% 
% this is just a wrapper for M. Sanner's MSMS program
%


global graph_settings
global model
global chemistry
global general
global hMain
global third_party

if nargin<3,
    radius=1.5;
end;

msg.error=0;
msg.text='No error.';

dospath=which('msms.exe');
if isempty(dospath),
    msg.text='ERROR: MSMS could not be found on the Matlab path.';
    msg.error=1;
    return
end;

indices=resolve_address(address);
indices=indices(indices>0);
[m,n]=size(indices);
if m==1,
    if n<=4,
        if n<2,
            adr=[mk_address(indices(1)) '(:){1}'];
        elseif n<3,
            adr=[mk_address(indices(1:2)) '{1}'];
        else
            adr=mk_address(indices);
        end;
        indices=resolve_address(adr);
    else
        msg.text='ERROR: Solvent accessible surface is not defined on atom level.';
        msg.error=2;
        return
    end;
else
    msg.text='ERROR: Exactly one object must be selected for computation of solvent accessible surface.';
    msg.error=3;
    return
end;

[msg1,coor]=get_object(indices,'xyz');
if iscell(coor),
    coor=cat(1,coor{:});
end;

[msg1,elements]=get_object(indices,'elements');
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
dstring=sprintf(' -density %i -probe_radius %4.2f',density,radius);

if isfield(model,'SAS'),
    poi=length(model.SAS)+1;
else
    poi=1;
end;

stag=strtrim(char(tag));
found=0;
if poi>1,
    for k=1:poi-1,
        if strcmp(stag,model.SAS(k).tag),
            found=1;
        end;
    end;
end;
if found,
    msg.text='ERROR: Surface tag already exists.';
    msg.error=4;
    return
end;

set(gcf,'Pointer','watch');
drawnow;
outfile=[general.tmp_files stag '.xyzr'];
fid=fopen(outfile,'w');
if fid==-1,
    msg.error=5;
    msg.text='ERROR: Coordinate file could not be opened for writing.';
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

[pathstr, name, ext, versn] = fileparts(outfile);
outfile=fullfile(pathstr,[name ext versn]);

msmsfile=[general.tmp_files stag];
[pathstr, name, ext, versn] = fileparts(msmsfile);
msmsfile=fullfile(pathstr,[name ext versn]);
cmd=[dospath ' -if ' outfile ' -of ' msmsfile dstring];
[s, w] = dos(cmd);
if s~=0,
    msg.text='ERROR: MSMS did not run successfully.';
    msg.error=6;
    set(hMain.MMM,'Pointer','arrow');
    return
else
    comments=textscan(w,'%s','Delimiter','\n');
    lines=comments{1};
    for k=1:length(lines),
        msg1=char(lines(k));
        add_msg_board(msg1);
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

add_msg_board('Now displaying surface.');

axes(hMain.axes_model);

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

set(gcf,'Pointer','arrow');


