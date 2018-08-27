function [fname,access]=accessibility_report(newrad,indices)

global model
global chemistry
global general
global hMain
global help_files
global third_party
global residue_defs

tolerance=1.2; % tolerance for surface distance in Å

fname='';

entry=strcat(help_files,'third_party.html#MSMS');

dospath=which('msms.exe');
if isempty(dospath),
    add_msg_board('Accessibility analysis requires MSMS by Michel Sanner');
    add_msg_board('ERROR: MSMS could not be found on the Matlab path.');
    add_msg_board('Please check whether MSMS is installed and the path set.');
    add_msg_board('(see also help browser)');
    webcall(entry,'-helpbrowser');
    fname='';
    return
end;

add_msg_board('Computing solvent-accessible surface...');
adr=[mk_address(model.current_structure) '(:){1}'];
sindices=resolve_address(adr);

[msg,coor]=get_object(sindices,'xyz');
if iscell(coor),
    coor=cat(1,coor{:});
end;

[msg,elements]=get_object(sindices,'elements');
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


if nargin<1,
    radius=1.5;
    title=sprintf('Accessibility of residues in %s',adr);
    quest='Radius of the probe sphere in Å:';
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(quest,title,1,{sprintf('%6.1f',radius)},options);
    newrad=str2double(answer);
    if isnan(newrad),
        newrad=radius;
    end;
end;
radstr=sprintf(' -probe_radius %4.1f',newrad);

set(hMain.MMM,'Pointer','watch');
drawnow;
outfile=[general.tmp_files 'tmp.xyzr'];
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

msmsfile=[general.tmp_files 'tmp'];
[pathstr, name, ext] = fileparts(msmsfile);
msmsfile=fullfile(pathstr,[name ext]);
cmd=[dospath ' -if ' outfile ' -of ' msmsfile dstring radstr];
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
    add_msg_board('Now importing solvent accessible surface...');
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

[tri,x,y,z,info]=rd_msms(msmsfile);


[mt,nt]=size(tri);

add_msg_board(sprintf('Solvent accessible surface with %i faces imported.',mt));
add_msg_board('Now computing face areas and midpoints...');

midpoints=zeros(mt,3);
area=zeros(1,mt);
for k=1:mt,
    vert=[x(tri(k,:))',y(tri(k,:))',z(tri(k,:))'];
    a=norm(vert(1,:)-vert(2,:));
    b=norm(vert(1,:)-vert(3,:));
    c=norm(vert(2,:)-vert(3,:));
    midpoints(k,:)=mean(vert);
    s=(a+b+c)/2;
    area(k)=sqrt(s*(s-a)*(s-b)*(s-c));
end;

add_msg_board('Face areas and midpoints computed.');
add_msg_board('Now scanning residues...');

if nargin<2,
    indices=resolve_address('*');
    if isempty(indices),
        snum=model.current_structure;
        ctag=model.current_chain;
        cnum=tag2id(ctag,model.chain_tags{snum});
        indices0=[snum cnum 1];
        [msg,indices]=get_object(indices0,'children');
        if iscell(indices),
            indices=vertcat(indices{:});
        end;
    end;

    if isempty(indices),
        add_msg_board('No residues addressed');
    end;
end;
[m,n]=size(indices);

fname=[general.tmp_files 'accessibility.dat'];
fid=fopen(fname,'w');
fprintf(fid,'%% MMM accessibility analysis based on MSMS by M. Sanner\n');
fprintf(fid,'%% (negative rel. accessibility means that total surface area is not known)\n');
fprintf(fid,'%% Acc. area (Å^2)   rel. acc.  %% residue address\n');

access=zeros(5000,6);
poi=0;
for k=1:m,
    comp=round(100*(k-1)/m);
    if mod(k,20)==0,
        add_msg_board(sprintf('%i%% scanned.',comp));
    end;
    cindices=indices(k,:);
    cindices=cindices(cindices>0);
    if length(cindices)==1,
        sadr=mk_address(cindices);
        adr=sprintf('%s(:){1}',sadr);
        cindices=resolve_address(adr);
        [msg,cindices]=get_object(cindices,'children');
        if iscell(cindices),
            cindices=vertcat(cindices{:});
        end;
    elseif length(cindices)==2,
        cadr=mk_address(cindices);
        adr=sprintf('%s{1}',cadr);
        cindices=resolve_address(adr);
        [msg,cindices]=get_object(cindices,'children');
        if iscell(cindices),
            cindices=vertcat(cindices{:});
        end;
    elseif length(cindices)==3,
        [msg,cindices]=get_object(cindices,'children');
        if iscell(cindices),
            cindices=vertcat(cindices{:});
        end;
    end;
    if ~isempty(cindices),
        [mm,n]=size(cindices);
        for kk=1:mm,
            myindices=cindices(kk,:);
            [msg,xyz]=get_residue(myindices,'xyz_paradigm');
            [msg,elements]=get_residue(myindices,'elements_paradigm');
            [message,tag]=get_residue(myindices,'name');
            if strcmpi(tag,'R1A'),
                tag='ILE';
            end;
            resnum=tag2id(tag,upper(residue_defs.restags));
            if ~isempty(resnum),
                info=residue_defs.residues(resnum);
            else
                info=[];
            end;
            mask=zeros(1,mt);
            for kkk=1:length(elements),
                if elements(kkk)>0 && elements(kkk)<length(chemistry.pse),
                    vdW=chemistry.pse(elements(kkk)).vdW;
                    rad=(vdW+tolerance)^2;
                    atom=repmat(xyz(kkk,:),mt,1);
                    dvec=midpoints-atom;
                    dist=sum(dvec.^2,2)-rad;
                    dist=dist<=0;
                    mask=mask+dist';
                end;
            end;
            myarea=sum(area(mask>0));
            if ~isempty(info),
                relacc=myarea/info.A;
            else
                relacc=-1;
            end;
            adr=mk_address(myindices);
            fprintf(fid,'%6.1f%19.3f      %% %s\n',myarea,relacc,adr);
            poi=poi+1;
            access(poi,1:4)=myindices;
            access(poi,5)=myarea;
            access(poi,6)=relacc;
        end;
    end;
end;

fclose(fid);
access=access(1:poi,:);

set(hMain.MMM,'Pointer','arrow');
