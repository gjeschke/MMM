function compile_labels
% Compiles a definition file for spin labels
% templates are taken from the PDB monomer library, or are self-defined
% templates in PDB style format
%

global rootdir

defpath='definitions/';
defpath=strcat(rootdir,defpath);
outfile=strcat(defpath,'labels.def');
outfile2=strcat(defpath,'labels.mat');

psefile=strcat(defpath,'PSE.mat');

load(psefile);

% definition of amino acids that are classified as heteroatoms in the PDB,
% although they are linked within chains
restags=':MTN:R1A:R1B:R1F:R7A:IA1:'; 
names=':MTSL:MTSL:Me-MTSL:Ph-MTSL:Br-MTSL:IA-Proxyl:';
Ntags=':N1:N1:N1:N1:N1:N1:'; % tag of nitroxide N atom
Otags=':O1:O1:O1:O1:O1:O1:'; % tag of nitroxide O atom
Ctags=':C2:C2:C2:C2:C2:C2:'; % tag of C atom bound to nitroxide N, defines xy plane of label frame

colors=[204,0,204;204,0,204;204,0,102;102,0,204;153,0,102;255,0,255;153,51,255];
radius=0.64; % half the N-O bond length
[ncol,k]=size(colors);

wfile=fopen(outfile,'wt');
nonsense=textscan(restags,'%s','Delimiter',':');
residue_names=nonsense{1};

for k=2:length(residue_names),
    tc=char(residue_names(k));
    disp(sprintf('Compiling label %s.',tc));
    id=tag2id(tc,restags);
    name=id2tag(id,names);
    Ntag=id2tag(id,Ntags);
    Otag=id2tag(id,Otags);
    Ctag=id2tag(id,Ctags);    
    monomer=upper(tc);
    template=rd_pdb_template(monomer);
    template.tc=tc;
    template.short_name=name;
    template.frame=sprintf(':%s:%s:%s:',Ntag,Otag,Ctag);
    cpoi=mod(id,ncol)+1;
    template.color=colors(cpoi,:);
    template.radius=radius;
    residues(id)=template;
    fprintf(wfile,'begin %s\n',template.tc);
    fprintf(wfile,'\t name %s\n',template.name);
    fprintf(wfile,'\t short %s\n',template.short_name);
    fprintf(wfile,'\t frame %s %s %s\n',Ntag,Otag,Ctag);
    fprintf(wfile,'\t color %i %i %i\n',template.color);
    fprintf(wfile,'\t radius %5.3f\n',template.radius);
    tags=textscan(template.atoms,'%s','Delimiter',':');
    for k=1:length(template.elements),
        tag=char(tags{1}(k+1));
        if template.elements(k)==0 || template.elements(k)>92,
            symb='X';
        else
            symb=pse(template.elements(k)).symbol;
        end;
        if strcmp(symb,'X'),
            iso=0;
        else
            [abundance,poi]=max(pse(template.elements(k)).isotopes(:,2));
            iso=pse(template.elements(k)).isotopes(poi,1);
        end;
        fprintf(wfile,'\t atom %i %s %s %i\n',k,tag,symb,iso);
    end;
    for k=1:length(template.elements),
        fprintf(wfile,'\t bonds %i ',k);
        conn=template.conn(k,:);
        for kk=1:length(conn),
            if conn(kk)>0,
                fprintf(wfile,'%i ',conn(kk));
            end;
        end;
        fprintf(wfile,'\n');
    end;
    fprintf(wfile,'end %s\n',template.tc);
end;
fclose(wfile);

label_defs.residues=residues;
label_defs.restags=restags;

save(outfile2,'label_defs');

