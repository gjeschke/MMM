function compile_amino_acids
% Compiles a definition file for all amino acids with the same numbering
% convention as in the Bioinformatics toolbox of matlab
% templates are taken from the PDB monomer library, special treatment for
% single letter code X, as XAA is defined as a different residue, X is
% encoded as backbone only, i.e. as GLY
%

global rootdir

defpath='/definitions/';
defpath=strcat(rootdir,defpath);
outfile=strcat(defpath,'residues.def');
outfile2=strcat(defpath,'residues.mat');

psefile=strcat(defpath,'PSE.mat');

load(psefile);

betaprop_file=strcat(defpath,'beta_barrel_propensities.dat');
betaprop=load(betaprop_file);

% load the additional definitions from the book D. Whitford, Proteins- Structure and Function,
% Wiley, Chichester, 2005

resdefs=zeros(50,13);
rfile=fopen(strcat(defpath,'amino_acids_Whitford.txt'),'rt');
comm='%';
while comm=='%',
    tline = fgetl(rfile);
    if ~ischar(tline), break, end
    if isempty(tline), continue; end;
    tline=strtrim(tline);
    comm=tline(1);
end;
poi=0;
while 1
    [res,rem]=strtok(tline);
    poi=poi+1;
    reslist{poi}=strtrim(res);
    resdefs(poi,:)=str2num(rem);
    tline = fgetl(rfile);
    if ~ischar(tline), break, end
    if isempty(tline), continue; end;
end;
fclose(rfile);
resdefs=resdefs(1:poi,:);

% definition of amino acids that are classified as heteroatoms in the PDB,
% although they are linked within chains
hettags=':R1A:IA1:CSE:MSE:R1B:R1F:R7A:V1A:MTN:'; 

residue_list=strcat(defpath,'amino_acids_toolbox.txt');
fid0=fopen(residue_list,'rt');
wfile=fopen(outfile,'wt');
fprintf(wfile,'%% Do not edit!\n%% This file exists only for information and as a template for additional definitions.\n');
fprintf(wfile,'%% Actual standard residue definitions are in the binary file residues.mat\n');
restags=':';
single_letter_code='';
while 1
    tline = fgetl(fid0);
    if ~ischar(tline), break, end
    if isempty(tline), continue; end;
    residue=textscan(tline,'%s');
    sc=char(residue{1}(1));
    single_letter_code=[single_letter_code sc];
    num=str2double(residue{1}(2));
    tc=char(residue{1}(3));
    addinfo=[];
    for k=1:length(reslist),
        if strcmpi(reslist{k},tc),
            addinfo=resdefs(k,:);
        end;
    end;
    restags=[restags tc ':'];
    disp(sprintf('Compiling residue %s.',tc));
    name=char(residue{1}(4));
    ncod=length(residue{1});
    codons='';
    for k=5:ncod,
        codons=[codons char(residue{1}(k)) ' '];
    end;
    if char(residue{1}(1))=='*',
        break,
    end;
    if char(residue{1}(1))=='X',
        monomer='GLY';
    else
        monomer=upper(tc);
    end;
    template=rd_pdb_template(monomer);
    template.sc=sc;
    template.tc=tc;
    template.name=name;
    template.codons=codons;
    if ~isempty(addinfo),
        template.pK1=addinfo(1);
        template.pK2=addinfo(2);
        if addinfo(3)>0,
            template.pKr=addinfo(3);
        else
            template.pKr=[];
        end;
        template.V=addinfo(4);
        template.A=addinfo(5);
        template.mass=addinfo(6);
        template.type=addinfo(7);
        template.helix_propensity=addinfo(8);
        template.hydropathy=addinfo(9);
        template.charge=addinfo(10);
        template.hydrophobicity=addinfo(11);
        template.TMLIPH=addinfo(12);
        template.TMLIPC=addinfo(13);
        template.TMBH=-betaprop(num,1);
        template.TMBC=-betaprop(num,2);
    else
        template.pK1=[];
        template.pK2=[];
        template.pKr=[];
        template.V=[];
        template.A=[];
        template.mass=[];
        template.type=[];
        template.helix_propensity=[];        
        template.hydropathy=[];
        template.charge=[];
        template.hydrophobicity=[];
        template.TMLIPH=[];
        template.TMLIPC=[];
        template.TMBH=[];
        template.TMBC=[];
    end;
    residues(num)=template;
    fprintf(wfile,'begin %s\n',template.tc);
    fprintf(wfile,'\t sc %s\n',template.sc);
    fprintf(wfile,'\t name %s\n',template.name);
    fprintf(wfile,'\t codons %s\n',template.codons);
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
fclose(fid0);
fclose(wfile);

residue_defs.residues=residues;
residue_defs.restags=restags;
residue_defs.hettags=hettags;
residue_defs.single_letter_code=single_letter_code;
residue_defs.nucleotide_tags=': DG: DC: DT: DA:  G:  C:  U:  A:';
residue_defs.cyana_nucleotide_tags=':GUA :CYT :THY :ADE :RGUA:RCYT:URA :RADE:';
residue_defs.nucleotide_slc='GCTAgcua';
residue_defs.nucleotides(1).slc='G';
residue_defs.nucleotides(1).backbone=':C4'':C3'':';
residue_defs.nucleotides(1).Hbonds=':O6:N1:N2:';
residue_defs.nucleotides(1).polygon=':N9:C4:N3:C2:N1:C6:C5:N7:C8:';
residue_defs.nucleotides(1).color='lightgreen';
residue_defs.nucleotides(2).slc='C';
residue_defs.nucleotides(2).backbone=':C4'':C3'':';
residue_defs.nucleotides(2).Hbonds=':O2:N3:N4:';
residue_defs.nucleotides(2).polygon=':N1:C2:N3:C4:C5:C6:';
residue_defs.nucleotides(2).color='lightcoral';
residue_defs.nucleotides(3).slc='T';
residue_defs.nucleotides(3).backbone=':C4'':C3'':';
residue_defs.nucleotides(3).Hbonds=':O4:N3:';
residue_defs.nucleotides(3).polygon=':N1:C2:N3:C4:C5:C6:';
residue_defs.nucleotides(3).color='gold';
residue_defs.nucleotides(4).slc='A';
residue_defs.nucleotides(4).backbone=':C4'':C3'':';
residue_defs.nucleotides(4).Hbonds=':N1:N6:';
residue_defs.nucleotides(4).polygon=':N9:C4:N3:C2:N1:C6:C5:N7:C8:';
residue_defs.nucleotides(4).color='mediumpurple';
residue_defs.nucleotides(5).slc='g';
residue_defs.nucleotides(5).backbone=':C4'':C3'':';
residue_defs.nucleotides(5).Hbonds=':O6:N1:N2:';
residue_defs.nucleotides(5).polygon=':N9:C4:N3:C2:N1:C6:C5:N7:C8:';
residue_defs.nucleotides(5).color='lightgreen';
residue_defs.nucleotides(6).slc='c';
residue_defs.nucleotides(6).backbone=':C4'':C3'':';
residue_defs.nucleotides(6).Hbonds=':O2:N3:N4:';
residue_defs.nucleotides(6).polygon=':N1:C2:N3:C4:C5:C6:';
residue_defs.nucleotides(6).color='lightcoral';
residue_defs.nucleotides(7).slc='u';
residue_defs.nucleotides(7).backbone=':C4'':C3'':';
residue_defs.nucleotides(7).Hbonds=':O4:N3:';
residue_defs.nucleotides(7).polygon=':N1:C2:N3:C4:C5:C6:';
residue_defs.nucleotides(7).color='gold';
residue_defs.nucleotides(8).slc='a';
residue_defs.nucleotides(8).backbone=':C4'':C3'':';
residue_defs.nucleotides(8).Hbonds=':N1:N6:';
residue_defs.nucleotides(8).polygon=':N9:C4:N3:C2:N1:C6:C5:N7:C8:';
residue_defs.nucleotides(8).color='mediumpurple';


save(outfile2,'residue_defs');

% Keskin et al. contact potentials for interfaces

data=load('contact_potentials.txt');
translate=data(end,:); % sequence of MMM code numbers corresponding to amonio acid type sequence in Table
solvent_mediated=zeros(20,20);
residue_mediated=zeros(20,20);
for k1=1:20,
    row0=translate(k1+1);
    for k2=1:21,
        if k2>k1,
            col=translate(k2);
            if col<row0,
                row=col;
                col=row0;
            else
                row=row0;
            end;
            solvent_mediated(row,col)=data(k1,k2);
        else
            col=translate(k2+1);
            if col<row0,
                row=col;
                col=row0;
            else
                row=row0;
            end;
            residue_mediated(row,col)=data(k1,k2);
        end;
    end;
end;

% Miyazawa et al. contact potentials

data=load('contact_energies_1996.txt');
translate=data(end,:); % sequence of MMM code numbers corresponding to amonio acid type sequence in Table
contact_energies=zeros(20,20);
for k1=1:20,
    row0=translate(k1);
    for k2=1:20,
        if k2>=k1,
            col=translate(k2);
            if col<row0,
                row=col;
                col=row0;
            else
                row=row0;
            end;
            contact_energies(row,col)=data(k1,k2);
        end;
    end;
end;
save contact_potentials contact_energies solvent_mediated residue_mediated 

