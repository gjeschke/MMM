function template=rd_pdb_template(fname)
% function template=rd_pdb_template(fname)
%
% Reads a single residue/molecule template in PDB monomer library format
% data are stored in a structure with name template
%
% the function first looks for a local file, otherwise it tries to
% establish an FTP connection and to download the file from the PDB web site
% an empty structure is returned, if none of these methods succeeds
%
% fname         filename for the template, if no extension is given, files
%               with no extension and extension .txt are tried in that
%               sequence
% template      structure containing fields:
%               resname: the PDB residue name
%               atoms: a string of colon-separated atom names, their
%                      sequence defines the numbering of the atoms in conn
%                      and elements
%               conn:  connection table, 
%               name:  the full name of the residue (lower case)
%               elements: element numbers of all atoms

global web_adr
global queries

PDB_ftp=web_adr.PDB_US;
monomer_lib=queries.PDB_monomers;

whereami=which('rd_pdb_template.m');
mypath=whereami(1:end-length('rd_pdb_template.m'));

defpath=[mypath '../definitions/'];

fname1=strcat(defpath,fname);

template={};

fid=fopen(fname1);
if fid==-1, % local file with given extension or without extension cannot be opened
    fname2=strcat(fname1,'.txt');
    fid=fopen(fname2);
    if fid==-1, % local file with added extension '.txt' cannot be opened either
        backdir=pwd;
        cd(defpath);
        ftp_obj=ftp(PDB_ftp);
        cd(ftp_obj,monomer_lib);
        ascii(ftp_obj);
        try
            f1=mget(ftp_obj,fname);
        catch
            disp(lasterror);
            try
                f2=mget(ftp_obj,strcat(fname,'.txt'));
            catch
                disp(lasterror);
            end
        end
        close(ftp_obj);
        cd(backdir);
        fid=fopen(fname1);
        if fid==-1,
            fid=fopen(fname2);
        end;
    end;
end;

if fid~=-1,
    template=analyze_pdb_template(fid);
    fclose(fid);
end;

function template=analyze_pdb_template(fid)

pse=' HHeLiBe B C N O FNeNaMgAlSi P SClAr KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMdNoLr';

template={};
atpoi=0;
atoms=':';
resname='';
hetname='';
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    [mytok,args]=strtok(tline);
    switch mytok
        case 'RESIDUE'
            [resname,arg2]=strtok(args);
            template.resname=resname;
            atnum=str2double(arg2);
            connects=zeros(atnum,12);
        case 'CONECT'
            atpoi=atpoi+1;
            [atname,arg2]=strtok(args);
            atoms=strcat(atoms,atname);
            atoms=strcat(atoms,':');
            connstr{atpoi}=arg2;
        case 'HETNAM'
            [arg2,arg3]=strtok(args);
            if strcmp(arg2,resname), hetname=lower(fliplr(deblank(fliplr(arg3))));
            else
                [arg2,arg3]=strtok(arg3);
                hetname=strcat(hetname,lower(arg3));
            end;
            template.name=hetname;
        case 'FORMUL'
            [arg2,formula]=strtok(args);
    end;
end
if length(atoms)==1,
    elm=[':' resname ':'];
    numelm=1;
    atoms=elm;
    template.atoms=atoms;
else
    template.atoms=atoms;
    % Elemental analysis according to formula
    composition=textscan(formula,'%s');
    [numelm,n]=size(composition{:}); % number of elements
    elm=':';
    for k=1:numelm,
        id=char(composition{1}(k));
        symbol_length=sum(isletter(id));
        elm=strcat(elm,id(1:symbol_length));
        elm=strcat(elm,':');
    end;
    elm=elm(1:end-1);
end;
% Determine elements
atsymb=textscan(atoms(2:end),'%s','Delimiter',':');
[numat,n]=size(atsymb{:}); % number of elements
template.elements=zeros(1,atpoi);
for k=1:numat;
    id=char(atsymb{1}(k));
    element=[':' id(1) ':'];
    k1=findstr(upper(element),upper(elm));
    k2=[];
    if length(id)>1,
        element2=[':' id(1:2) ':'];
        k2=findstr(upper(element2),upper(elm));
    end;
    if ~isempty(k1) && ~isempty(k2)
        disp(sprintf('Element assignment is ambiguous for symbol %s',id'));
    else
        if ~isempty(k2),
            element=element2;
        end;
        element=element(2:end-1);
        if length(element)<2,
            element=[' ' element];
        end;
        if strcmp(element,' X')
            template.elements(k)=0;
        else
            atnum=findstr(upper(pse),element)+1;
            if length(atnum)~=1,
                atnum=0;
            end;
            template.elements(k)=atnum/2;
        end;
    end;
end;
% Make connections
maxbonds=0;
for k=1:atpoi,
    conny=char(connstr{k});
    [arg1,arg2]=strtok(conny);
    bonds=str2double(arg1);
    if bonds>maxbonds, maxbonds=bonds; end;
    conns=textscan(arg2,'%s');
    for kk=1:bonds,
        partner=[':'  char(conns{1}(kk)) ':'];
        connects(k,kk)=1+sum(find(atoms==':')<findstr(atoms,partner));
    end;
end;  
template.conn=connects(1:atpoi,1:maxbonds);

