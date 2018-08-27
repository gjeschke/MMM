function [defs,links,sites,failed] = rd_definition_stemloop_library(fname)
% function [defs,links,sites,failed] = rd_definition_stemloop_library(fname)
%
% reads definitions for an RNA stemloop library based on Rosetta (Rosie
% FarFar) or other models
%
% there are almost no syntax checks, as this has to be done very rarely
%
% format of the definition file
%
% - comments start with %
% - key # STEMLIB nti nte rrm motif
%   nti    number of the first nt in a stemloop model
%   nte    number of the last nt in a stemloop model
%   rrm    chain identifier of the RRM to which the RNA binds
%   motif  chain identifier of the RNA binding motif
% - # STEMLIB is a block key with following lines corresponding to model
%   sets
% - model set lines have the following format
%   directory idef edef iused eused ibind ebind
%   directory   name of the folder containing models as PDB
%   idef, edef  first and last nt defined in a model (numbering as in the 
%               model PDB), 
%   isued,eused first and last nt contained in the library conformations 
%               (numbering as in the model PDB)
%   ibind,ebind first and last nt that bind to the protein
%   opt         the optional argument 'opt' requests AMBER optimization of
%               the binding pose, this does not currently work well 
% 
% - key # RLINK, block key specifying links between RNA pieces
%       block key with lines of the form
%   anchi   initial anchor nt
%   anche   final anchor nt
%
% - key # SITES, block key specifying spin-labelled sites, lines of the
%               form
%   adr     address of the spin-labelled site
%   label   label type
%
% fname     file name of the definition file
%
% defs      definitions, array of structs
% links     list of RNA-RNA links
% sites     list of spin-labelled sites
% failed    flag that indicates failure
%
% G. Jeschke, 19.3.2018

failed = true;

fid=fopen(fname);
if fid==-1
    add_msg_board('ERROR: Restraint file does not exist');
    return;
end

nl = 0;
defpoi = 0;
linkpoi = 0;
sitepoi = 0;
mode = 0;
while 1
    tline = fgetl(fid);
    nl = nl + 1;
    if ~ischar(tline) || mode < 0, break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            if k(1)>1
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end
        end
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#')
            switch upper(char(args(2)))
                case 'STEMLIB'
                    mode=1;
                    if defpoi > 0
                        defs(defpoi).offset = defs(defpoi).offset(1:dirpoi);
                        defs(defpoi).defined = defs(defpoi).defined(1:dirpoi,:);
                        defs(defpoi).used = defs(defpoi).used(1:dirpoi,:);
                        defs(defpoi).binding = defs(defpoi).binding(1:dirpoi,:);
                        defs(defpoi).optimize = defs(defpoi).optimize(1:dirpoi);
                    end
                    defpoi = defpoi + 1;
                    defs(defpoi).name = '';
                    defs(defpoi).range = [];
                    defs(defpoi).protein = '';
                    defs(defpoi).motif = '';
                    defs(defpoi).directory{1} = '';
                    defs(defpoi).defined = zeros(100,2);
                    defs(defpoi).used = zeros(100,2);
                    defs(defpoi).binding = zeros(100,2);
                    defs(defpoi).offset = zeros(1,100);
                    defs(defpoi).optimize = zeros(1,100);
                    dirpoi = 0;
                    irange = str2double(char(args(4)));
                    erange = str2double(char(args(5)));
                    defs(defpoi).name = char(args(3));
                    defs(defpoi).range = [irange,erange];
                    defs(defpoi).protein = char(args(6));
                    defs(defpoi).motif = char(args(7));
                case 'RLINK'
                    mode = 2;
                case 'SITES'
                    mode = 3;
                otherwise
                    mode=0;
                    add_msg_board('Warning: Unknown definition mode');
            end
        elseif mode>0 && ~strncmp(strtrim(char(args(1))),'%',1) % lines inside the block key
            switch mode
                case 1 % STEMLIB
                    dirpoi = dirpoi + 1;
                    defs(defpoi).directory{dirpoi} = char(args(1));
                    idef = str2double(char(args(2)));
                    edef = str2double(char(args(3)));
                    iused = str2double(char(args(4)));
                    eused = str2double(char(args(5)));
                    ibind = str2double(char(args(6)));
                    ebind = str2double(char(args(7)));
                    defs(defpoi).offset(dirpoi) = defs(defpoi).range(1) - iused;
                    defs(defpoi).defined(dirpoi,:) = [idef,edef];
                    defs(defpoi).used(dirpoi,:) = [iused,eused];
                    defs(defpoi).binding(dirpoi,:) = [ibind,ebind];
                    if length(args) > 7 && strcmpi(char(args(8)),'opt')
                        defs(defpoi).optimize(dirpoi) = true;
                    else
                        defs(defpoi).optimize(dirpoi) = false;
                    end
                case 2 % RLINK
                    linkpoi = linkpoi + 1;
                    links(linkpoi).adr1 = char(args(1));
                    links(linkpoi).adr2 = char(args(2));
                case 3 % SITES
                    sitepoi = sitepoi + 1;
                    sites(sitepoi).adr = char(args(1));
                    sites(sitepoi).label = char(args(2));
            end
        end
    end
end
if defpoi > 0 && dirpoi > 0
    failed = false;
    defs(defpoi).offset = defs(defpoi).offset(1:dirpoi);
    defs(defpoi).defined = defs(defpoi).defined(1:dirpoi,:);
    defs(defpoi).used = defs(defpoi).used(1:dirpoi,:);
    defs(defpoi).binding = defs(defpoi).binding(1:dirpoi,:);
    defs(defpoi).optimize = defs(defpoi).optimize(1:dirpoi);
end

fclose(fid);

