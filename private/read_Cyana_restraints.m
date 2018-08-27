function restab = read_Cyana_restraints(fname,trans)
% restab = read_Cyana_restraints(fname,trans)
% 
% reads upper or lower bound restraints in CYANA format with optional
% translation of residue numbers to MMM addresses and optional checking of
% consistency
%
% fname     name of the restraint file with extension, caller is
%           responsible for interpreting restraint type (upper or lower
%           bound)
% trans     optional translation table of the form
%           trans(k).range = [ra,re]    ra, re first and last residue
%                                       number in a chain section
%           trans(k).chain              chain identifier (string) for
%                                       this range, MMM address form
%                                       [stag](cid), where stag is the
%                                       structure tag and cid the chain tag
%           trans(k).offset             residue offset between restraint
%                                       file and structure file
%           if a translation table exists, then all residue numbers
%           in the restraint file must fall into exactly onbe of the
%           specified ranges
%
% restab    restraint table in the form
%           restab(k).adr1  address of first residue
%           restab(k).adr2  address of second residue
%           restab(k).r     distance restraint [Å]
%           restab(k).type1 residue type (e.g. TYR, C) of first residue
%           restab(k).type2 residue type of second residue
%
% G. Jeschke, 16.8.2017

if ~exist('trans','var') || isempty(trans)
    translation = false;
else
    translation = true;
end

if ~exist('check','var') || isempty(check)
    check = false;
end;

restab = [];

fid=fopen(fname);
if fid==-1,
    fprintf(2,'ERROR: Restraint file does not exist');
    return;
end;

nl = 0;
resnum = 0;
while 1
    tline = fgetl(fid);
    nl = nl + 1;
    if ~ischar(tline) break, end
    if ~isempty(tline),
        k = strfind(tline,'#'); % remove comments
        if ~isempty(k),
            if k(1)>1,
                tline = tline(1:k(1)-1);
            else
                tline = '';
            end;
        end;
        if ~isempty(tline)
            myline = textscan(tline,'%s');
            args = myline{1};
            narg = length(args);
            if narg < 7
                fprintf(2,'Warning: Wrong number of arguments (%i) in line %i of the restraint file (7 expected).\n',narg,nl);
            else
                resnum = resnum + 1;
                restab(resnum).r = str2double(char(args(7)));
                restab(resnum).type1 = Cyana2PDB(char(args(2)));
                restab(resnum).type2 = Cyana2PDB(char(args(5)));
                r1 = str2double(char(args(1)));
                if translation
                    chain = 0;
                    for k = 1:length(trans)
                        if r1 >= trans(k).range(1) && r1 <= trans(k).range(2)
                            if chain > 0
                                fprintf(2,'ERROR: A residue appears in two ranges (%i,%i) within the translation table. Aborting.\n',chain,k);
                                fclose(fid);
                                return
                            else
                                chain = k;
                                adstub1 = trans(k).chain;
                                r1 = r1 + trans(k).offset;
                            end
                        end
                    end
                    if chain == 0
                        fprintf(2,'Warning: Residue %i not found in any range within the translation table. Restraint ignored.\n',r1);
                    end
                else
                    adstub1 = '';
                end
                restab(resnum).adr1 = sprintf('%s%i.%s',adstub1,r1,char(args(3)));
                r2 = str2double(char(args(4)));
                if translation
                    chain = 0;
                    for k = 1:length(trans)
                        if r2 >= trans(k).range(1) && r2 <= trans(k).range(2)
                            if chain > 0
                                fprintf(2,'ERROR: A residue appears in two ranges (%i,%i) within the translation table. Aborting.\n',chain,k);
                                fclose(fid);
                                restab = [];
                                return
                            else
                                chain = k;
                                adstub2 = trans(k).chain;
                                r2 = r2 + trans(k).offset;
                            end
                        end
                    end
                    if chain == 0
                        fprintf(2,'Warning: Residue %i not found in any range within the translation table. Restraint ignored.\n',r2);
                    end
                else
                    adstub2 = '';
                end
                restab(resnum).adr2 = sprintf('%s%i.%s',adstub2,r2,char(args(6)));
            end
        end
    end
end

fclose(fid);

function type = Cyana2PDB(type)

switch type % translate CYANA nucleotide names to PDB convention
    case 'RADE'
        type = '  A';
    case 'RCYT'
        type = '  C';
    case 'RGUA'
        type = '  G';
    case 'URA'
        type = '  U';
    case 'ADE'
        type = ' DA';
    case 'CYT'
        type = ' DC';
    case 'GUA'
        type = ' DG';
    case 'THY'
        type = ' DT';
end
