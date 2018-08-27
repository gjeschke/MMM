function chain=mk_internal_bonds(chain,resnum,residue_defs,non_standard)
% function chain=mk_internal_bonds(chain,resnum,residue_defs)
%
% makes intramolecular bonds of standard residues from template, very
% defensive programming, slow
% can also be applied to non-standard residues, although these should
% usually have CONECT records, in this case flag non_standard must be
% present and true, bonds for non-standard residues are only made, if these
% residues are defined in the PDB monomer directory or an equivalent file
% is in the definition directory
%
% chain         Matlab structure for a single chain as defined in rd_pdb
% resnum        internal residue number
% residue_defs  standard residue definitions
% non_standard  optional flag that requests bonds for a non-standard
%               residue, if set (1), if zero, no effect
%
% G. Jeschke, 2009

if nargin<4, % if non-standard flag is missing, only standard residues are treated 
    non_standard=0;
end;

if chain.residues{1}.info(resnum).connected, % if it's already done, return immediately
    return
end;

template={};
info=chain.residues{1}.info;
myname=strtrim(info(resnum).name);
res_id=tag2id(myname,upper(residue_defs.restags)); % check if this is an amino acid
if isempty(res_id), % check if this is something else (including nucleic acids)
    if non_standard % if requested, get template for non-standard residue
        template=rd_pdb_template(myname);
    end;
else % get amino acid template
    template=residue_defs.residues(res_id); % template from preloaded residue definitions
end;
if ~isempty(template) % make intramolecular connections if template exists
    atoms=textscan(info(resnum).atom_tags(2:end),'%s','Delimiter',':');
    for k3=1:length(atoms{1}),
        myatom=char(atoms{1}(k3));
        myatom_id=tag2id(myatom,info(resnum).atom_tags);
        myatom_num=info(resnum).atom_numbers{myatom_id};
        template_id=tag2id(myatom,template.atoms);
        if ~isempty(template_id) && ~isempty(template.conn),
            bonds=template.conn(template_id,:);
            for k4=1:length(bonds),
                if bonds(k4)>0,
                    partner_tag=id2tag(bonds(k4),template.atoms);
                    partner=tag2id(partner_tag,info(resnum).atom_tags);
                    if ~isempty(partner),
                        partner_num=info(resnum).atom_numbers{partner};
                        [m1,n1]=size(myatom_num); % necessary because of alternate positions
                        [m2,n2]=size(partner_num); % necessary because of alternate positions
                        for o=1:m1,
                            for oo=1:m2,
                                if n1<3 || n2<3 || myatom_num(o,3)==partner_num(oo,3) % if at least one atom has no alternate locations
                                                                                      % or both alternate location tags are the same
                                    chain=make_bond(chain,myatom_num(o,1),partner_num(oo,1));
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;
chain.residues{1}.info(resnum).connected=1; % store information that it is done
