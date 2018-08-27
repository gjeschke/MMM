function type = selection_type(indices)
% type = selection_type(indices)
%
% Determines the type of selected chains, chain models, or residues
%
% indices   MMM index array m-by-n, n = 2 chains selected, n = 3 chain
%           models selected, n = 4 residues selected, otherwise result is
%           'undetermined', n = 1 structure selected
%
% type      'peptide' for all chains being peptide chains or all residues
%           being peptide residues, HET peptides are allowed
%           'nucleotide' for all chains/residues being nucleotide residues
%           'mixed' for all chains/residues being either peptide or
%           nucleotide residues
%           'nonstandard' for at least one chain not being a peptide or 
%           nucleic acid chain or at least one residue being neither 
%           peptide nor nucleotide
%           'undetermined' for atoms or locations being selected
%
% G. Jeschke, 2014

global model

type = 'undetermined';

[m,n] = size(indices);

if n > 4, % atoms or atom locations
    return;
end;

pep_found = false;
nuc_found = false;
nstd_found = false;

for k = 1:m, % all selected objects
    switch n,
        case 1 % structures
            for kk = 1:length(model.structures{indices(k,1)}),
                ctype = get_chain_type(model,[indices(k,1),kk]);
                switch ctype
                    case {-1,0}
                        nstd_found = true;
                    case 1
                        pep_found = true;
                    case {2,3}
                        nuc_found = true;
                    otherwise
                        nstd_found = true;
                end;
            end;
        case 2 % chains
            ctype = get_chain_type(model,indices(k,:));
            switch ctype
                case {-1,0}
                    nstd_found = true;
                case 1
                    pep_found = true;
                case {2,3}
                    nuc_found = true;
                otherwise
                    nstd_found = true;
            end;
        case 3 % chain models
            ctype = get_chain_type(model,indices(k,1:2));
            switch ctype
                case {-1,0}
                    nstd_found = true;
                case 1
                    pep_found = true;
                case {2,3}
                    nuc_found = true;
                otherwise
                    nstd_found = true;
            end;
        case 4 % residues
            rtype = get_residue_type(model,indices(k,:));
            switch rtype
                case 0
                    nstd_found = true;
                case 1
                    pep_found = true;
                case 2
                    nuc_found = true;
                otherwise
                    nstd_found = true;
            end;
    end;
end;

if pep_found && ~nuc_found && ~nstd_found,
    type = 'peptide';
end;
if nuc_found && ~pep_found && ~nstd_found,
    type = 'nucleotide';
end;
if nuc_found && pep_found,
    type = 'mixed';
end;
if nuc_found && nstd_found,
    type = 'mixed';
end;
if pep_found && nstd_found,
    type = 'mixed';
end;
if nstd_found,
    type = 'nonstandard';
end;

% for testing
% fprintf(1,'Selection type is: %s\n',type);

function ctype = get_chain_type(model,indices)
% Determines the type of a chain
%
% ctype     -1 undetermined sequence
%           0  Heterosequence
%           1  peptide sequence
%           2  DNA sequence
%           3  RNA sequence

snum = indices(1);
cnum = indices(2);
sequence = model.structures{snum}(cnum).sequence;
if isfield(model.structures{snum}(cnum),'seqtype'), % determine sequence type, information is present
    ctype = model.structures{snum}(cnum).seqtype;
    seqtest=double(sequence)>=double('a');
    if ctype == 2 && sum(seqtest),
        ctype=3;
    end;
end;

function rtype = get_residue_type(model,indices)
% Dtermines type of residue
%
% rtype     0 non-standard
%           1 peptide
%           2 nucleotide

rtype = model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).type;
