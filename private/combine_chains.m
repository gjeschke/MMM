function [chain,warnings] = combine_chains(indices,ctag,dbref,header)
% chain = combine_chains(indices)
%
% Combines several chains, which may stem from different structures, 
% into a single chain
% residue numbering in the combined chain is the same as in the input
% chains, the caller is responsible for it not to overlap
% if not all models of an input chain have the same number of atoms,
% connection information will be wrong and an error is set
% the .seqres field is set empty, as it is not used by MMM in writing out
% the structure or for internal purposes
% the .fasta field is set empty as it may lead to inconsistencies for
% chimera chains
% the .seqexist field is not set (it is not set for newly read structures)
% 
% indices   cell array of chain indices, each element must be a vector of
%           either 2 or 3 indices
%           for 2 indices, all models of the chain are copied
%           for 3 indices, the given model is combined with all models of
%           other chains
%           all chains with 2 indices must have the same number of models
%           or exactly one model
% ctag      chain tag (single uppercase letter) for the combined chain
% dbref     optional data base reference for the sequence of the whole
%           chain, chain.dbref will be empty if the parameter is missing
% header    optional header for the sequence in FASTA format, initialized
%           to unknown if it is not specified or empty
%
% chain     chain record for an MMM model; model.structures{s}(c)
%           if an error occurs, chain.name reports the error text and
%           chain.xyz is empty;
% warnings  cell array of warning messages 
%
% G. Jeschke, 17.5.2017

global model

warnings = {};

if~exist('dbref','var')
    dbref = '';
end

if~exist('header','var') || isempty(header)
    header = 'unknown';
end

chain.xyz = cell(1,0);
maxat = 10000;

% Determine whether the number of chain models is consistent between chains

nm = 0; % number of models
for kc = 1:length(indices)
    cind = indices{kc};
    if length(cind) == 3
        cnm = 1; % a single model is specified
        coor = model.structures{cind(1)}(cind(2)).xyz{cind(3)};
        [na,~] = size(coor);
        maxat = maxat + na;
    else
        cnm = length(model.structures{cind(1)}(cind(2)).xyz); % there are as many models as coordinate arrays
        coor = model.structures{cind(1)}(cind(2)).xyz{1};
        [na,~] = size(coor);
        maxat = maxat + na;
    end
    if nm == 0
        nm = cnm;
    else
        if cnm > 1
            if nm == 1
                nm = cnm;
            else
                if nm ~= cnm
                    chain.name = 'Mismatch in number of models per chain';
                    return
                end
            end
        end
    end
end

% initialize fields of the chain structure

chain.name = ctag;
chain.modified = 0;
chain.nonstandard = 0;
chain.helices = 0;
chain.strands = 0;
chain.xyz = cell(1,nm);
chain.isotopes = [];
chain.Bfactor = cell(1,nm);
chain.Btensor = cell(1,nm);
chain.atoms = cell(1,nm);
for km = 1:nm
    chain.atoms{km} = 0;
end;
chain.conn = zeros(maxat,30,'int32');
chain.maxconn = 0;
chain.het = [];
het = 0;
helices = 0;
strands = 0;
chain.sequence = '';
chain.seqres = [];
chain.seqtype = -1;
chain.fasta = [];
chain.restags = ':';
chain.dbref = dbref;
chain.header = header;
chain.loop_defs = {};
loop_defs = 0;
chain.helix_defs = {};
helix_defs = 0;
chain.sheet_defs = {};
sheet_defs = 0;
chain.seqexist = [];
chain.residues = cell(1,nm);

warnpoi = 0;

% redistribute and reindex the data

isotopes = zeros(maxat,2,'single');
cat0 = zeros(1,length(indices));
for km = 1:nm % all output models
    atpoi = 0;
    respoi = 0;
    absresnum = 0;
    xyz = zeros(maxat,3);
    Bfactor = zeros(1,maxat);
    Btensor = zeros(maxat,6,'int32');
    residues.residue_tags = ':';
    for kc = 1:length(indices) % all input chains
        cind = indices{kc};
        cchain = model.structures{cind(1)}(cind(2));
        % determine what chain model should be used
        if length(cind) == 2
            cnm = length(cchain.xyz);
            if cnm == 1
                cind = [cind 1];
            else
                cind = [cind km];
            end;
        end;
        iso = cchain.isotopes;
        [cat,~] = size(iso);
        if km == 1
            cat0(kc) = cat;
            chain.modified = chain.modified + cchain.modified;
            chain.nonstandard = chain.nonstandard + cchain.nonstandard;
            chain.helices = chain.helices + cchain.helices;
            chain.strands = chain.strands + cchain.strands;
            isotopes(atpoi+1:atpoi+cat,:) = cchain.isotopes;
            [~,clen] = size(cchain.conn);
            if clen > chain.maxconn
                chain.maxconn = clen;
            end
            chain.conn(atpoi+1:atpoi+cat,1:clen) = cchain.conn;
            chain.seqtype = cchain.seqtype;
        else
            if cat~=cat0(kc)
                warnpoi = warnpoi + 1;
                warnings{warnpoi} = sprintf('Different numbers of atoms in models of input chain %i',kc);
            end
            if chain.seqtype~=cchain.seqtype
                warnpoi = warnpoi + 1;
                warnings{warnpoi} = sprintf('Sequence type mismatch of input chain %i with previous chains',kc);
            end
        end
        xyz(atpoi+1:atpoi+cat,:) = cchain.xyz{cind(3)};
        Bfactor(atpoi+1:atpoi+cat) = cchain.Bfactor{cind(3)};
        Btensor(atpoi+1:atpoi+cat,:) = cchain.Btensor{cind(3)};
        chain.atoms{km} = chain.atoms{km} + cat;
        if km == 1
            if isfield(cchain,'het')
                chain.het(het+1:het+length(cchain.het)) = cchain.het;
                het = het + length(cchain.het);
                for kh = het+1:het+length(cchain.het)
                    chain.het(kh).number = chain.het(kh).number + atpoi;
                end
            end
            if isfield(cchain,'helix')
                if ~isempty(cchain.helix)
                    chain.helix(helices+1:helices+length(cchain.helix)) = cchain.helix;
                end
                for kh = helices+1:helices+length(cchain.helix)
                    chain.helix(kh).start = chain.helix(kh).start + respoi;
                    chain.helix(kh).end = chain.helix(kh).end + respoi;
                end
                helices = helices + length(cchain.helices);
            end
            if isfield(cchain,'strand')
                if ~isempty(cchain.strand)
                    chain.strand(strands+1:strands+length(cchain.strand)) = cchain.strand;
                end
                for kh = strands+1:strands+length(cchain.strand)
                    chain.strand(kh).start = chain.strand(kh).start + respoi;
                    chain.strand(kh).end = chain.strand(kh).end + respoi;
                end
                strands = strands + length(cchain.strands);
            end
            chain.sequence = strcat(chain.sequence,restag2seq(cchain.restags));
            chain.restags = strcat(chain.restags,cchain.restags(2:end));
            if isfield(cchain,'loop_defs')
                chain.loop_defs(loop_defs+1:loop_defs+length(cchain.loop_defs)) = cchain.loop_defs;
                for kh = loop_defs+1:loop_defs+length(cchain.loop_defs)
                    chain.loop_defs{kh}.range = chain.loop_defs{kh}.range + respoi;
                end
                loop_defs = loop_defs + length(cchain.loop_defs);
            end
            if isfield(cchain,'helix_defs')
                chain.helix_defs(helix_defs+1:helix_defs+length(cchain.helix_defs)) = cchain.helix_defs;
                for kh = helix_defs+1:helix_defs+length(cchain.helix_defs)
                    chain.helix_defs{kh}.range = chain.helix_defs{kh}.range + respoi;
                end
                helix_defs = helix_defs + length(cchain.helix_defs);
            end
            if isfield(cchain,'sheet_defs')
                chain.sheet_defs(sheet_defs+1:sheet_defs+length(cchain.sheet_defs)) = cchain.sheet_defs;
                for kh = sheet_defs+1:sheet_defs+length(cchain.sheet_defs)
                    chain.sheet_defs{kh}.range = chain.sheet_defs{kh}.range + respoi;
                end
                sheet_defs = sheet_defs + length(cchain.sheet_defs);
            end
        end
        residues.residue_tags = strcat(residues.residue_tags,cchain.residues{cind(3)}.residue_tags(2:end));
        % correct atom numbers
        info = cchain.residues{cind(3)}.info;
        for kr = 1:length(info)
            ratnum = info(kr).atom_numbers;
            absresnum = absresnum + 1;
            info(kr).absresnum = absresnum;
            for ka = 1:length(ratnum)
                aatnum = ratnum{ka};
                [m,n] = size(aatnum);
                if n == 3
                    for kl = 1:m
                        aatnum(kl,1) = aatnum(kl,1) + atpoi;
                    end
                else
                    aatnum = aatnum + atpoi;
                end
                ratnum{ka} = aatnum;
            end
            info(kr).atom_numbers = ratnum;
            if kc ~= length(indices)
                info(kr).terminal = [];
            end
        end;
        if kc == 1
            residues.info = info;
        else
            inames = fieldnames(residues.info(1));
            clear infoc
            for kf = 1:length(inames)
                for kr = 1:length(info)
                    if isfield(info(kr),inames{kf})
                        infoc(kr).(inames{kf}) = info(kr).(inames{kf});
                    else
                        infoc(kr).(inames{kf}) = [];
                    end
                end
            end
            residues.info = [residues.info infoc];
        end
        atpoi = atpoi + cat;
        respoi = respoi + length(cchain.sequence);
    end
    chain.xyz{km} = xyz(1:atpoi,:);
    chain.Bfactor{km} = Bfactor(1:atpoi);
    chain.Btensor{km} = Btensor(1:atpoi,:);
    chain.residues{km} = residues;
end
chain.isotopes = isotopes(1:atpoi,:);
chain.conn = chain.conn(1:atpoi,1:chain.maxconn);