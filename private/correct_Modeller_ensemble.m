function correct_Modeller_ensemble(snum,ares,chain_tags)
% corrects residue numbers and chain tags in a structure ensemble (or
% single structure) that was created by Modeller and read via
% rd_pdb_ensemble.m
%
% snum          number of the structure to be corrected
% ares          (optional) start residue, if specified and not empty, this
%               must be either the integer residue number of the first
%               residue in the first chain or a vector of first residue
%               numbers for all chains, empty ares or missing ares defaults
%               to residue numer 1 for the first residue of all chains
% chain_tags    (optional) list of chain tags, defaults to :A:B:C: and so
%               on, default is also used if argument is present, but empty
%
% Remarks:
% - returns without action if structure snum does not exist
%

global model

if snum<1 || snum>length(model.structures),
    return
end;

if nargin<2,
    ares=[];    
end;

chains=length(model.chain_ids{snum});

if nargin<3 || isempty(chain_tags),
    chain_tags=':';
    for k=1:chains,
        chain_tags=[chain_tags char(double('A')-1+k) ':'];
    end;
end;

model.chain_tags{snum}=chain_tags;
currchain=id2tag(1,chain_tags);
model.current_chain=currchain;

for kc=1:chains,
    offset=0;
    if length(ares)>=kc,
        offset=ares(kc)-1;
    end;
    ctag=id2tag(kc,chain_tags);
    model.structures{snum}(kc).name=ctag;
    seq=model.structures{snum}(kc).sequence;
    poi=offset;
    while poi>0,
        poi=poi-1;
        seq=['?' seq];
    end;
    model.structures{snum}(kc).sequence=seq;
    restags=model.structures{snum}(kc).restags;
    poi=offset;
    while poi>0,
        poi=poi-1;
        restags=[':UNK' restags];
    end;
    model.structures{snum}(kc).restags=restags;
    restags=':';
    for km=1:length(model.structures{snum}(kc).residues),
        info=model.structures{snum}(kc).residues{km}.info;
        coffset=offset-info(1).number+1;
        for kr=1:length(info),
            info(kr).number=info(kr).number+coffset;
            restags=sprintf('%s%i:',restags,info(kr).number);
        end;
        model.structures{snum}(kc).residues{km}.info=info;
        model.structures{snum}(kc).residues{km}.residue_tags=restags;
    end;
    if offset~=0,
        if isfield(model.structures{snum}(kc),'helix_defs'),
            for kh=1:length(model.structures{snum}(kc).helix_defs),
                model.structures{snum}(kc).helix_defs{kh}.range=...
                    model.structures{snum}(kc).helix_defs{kh}.range+offset;
            end;
        end;
        if isfield(model.structures{snum}(kc),'loop_defs'),
            for kh=1:length(model.structures{snum}(kc).loop_defs),
                model.structures{snum}(kc).loop_defs{kh}.range=...
                    model.structures{snum}(kc).loop_defs{kh}.range+offset;
            end;
        end;
        if isfield(model.structures{snum}(kc),'sheet_defs'),
            for kh=1:length(model.structures{snum}(kc).sheet_defs),
                model.structures{snum}(kc).sheet_defs{kh}.range=...
                    model.structures{snum}(kc).sheet_defs{kh}.range+offset;
            end;
        end;
    end;
end;