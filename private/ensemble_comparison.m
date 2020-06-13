function [sigma,pair_rmsd,ordering,densities] = ensemble_comparison(stag1,chains1,pop1,stag2,chains2,pop2,options)
% [sigma,pair_rmsd] = ensemble_comparison(stag1,chains1,stag2,chains2,options)
%
% Analyzes internal variation of conformations in ensemble structures as
% well as variation between two ensembles, this can be used to assess if
% two ensemble structures are more diverse between themselves than
% internally
%
% stag1     structure tag of first ensemble
% chains1   vector of single-character chain identifiers of chains to be
%           included in analysis of the first ensemble, defaults to all
%           chains
% pop1      population vector, if present, must have as many elements as
%           there are ensemble members (conformations), defaults to uniform
%           distribution (for first ensemble)
% stag2     structure tag of second ensemble, if emtpy or missing, only
%           internal analysis of first ensemble is performed
% chains2   vector of single-character chain identifiers of chains to be
%           included in analysis of the second ensemble, defaults to all
%           chains, if present, length must coincide with length of chains1
% pop2      population vector, if present, must have as many elements as
%           there are ensemble members (conformations), defaults to uniform
%           distribution (for second ensemble)
% options   computation options, structure with fields
%           .mode   'backbone' or 'heavy' for backbone or heavy atoms,
%           defaults to 'backbone'
%
% sigma     vector of ensemble standard deviations, the first two elements
%           correspond to ensmebles 1 and 2, element 3 is cross-ensemble
% pair_rmsd cell of matrices [n1,n1], [n2,n2], [n1,n2] of pairwise rmsd of
%           conformers in the two ensembles and between ensembles, where n1
%           and n2 are numbers of conformers in the first and second
%           ensemble, respectively
% ordering  cell of vectors [1,n1], [1,n2], [1,n1], [1,n2] that assign
%           original conformer numbers to entries in the sorted matrices
% densities ensemble densities and similarities, single number, if only one
%           ensemble is input, otherwise vector(1,4) density1, density2,
%           similarity12, similarity21
%
% G. Jeschke, 25.9.2019-04.03.2020

global model

hfig = gcf;
set(hfig,'Pointer','watch');
drawnow

if ~exist('options','var') || isempty(options)
    options.mode = 'backbone';
end

if ~isfield(options,'cluster_sorting')
    options.cluster_sorting = false;
end

if exist('stag2','var') && ~isempty(stag2)
    sigma = zeros(1,3);
    pair_rmsd = cell(1,3);
    options.ensembles = 2;
    [sind2,msg] = resolve_address(sprintf('[%s]',stag2));
    if isempty(sind2) || msg.error
        add_msg_board('ERROR: Second ensemble structure does not exist. Aborting.');
        add_msg_board(sprintf('%s',msg.text));
        set(hfig,'Pointer','arrow');
        return
    end
else
    sigma = 0;
    pair_rmsd = {[]};
    options.ensembles = 1;
end

[sind1,msg] = resolve_address(sprintf('[%s]',stag1));
if isempty(sind1) || msg.error
    add_msg_board('ERROR: First ensemble structure does not exist. Aborting.');
    add_msg_board(sprintf('%s',msg.text));
    set(hfig,'Pointer','arrow');
    return
end

if ~exist('chains1','var') || isempty(chains1)
    chains1 = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';
    pch = 0;
    for kc = 1:length(model.structures{sind1})
        if model.structures{sind1}(kc).seqtype > 0
            pch = pch + 1;
            chains1(pch) = model.structures{sind1}(kc).name(1);
        end
    end
    chains1 = chains1(1:pch);
end

if options.ensembles > 1
    if ~exist('chains2','var') || isempty(chains2)
        chains2 = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';
        for kc = 1:length(model.structures{sind2})
            chains2(kc) = model.structures{sind2}(kc).name(1);
        end
        chains2 = chains2(1:length(model.structures{sind2}));
    end
end

if options.ensembles > 1 && length(chains1) ~= length(chains2)
    add_msg_board('ERROR: Mismatch in the (specified) numbers of chains for the two ensembles. Aborting.');
    set(hfig,'Pointer','arrow');
    return
end

adr0 = sprintf('[%s](%s)',stag1,chains1(1));
[cindices,msg] = resolve_address(adr0);
if msg.error
    msg.text = sprintf('Address %s cannot be resolved (%s)',adr0,msg.text);
    set(hfig,'Pointer','arrow');
    return
end
nconf1 = length(model.structures{cindices(1)}(cindices(2)).xyz); % number of conformers
if ~exist('pop1','var') || isempty(pop1)
    pop1 = ones(1,nconf1);
end
pop1 = pop1/sum(pop1);

if options.ensembles > 1
    adr0 = sprintf('[%s](%s)',stag2,chains2(1));
    [cindices,msg] = resolve_address(adr0);
    if msg.error
        msg.text = sprintf('Address %s cannot be resolved (%s)',adr0,msg.text);
        set(hfig,'Pointer','arrow');
        return
    end
    nconf2 = length(model.structures{cindices(1)}(cindices(2)).xyz); % number of conformers
    if ~exist('pop2','var') || isempty(pop2)
        pop2 = ones(1,nconf2);
    end
    pop2 = pop2/sum(pop2);
end

if length(pop1) ~= nconf1
    add_msg_board('ERROR: Mismatch in the numbers of conformers and length of population vector for ensemble 1. Aborting.');
    set(hfig,'Pointer','arrow');
    return
end

if options.ensembles > 1 && length(pop2) ~= nconf2
    add_msg_board('ERROR: Mismatch in the numbers of conformers and length of population vector for ensemble 2. Aborting.');
    set(hfig,'Pointer','arrow');
    return
end

options.mode = 'backbone';
[pair_rmsd1,msg] = ensemble_pair_rmsd(stag1,chains1,options);
if msg.error
    add_msg_board('ERROR: Pair rmsd of ensemble 1 could not be computed. Aborting.');
    add_msg_board(sprintf('%s',msg.text));
    set(hfig,'Pointer','arrow');
    return
end
rmsdsum = sum(pair_rmsd1);
[~,sorting] = sort(rmsdsum);
pair_rmsd{1} = pair_rmsd1(sorting,sorting);
ordering{1} = sorting;

sigma(1) = get_sigma(pair_rmsd{1},pop1);

if options.ensembles == 1 && options.cluster_sorting
    options2.fname = sprintf('cluster_analysis_ensemble_%s.dat',stag1);
    [new_pair_rmsd,new_ordering] = cluster_analysis(pair_rmsd{1},pop1,sorting,options2);
    ordering{1} = new_ordering;
    pair_rmsd{1}= new_pair_rmsd;
end

if options.ensembles > 1
    [pair_rmsd2,msg] = ensemble_pair_rmsd(stag2,chains2,options);
    if msg.error
        add_msg_board('ERROR: Pair rmsd of ensemble 2 could not be computed. Aborting.');
        add_msg_board(sprintf('%s',msg.text));
        set(hfig,'Pointer','arrow');
        return
    end
    rmsdsum = sum(pair_rmsd2);
    [~,sorting] = sort(rmsdsum);
    pair_rmsd{2} = pair_rmsd2(sorting,sorting);
    ordering{2} = sorting;
    sigma(2) = get_sigma(pair_rmsd{2},pop2);
    
    [cross_rmsd,msg] = cross_ensemble_pair_rmsd(stag1,chains1,stag2,chains2,options);
    if msg.error
        add_msg_board('ERROR: Cross-ensemble pair rmsd could not be computed. Aborting.');
        add_msg_board(sprintf('%s',msg.text));
        set(hfig,'Pointer','arrow');
        return
    end
    rmsdsum1 = sum(cross_rmsd,1);
    rmsdsum2 = sum(cross_rmsd,2);
    [~,sorting1] = sort(rmsdsum1);
    ordering{3} = sorting1;
    [~,sorting2] = sort(rmsdsum2);
    ordering{4} = sorting2;
    pair_rmsd{3} = cross_rmsd(sorting2,sorting1);
    sigma(3) = get_sigma(pair_rmsd{3},pop1,pop2);
end

% compute ensemble densities
densities = zeros(1,4);
for ke = 1:options.ensembles
    distmat = pair_rmsd{ke};
    [m,~] = size(distmat);
    msq = zeros(1,m);
    for k = 1:m
        cmsd = distmat(k,:).^2;
        cmsd(k) = 1e12;
        msq(k) = min(cmsd);
    end
    densities(ke) = sqrt(mean(msq));
end
if options.ensembles < 2
    densities = densities(1);
else
    distmat = pair_rmsd{3};
    [m1,m2] = size(distmat);
    msq = zeros(1,m1);
    for k = 1:m1
        cmsd = distmat(k,:).^2;
        msq(k) = min(cmsd);
    end
    densities(3) = sqrt(mean(msq));
    msq = zeros(1,m2);
    for k = 1:m2
        cmsd = distmat(:,k).^2;
        msq(k) = min(cmsd);
    end
    densities(4) = sqrt(mean(msq));
end

set(hfig,'Pointer','arrow');

function sigma = get_sigma(pair_rmsd,pop1,pop2)

sigma = [];
[m,n] = size(pair_rmsd);

if ~exist('pop2','var') || isempty(pop2) % internal variation
    if m ~= n
        add_msg_board('ERROR: Requested internal variation for non-square matrix.');
        return
    end
    if m ~= length(pop1)
        add_msg_board('ERROR: Population vector length does not match rmsd matrix size.');
        return
    end
    msdsum = 0;
    popsum = 0;
    for k1 = 1:m-1
        for k2 = k1+1:m
            msdsum = msdsum + pop1(k1)*pop1(k2)*pair_rmsd(k1,k2)^2;
            popsum = popsum + pop1(k1)*pop1(k2);
        end
    end
    sigma = sqrt(msdsum/popsum);
else
    if m ~= length(pop1) || n ~= length(pop2)
        add_msg_board('ERROR: Population vector lengths do not match rmsd matrix size.');
        return
    end
    msdsum = 0;
    popsum = 0;
    for k1 = 1:m
        for k2 = 1:n
            msdsum = msdsum + pop1(k1)*pop2(k2)*pair_rmsd(k1,k2)^2;
            popsum = popsum + pop1(k1)*pop2(k2);
        end
    end
    sigma = sqrt(msdsum/popsum);
end

