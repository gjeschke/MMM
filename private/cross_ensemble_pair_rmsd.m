function [pair_rmsd,msg] = cross_ensemble_pair_rmsd(stag1,chains1,stag2,chains2,options)
% function [pair_rmsd,msg] = cross_ensemble_pair_rmsd(stag,chains)
%
% Computes a matrix of pairwise rmsd between conformations in an ensemble
%
% stag1     structure tag of the first ensemble
% chains1   vector of single-character chain identifiers of chains to be
%           included in analysis for the first ensemble
% stag2     structure tag of the second ensemble
% chains2   vector of single-character chain identifiers of chains to be
%           included in analysis for the second ensemble
% options   computation options, structure with fields
%           .mode   'backbone' or 'heavy' for backbone or heavy atoms,
%           defaults to 'backbone'
%
% pair_rmsd matrix of pairwise rmsd between conformers (models), symmetric
% msg       error message with fields .error and .text, msg.error = 0
%           indicates normal execution
%
% G. Jeschke, 25.9.2019

global model

maxat = 200000;

pair_rmsd = [];

if ~exist('options','var') || isempty(options)
    options.mode = 'backbone';
end

adr0 = sprintf('[%s](%s)',stag1,chains1(1));
[cindices,msg] = resolve_address(adr0);
if msg.error
    msg.text = sprintf('Address %s cannot be resolved (%s)',adr0,msg.text);
    return
end
n1 = length(model.structures{cindices(1)}(cindices(2)).xyz); % number of conformers
coor1 = cell(1,n1);
for km = 1:n1
    atnum = 0;
    ccoor = zeros(maxat,3);
    for kc = 1:length(chains1)
        atnumc = 0;
        adr = sprintf('[%s](%s)',stag1,chains1(kc));
        [cindices,msg] = resolve_address(adr);
        if msg.error
            msg.text = sprintf('Address %s cannot be resolved (%s)',adr,msg.text);
            return
        end
        if length(cindices) ~= 2
            msg.text = sprintf('Address %s is not a chain address)',adr);
            msg.error = 101;
            return
        end
        if options.mode == 'backbone'
            if model.structures{cindices(1)}(cindices(2)).seqtype == 1
                adr = sprintf('%s{%i}.N,CA,C,O',adr,km);
            elseif model.structures{cindices(1)}(cindices(2)).seqtype == 2
                adr = sprintf('%s{%i}.P,O5'',C5'',C4'',C3'',O3''',adr,km);
            end
        else
            adr = sprintf('%s{%i}',adr,km);
        end
        [indices,msg] = resolve_address(adr);
        if msg.error
            msg.text = sprintf('Address %s cannot be resolved (%s)',adr,msg.text);
            return
        end
        [m,adepth]=size(indices);
        for k = 1:m
            switch adepth
                case 3 % chain model
                    [msg,cc0]=get_chain_model(indices,'xyz_heavy');
                    if msg.error
                        msg.text = 'Chain heavy-atom coordinates could not be retrieved';
                        return
                    end
                    [nat,~] = size(cc0);
                    ccoor(atnum+1:atnum+nat,:) = cc0;
                    atnum = atnum + nat;
                    atnumc = atnumc + nat;
                case 5 % atoms
                    [msg,cc0]=get_atom(indices(k,:),'coor');
                    if msg.error
                        msg.text = 'Atom coordinates could not be retrieved';
                        return
                    end
                    atnum = atnum + 1;
                    atnumc = atnumc + 1;
                    ccoor(atnum,:) = cc0;
                otherwise
                    msg.error = 102;
                    msg.text = 'Address relates neither to chain model nor to atoms';
                    return
            end
        end
        if km == 1
            fprintf(1,'Chain %s in ensemble 1 has %i atoms.\n',chains1(kc),atnumc);
        end
    end
    coor1{km} = ccoor(1:atnum,:);
end


adr0 = sprintf('[%s](%s)',stag2,chains2(1));
[cindices,msg] = resolve_address(adr0);
if msg.error
    msg.text = sprintf('Address %s cannot be resolved (%s)',adr0,msg.text);
    return
end
n2 = length(model.structures{cindices(1)}(cindices(2)).xyz); % number of conformers
coor2 = cell(1,n2);
for km = 1:n2
    atnum = 0;
    ccoor = zeros(maxat,3);
    for kc = 1:length(chains2)
        atnumc = 0;
        adr = sprintf('[%s](%s)',stag2,chains2(kc));
        [cindices,msg] = resolve_address(adr);
        if msg.error
            msg.text = sprintf('Address %s cannot be resolved (%s)',adr,msg.text);
            return
        end
        if length(cindices) ~= 2
            msg.text = sprintf('Address %s is not a chain address)',adr);
            msg.error = 101;
            return
        end
        if options.mode == 'backbone'
            if model.structures{cindices(1)}(cindices(2)).seqtype == 1
                adr = sprintf('%s{%i}.N,CA,C,O',adr,km);
            elseif model.structures{cindices(1)}(cindices(2)).seqtype == 2
                adr = sprintf('%s{%i}.P,O5'',C5'',C4'',C3'',O3''',adr,km);
            end
        else
            adr = sprintf('%s{%i}',adr,km);
        end
        [indices,msg] = resolve_address(adr);
        if msg.error
            msg.text = sprintf('Address %s cannot be resolved (%s)',adr,msg.text);
            return
        end
        [m,adepth]=size(indices);
        for k = 1:m
            switch adepth
                case 3 % chain model
                    [msg,cc0]=get_chain_model(indices,'xyz_heavy');
                    if msg.error
                        msg.text = 'Chain heavy-atom coordinates could not be retrieved';
                        return
                    end
                    [nat,~] = size(cc0);
                    ccoor(atnum+1:atnum+nat,:) = cc0;
                    atnum = atnum + nat;
                    atnumc = atnumc + nat;
                case 5 % atoms
                    [msg,cc0]=get_atom(indices(k,:),'coor');
                    if msg.error
                        msg.text = 'Atom coordinates could not be retrieved';
                        return
                    end
                    atnum = atnum + 1;
                    atnumc = atnumc + 1;
                    ccoor(atnum,:) = cc0;
                otherwise
                    msg.error = 102;
                    msg.text = 'Address relates neither to chain model nor to atoms';
                    return
            end
        end
        if km == 1
            fprintf(1,'Chain %s in ensemble 2 has %i atoms.\n',chains2(kc),atnumc);
        end
    end
    coor2{km} = ccoor(1:atnum,:);
end


pair_rmsd = zeros(n1,n2);

for k1 = 1:n1
    for k2 = 1:n2
        rmsd = rmsd_superimpose(coor1{k1},coor2{k2});
        pair_rmsd(k1,k2) = rmsd;
    end
end