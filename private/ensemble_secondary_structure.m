function secstruct = ensemble_secondary_structure(file_list,pop)


if ~exist('pop','var') || isempty(pop)
    pop = ones(1,length(file_list));
    pop = pop/sum(pop);
end


for k = 1:length(file_list)
    dssp = get_dssp(file_list{k});
    css = process_dssp(dssp);
    if k == 1
        secstruct.chain = css.chain;
        secstruct.sequence = css.sequence;
        secstruct.resnum = css.resnum;
        secstruct.H = pop(1)*css.H;
        secstruct.B = pop(1)*css.B;
        secstruct.E = pop(1)*css.E;
        secstruct.G = pop(1)*css.G;
        secstruct.I = pop(1)*css.I;
        secstruct.T = pop(1)*css.T;
        secstruct.S = pop(1)*css.S;
        secstruct.none = pop(1)*css.none;
        n = length(css.none);
        secstruct.phi = zeros(n,361);
        secstruct.psi = zeros(n,361);
        secstruct.acc = pop(1)*css.acc;
        for kr = 1:n
            pphi = round(css.phi(kr))+181;
            secstruct.phi(kr,pphi) = pop(1);
            ppsi = round(css.psi(kr))+181;
            secstruct.psi(ppsi) = pop(1);
        end
    else
        secstruct.H = secstruct.H + pop(k)*css.H;
        secstruct.B = secstruct.B + pop(k)*css.B;
        secstruct.E = secstruct.E + pop(k)*css.E;
        secstruct.G = secstruct.G + pop(k)*css.G;
        secstruct.I = secstruct.I + pop(k)*css.I;
        secstruct.T = secstruct.T + pop(k)*css.T;
        secstruct.S = secstruct.S + pop(k)*css.S;
        secstruct.none = secstruct.none + pop(1)*css.none;
        secstruct.acc = secstruct.acc + pop(k)*css.acc;       
        for kr = 1:n
            pphi = round(css.phi(kr))+181;
            secstruct.phi(kr,pphi) = secstruct.phi(kr,pphi) + pop(k);
            ppsi = round(css.psi(kr))+181;
            secstruct.psi(ppsi) = secstruct.psi(ppsi) + pop(k);
        end
    end
end