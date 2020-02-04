function secstruct = process_dssp(dssp)
% process_dssp(dssp)
%
% Reforms DSSP output obtained by get_dssp.m into a data structure that
% makes further processing easier
%

n = length(dssp);

secstruct.H = zeros(1,n);
secstruct.B = zeros(1,n);
secstruct.E = zeros(1,n);
secstruct.G = zeros(1,n);
secstruct.I = zeros(1,n);
secstruct.T = zeros(1,n);
secstruct.S = zeros(1,n);
secstruct.none = zeros(1,n);
secstruct.resnum = zeros(1,n);
secstruct.phi = zeros(1,n);
secstruct.psi = zeros(1,n);
secstruct.acc = zeros(1,n);
secstruct.sequence = char(32*ones(1,n));
secstruct.chain = char(32*ones(1,n));

for k = 1:n
    secstruct.resnum(k) = str2double(dssp(k).tag);
    secstruct.chain(k) = dssp(k).chain;
    secstruct.sequence(k) = dssp(k).slc;
    secstruct.acc(k) = dssp(k).acc;
    secstruct.phi(k) = dssp(k).phi;
    secstruct.psi(k) = dssp(k).psi;   
    switch dssp(k).sec
        case 'H'
            secstruct.H(k) = 1;
        case 'B'
            secstruct.B(k) = 1;
        case 'E'
            secstruct.E(k) = 1;
        case 'G'
            secstruct.G(k) = 1;
        case 'I'
            secstruct.I(k) = 1;
        case 'T'
            secstruct.T(k) = 1;
        case 'S'
            secstruct.S(k) = 1;
        otherwise
            secstruct.none(k) = 1;
    end
end



