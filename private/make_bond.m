function chain=make_bond(chain,atnum1,atnum2)
% private function that makes a bond in rd_pdb
% chain  structure describing a single chain
% atnum1 atom for which the bond is stored
% atnum2 binding partner

bonds=chain.conn(atnum1,:);
n = length(bonds);
maxbonds=sum(bonds~=0);
if isempty(find(bonds==atnum2,1)) % check if this bond was already stored
    maxbonds=maxbonds+1;
    if maxbonds > chain.maxconn, % update maximum size of connection table
        chain.maxconn=maxbonds;
    end;
    bonds(maxbonds)=atnum2;
end;
if length(bonds) > n,
    [m,n] = size(chain.conn);
    conn = zeros(m,length(bonds));
    conn(1:m,1:n) = chain.conn;
    conn(atnum1,:) = bonds;
    chain.conn = conn;
else
    chain.conn(atnum1,:)=bonds;
end;
