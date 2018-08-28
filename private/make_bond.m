function chain = make_bond(chain,atnum1,atnum2)
% private function that makes a bond in rd_pdb
% chain  structure describing a single chain
% atnum1 atom for which the bond is stored
% atnum2 binding partner

bonds = chain.conn(atnum1,:);

% If this bond was not already stored, add bond
if ~any(bonds==atnum2)
    
    bondidx = sum(bonds~=0) + 1;
    
    % update maximum size of connection table
    if bondidx > chain.maxconn
        chain.maxconn = bondidx;
    end
    
    % add bond
    chain.conn(atnum1,bondidx) = atnum2;
    
end
