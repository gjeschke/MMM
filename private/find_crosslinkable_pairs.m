function [KK_pairs,KE_pairs] = find_crosslinkable_pairs(indices,thr)
% [KK_pairs,KE_pairs] = find_crosslinkable_pairs(indices)
%
% Finds crosslinkable residue pairs
% the CA-CA distance is provided for all LYS-LYS and LYS-GLU pairs
%
% indices   array of chain (m,2) or chain model (m,3) indices for analysis,
%           if only chain indices are given, model 1 is assumed for all
%           chains
% thr       upper threshold for the CA-CA distance of pairs to be reported 
%           (Å), defaults to (virtually) infinite distance 
%
% KK_pairs  array (nKK,9) of residue indices (:,1:4), (:5:8) and the CA-CA
%           distance (:,9) (Å) for all LYS-LYS pairs
% KE_pairs  array (nKE,9) of residue indices (:,1:4), (:5:8) and the CA-CA
%           distance (:,9) (Å) for all LYS-GLU pairs
%
% G. Jeschke, 5.11.2018

KK_pairs = [];
KE_pairs = [];

if ~exist('thr','var') || isempty(thr)
    thr = 1e6;
end

[m,n] = size(indices);
if n<2 || n>3
    add_msg_board(sprintf('ERROR in crosslink finder. Wrong index depth (%i)',n));
end

if n == 2
    indices = [indices ones(m,1)];
end

resindices = zeros(10000,4);
rpoi = 0;
for k = 1:m
    [msg,newindices] = get_chain_model(indices(k,:),'children');
    if msg.error
        add_msg_board(sprintf('Warning: No residues in chain model %i',k));
        continue
    end
    [mr,~] = size(newindices);
    resindices(rpoi+1:rpoi+mr,:) = newindices;
    rpoi = rpoi + mr;
end
resindices = resindices(1:rpoi,:);

KKpoi = 0;
KEpoi = 0;
KK_pairs = zeros(rpoi*(rpoi-1)/2,9);
KE_pairs = zeros(rpoi*(rpoi-1)/2,9);

for k1 = 1:rpoi-1
    [~,name1] = get_residue(resindices(k1,:),'name');
    if strcmpi(name1,'LYS') || strcmpi(name1,'GLU')
        for k2 = k1+1:rpoi
            [~,name2] = get_residue(resindices(k2,:),'name');
            if strcmpi(name2,'LYS') || strcmpi(name2,'GLU')
                adr1 = mk_address(resindices(k1,:));
                [~,coor1] = get_object([adr1 '.CA'],'coor');
                adr2 = mk_address(resindices(k2,:));
                [~,coor2] = get_object([adr2 '.CA'],'coor');
                r = norm(coor1-coor2);
                if r <= thr
                    if strcmpi(name1,'LYS') && strcmpi(name2,'LYS')
                        KKpoi = KKpoi + 1;
                        KK_pairs(KKpoi,1:4) = resindices(k1,:);
                        KK_pairs(KKpoi,5:8) = resindices(k2,:);
                        KK_pairs(KKpoi,9) = r;
                    else
                        KEpoi = KEpoi + 1;
                        KE_pairs(KEpoi,1:4) = resindices(k1,:);
                        KE_pairs(KEpoi,5:8) = resindices(k2,:);
                        KE_pairs(KEpoi,9) = r;
                    end
                end
            end
        end
    end
end

KK_pairs = KK_pairs(1:KKpoi,:);
KE_pairs = KE_pairs(1:KEpoi,:);

