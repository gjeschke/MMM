function anchor = get_anchor_pseudo_torsion(indices,O35p)
% anchor = get_anchor_pseudo_torsion(indices)
% 
% retrieves a coordinate array for the pseudo-torsion defining atom
% coordinates of an anchor nucleotide, C4'(i-1), P(i), C4'(i), P(i+1),
% and, optionally, O5'(i), O3'(i)
% non-existing atoms have NaN coordinates
%
% indices   MMM residue indices for the anchor nucleotide
% O35p      optional flag that requests the O3' and O5' coordinates,
%           defaults to false
% 
% anchor    4x3 or 6x3 coordinate array
%
% G. Jeschke, 1.1.2018

[stag,ctag,modelnum,resnum] = mk_address_parts(indices);

if exist('O35p','var') && ~isempty(O35p) && O35p
    anchor = nan(6,3);
    adr = sprintf('[%s](%s){%i}%i.O5''',stag,ctag,modelnum,resnum);
    [msg,O5p] = get_object(adr,'coor');
    if ~msg.error
        anchor(5,:) = O5p;
    end
    adr = sprintf('[%s](%s){%i}%i.O3''',stag,ctag,modelnum,resnum);
    [msg,O3p] = get_object(adr,'coor');
    if ~msg.error
        anchor(6,:) = O3p;
    end
else
    anchor = nan(4,3);
end


adr = sprintf('[%s](%s){%i}%i.C4''',stag,ctag,modelnum,resnum-1);
[msg,C4pp] = get_object(adr,'coor');
if ~msg.error
    anchor(1,:) = C4pp;
end
adr = sprintf('[%s](%s){%i}%i.P',stag,ctag,modelnum,resnum);
[msg,P] = get_object(adr,'coor');
if ~msg.error
    anchor(2,:) = P;
end
adr = sprintf('[%s](%s){%i}%i.C4''',stag,ctag,modelnum,resnum);
[msg,C4p] = get_object(adr,'coor');
if ~msg.error
    anchor(3,:) = C4p;
end
adr = sprintf('[%s](%s){%i}%i.P',stag,ctag,modelnum,resnum+1);
[msg,Pn] = get_object(adr,'coor');
if ~msg.error
    anchor(4,:) = Pn;
end

