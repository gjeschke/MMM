function [acodes,transmat] = get_RNA_initial_anchor(fragments,anchor,base)
% function [acodes,transmat] = get_RNA_initial_anchor(fragments,anchor,base)
% 
% Determines which fragments fit to the initial anchor geometry with reasonable
% accuracy and provides a transformation matrix from the fragment standrad
% frame to the initial anchor frame
%
% fragments     fragment library, as generated by recompile_library.m
% anchor        coordinates of the C4'(i-1), P(i), and C4'(i) atoms of the
%               initial anchor nucleotide
% base          single-letter code for the initial anchor nucleotide
% 
% acodes        vector of fragment numbers in the library that are allowed
%               for the initial anchor
% transmat      4x4 affine transformation matrix for rotation and
%               translation from the standard frame to the initial
%               anchor frame
%
% G. Jeschke, 24.1.2017

anchor_acc = 0.5; % RMSD accuracy of anchor coordinate fit
O3p_acc = 0.5; % RMSD accuracy of O3' anchor coordinate superposition

[m,~] = size(anchor);
if m >= 6
    aO3p = anchor(6,:);
    anchor = anchor(1:3,:);
end

% if there is no anchor, an empty transformation matrix is returned 
% this allows mk_RNA_loop_backbone to recognize this situation
% all initial fragments are allowed
if isempty(anchor)
    transmat = [];
    acodes = 1:length(fragments);
    return
end

[Ra,origa] =  get_trafo(anchor);
transmat = zeros(4);
transmat(1:3,1:3) = Ra';
transmat(1:3,4) = origa';
transmat(4,4) = 1;

rmsd_vec = zeros(1,length(fragments));
rmsd2_vec = zeros(1,length(fragments));
for k = 1:length(fragments)
    fanchor = fragments(k).(base).coor(fragments(k).(base).assign.previous,:);
    fanchor = [fanchor ones(3,1)];
    fanchor = fanchor*transmat';
    diff = anchor-fanchor(:,1:3);
    rmsd = sqrt(sum(sum(diff.^2))/3);
    rmsd_vec(k) = rmsd;
    if m >= 6
        atags = fragments(k).(base).atomtags;
        for ka = 1:length(atags)
            if strcmp(atags{ka},'O3''')
                pO3p = ka;
            end
        end
        O3p = fragments(k).(base).coor(pO3p,:);
        O3p = [O3p 1]*transmat';
        rmsd2_vec(k) = sqrt(sum((O3p(1:3)-aO3p).^2));
    end
end
codes = 1:length(fragments);
acodes = codes(rmsd_vec < anchor_acc); % consider only fragments that fit initial anchor coordinates with RMSD < anchor_acc
if isnan(min(rmsd_vec)) || isempty(acodes)
    fprintf(1,'Minimal anchor rmsd: %5.2f �\n',min(rmsd_vec));
end
na = length(acodes);
if na < 5
    fprintf(2,'Warning: Only %i fragments fit initial anchor coordinates with requested accuracy of %4.1f �\n',na,anchor_acc);
end
if m >= 6
    rmsd2_vec = rmsd2_vec(rmsd_vec < anchor_acc);
    acodes = acodes(rmsd2_vec < O3p_acc);
    na2 = length(acodes);
    if na2 < 5
        fprintf(2,'Warning: Only %i fragments fit O3'' anchor coordinate with requested accuracy of %4.1f �\n',na,O3p_acc);
    end
end

function [Rp,orig] = get_trafo(coor)

orig = coor(2,:);
coor = coor - repmat(orig,3,1);
x = coor(1,:)-coor(2,:); 
x = x/norm(x);    % unit vector along x
yp = coor(3,:)-coor(2,:); 
yp = yp/norm(yp);
z = cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z = z/norm(z);
y = cross_rowvec(z,x); % real (corrected) y axis
Rp = [x;y;z];

function c=cross_rowvec(a,b)
% A fast cross product that works only for two three-element row vectors 

c = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];
