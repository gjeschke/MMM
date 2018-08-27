function links = analyze_RNA_links(snum,shortfrag,cid,fragments,istrand1,estrand1)
% links = analyze_RNA_links(snum,shortfrag,cid,fragments,strand1)
%
% Provides information which fragments of a nucleotide rotamer library
% superimpose at pseudo-torsion angle defining atoms with the first and
% last nucleotide of an RNA piece, 
% corresponding target coordinates for the pseudo-torsion angle defining
% atom C4'(i-1) of the preceding nucleotide and P(j+1) of the following
% nucleotide are generated for each fitting fragment
% see: E. Humphris-Narayanan, A. M. Pyle, J. Mol. Biol. 2012, 421, 6-26
%
% a fit is defined as three-atom rmsd below thresh (currently 0.25
% Angstroem)
%
% snum      number of the structure that contains the RNA piece
% shortfrag pseudo-torsion angle defining atom coordinates for the fragment
%           library
% cid       chain ID of the RNA piece, optional, defaults to 1
% fragments fragments from library
% istrand1  flag, if true, the nucleotide for start attachment is the first
%           nt of the first strand, if false, it is the first nt of the
%           second strand
% estrand1  flag, if true, the nucleotide for end attachment is the last nt
%           of the first strand, if false, it is the last nt of the second
%           strand
%
% links     link definition, struct with fields
%           .ptaf       (3,3) coordinate array of pseudo-torsion defining
%                       atoms of the first nucleotide, P(i), C4'(i), P(i+1)
%           .ptah       (3,3) coordinate array of pseudo-torsion defining
%                       atoms of the last nucleotide of the first strand,
%                       P(i), C4'(i), P(i+1)
%           .ptal       (3,3) coordinate array of pseudo-torsion defining
%                       atoms of the last nucleotide, C4'(j-1), P(j), C4'(j)
%           .previous   vector of fragment library indices of p fitting i-1
%                       fragments
%           .next       vector of fragment library indices of n fitting j+1
%                       fragments
%           .prev_coor  array (p,3) of C4'(i-1) atom coordinates per
%                       fitting i-1 fragment
%           .next_coor  array (n,3) of P(j+1) atom coordinates per fitting
%                       j+1 fragment
%           .ptransmat	affine transformation matrices from standard 
%                       fragment frame to C5'-terminal anchor frame, cell 
%                       vector for p allowed fragments
%           .ntransmat	affine transformation matrix from standard fragment
%                       frame to C3'-terminal anchor frame
%
% G. Jeschke, 24.12.2017

global model

thresh = 0.5; % threshold for coordinate rmsd of pseudotorsion atoms to declare a fragment acceptable

if ~exist('cid','var') || isempty(cid)
    cid = 1;
end

ntnum = length(model.structures{snum}(cid).sequence);
hnum = floor(ntnum/2); 

previous = zeros(1,length(shortfrag));
next = zeros(1,length(shortfrag));
prev_coor = zeros(length(shortfrag),3);
% next_coor = zeros(length(shortfrag),3);
next_P = zeros(length(shortfrag),3);
next_C4p = zeros(length(shortfrag),3);

stag = mk_address_parts(snum);

if istrand1
    adr = sprintf('[%s]1.P',stag);
    [msg,P] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]1.C4''',stag);
    [msg,C4p] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]2.P',stag);
    [msg,Pp] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]2.C4''',stag);
    [msg,C4pp] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
else % initial nucleotide is the first one of strand 2
    adr = sprintf('[%s]%i.P',stag,hnum+1);
    [msg,P] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.C4''',stag,hnum+1);
    [msg,C4p] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.P',stag,hnum+2);
    [msg,Pp] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.C4''',stag,hnum+2);
    [msg,C4pp] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
end
links.ptaf = [P;C4p;Pp];
fptaf = [P;C4p;Pp;C4pp];

minrmsd = 1e6;
for k = 1:length(shortfrag)
    coor = shortfrag{k};
    rmsd = superimpose_3points(links.ptaf,coor(2:4,:));
    if rmsd < minrmsd
        links.ifragment = k;
        minrmsd = rmsd;
    end
end
fprintf(1,'Initial fragment is %i with rmsd %6.3f\n',links.ifragment,minrmsd);


if estrand1
    adr = sprintf('[%s]%i.C4''',stag,hnum-1);
    [msg,C4pp] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.P',stag,hnum);
    [msg,Ppn] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.C4''',stag,hnum);
    [msg,C4pn] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end

    links.ptal = [C4pp;Ppn;C4pn];

    [Ra,origa] =  get_trafo([C4pp;Ppn;C4pn]);
    transmat = zeros(4);
    transmat(1:3,1:3) = Ra';
    transmat(1:3,4) = origa';
    transmat(4,4) = 1;

    links.ntransmat = transmat;
else % (not) strand1 end attachment
    adr = sprintf('[%s]%i.C4''',stag,ntnum-1);
    [msg,C4pp] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.P',stag,ntnum);
    [msg,Ppn] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end
    adr = sprintf('[%s]%i.C4''',stag,ntnum);
    [msg,C4pn] = get_object(adr,'coor');
    if msg.error
        add_msg_board(sprintf('ERROR: Coordinates could not be retrieved for %s (%s). Aborting.',adr,msg.text));
        cd(mydir);
        set(gcf,'Pointer','arrow');
        return
    end

    links.ptal = [C4pp;Ppn;C4pn];

    [Ra,origa] =  get_trafo([C4pp;Ppn;C4pn]);
    transmat = zeros(4);
    transmat(1:3,1:3) = Ra';
    transmat(1:3,4) = origa';
    transmat(4,4) = 1;

    links.ntransmat = transmat;
end;

minrmsd = 1e6;
for k = 1:length(shortfrag)
    coor = shortfrag{k};
    rmsd = superimpose_3points(links.ptal,coor(1:3,:));
    if rmsd < minrmsd
        links.efragment = k;
        minrmsd = rmsd;
    end
end
fprintf(1,'End fragment is %i with rmsd %6.3f\n',links.efragment,minrmsd);

possible = 0;

ptransmat = cell(1,length(shortfrag));
for k = 1:length(shortfrag)
    cfrag = shortfrag{k};
    sfrag = [fragments(k).A.coor(fragments(k).A.assign.previous(1:3),:);...
            fragments(k).A.coor(fragments(k).A.assign.next(2:3),:);];
    [rmsd, ~, transmat] = rmsd_superimpose(fptaf,sfrag(2:5,:));
    cfragb = [sfrag(1:4,:) ones(4,1)]*transmat';
    if rmsd < thresh
%         fprintf(1,'rmsd(%i) %6.2f\n',k,rmsd);
        possible = possible + 1;
        previous(possible) = k;
        prev_coor(possible,:) = cfragb(1,1:3);
        [Ra,origa] =  get_trafo(cfragb(1:3,1:3));
        transmat = zeros(4);
        transmat(1:3,1:3) = Ra';
        transmat(1:3,4) = origa';
        transmat(4,4) = 1;
        ptransmat{possible} = transmat;
    else
%         fprintf(2,'rmsd(%i) %6.2f\n',k,rmsd);
    end
end

links.previous = previous(1:possible);
links.prev_coor = prev_coor(1:possible,:);
links.ptransmat = ptransmat(1:possible);

possible = 0;

for k = 1:length(fragments)
    cfrag = [fragments(k).A.coor(fragments(k).A.assign.previous(1:3),:);...
            fragments(k).A.coor(fragments(k).A.assign.next(2:3),:);];
    [rmsd, ~, transmat] = superimpose_3points(links.ptal,cfrag(1:3,:));
    cfragb = [cfrag ones(5,1)]*transmat';
    if rmsd < thresh
%         fprintf(1,'rmsd(%i) %6.2f\n',k,rmsd);
        possible = possible + 1;
        next(possible) = k;
        next_P(possible,:) = cfragb(4,1:3);
        next_C4p(possible,:) = cfragb(5,1:3);
    else
%         fprintf(2,'rmsd(%i) %6.2f\n',k,rmsd);
    end
end

links.next = next(1:possible);
links.next_P = next_P(1:possible,:);
links.next_C4p = next_C4p(1:possible,:);

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
