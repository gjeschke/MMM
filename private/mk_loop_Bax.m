function [coor,errcode,clashes,nofit] = mk_loop_Bax(sequence, anchorN, anchorNp, prot_coor, loc126)
% backbone based on: H. Sugeta, T. Miyazawa, Biopolymers, 1967, 5, 673-679.
%
% (c) G. Jeschke, 2007-2013

global Ramachandran
global residue_defs
global lovell

clashes = 0;
nofit = 0;

NO_std = [-1.1496,-2.4028,2.8125]; % NO midpoint coordinates in standard frame

clash_threshold = 2.8; % (was 2.5) minimum distance between atoms in non-consecutive residues
constraint_threshold = 6.0; % maximum distance [Å] between projected and experimental label position at residue 126 

if ~exist('maxattempts','var'),
    maxattempts = 1000; % attempts for Ramachandran closure
end;

% Default bond lengths
b1=1.416; % N-C_alpha
b1vec=[b1,0,0];
b2=1.529; % C_alpha-C
b2vec=[b2,0,0];
b3=1.320; % C-N
b3vec=[b3,0,0];

% Default bond angles
fi1=110*pi/180; % N-C_alpha-C
sfi1 = sin(fi1); cfi1 = cos(fi1);
fi2=114.2*pi/180; % C_alpha-C-N
sfi2 = sin(fi2); cfi2 = cos(fi2);
fi3=121*pi/180; % C-N-C_alpha
sfi3 = sin(fi3); cfi3 = cos(fi3);
omega = pi;
so = sin(omega); co = cos(omega);

ngap = length(sequence)-3;
rescodes = zeros(1,ngap+2);

for k = 1:ngap+3,
    rescodes(k) = strfind(residue_defs.single_letter_code,sequence(k));
end;

backbone = zeros(4*(ngap+2)+1,3);
backbone(1:4,:) = anchorN;

phivec = zeros(1,ngap);
psivec = zeros(1,ngap);


% bootstrapping, psi of N terminal anchor not yet defined
C_back = anchorNp(3,:);
N = anchorN(1,:);
CA = anchorN(2,:);
C = anchorN(3,:);
phi = dihedral(C_back,N,CA,C);

% local frame of Calpha of the N-terminal anchor residue
x = anchorN(3,:)-anchorN(2,:);
bondx = norm(x);
x = x/norm(x);
yp = anchorN(2,:)-anchorN(1,:);
yp = yp/norm(yp);
z = cross(x,yp);
z = z/norm(z);
y = cross(z,x);
y = y/norm(y);
A = [x;y;z];
% transformation matrix into that local frame
A = A';
% coordinates of Calpha of the N-terminal anchor residue
acoor = anchorN(2,:);
acoor = acoor';
acoor = acoor + A*[bondx;0;0]; % coordinates of C of the N-terminal anchor residue

failed = true;
attempts = 0;

acoor00 = acoor;
A00 = A;

while failed && attempts < maxattempts,
    attempts = attempts + 1;
    allowed = false;
    while ~allowed,
        indy = rand(1,length(Ramachandran.ephi{rescodes(1)})); % scramble order
        [~,poi] = sort(indy); % since the first fitting phi will be selected
        ephi = Ramachandran.ephi{rescodes(1)}(poi);
        epsi = Ramachandran.epsi{rescodes(1)}(poi);
        [~,poi] = min(abs(ephi-phi)); % find close match of N-terminal anchor phi
        psi = epsi(poi); % set a corresponding psi value
        allowed = Ramachandran_test(phi,psi,sequence(1:2),lovell);
    end;

    ctau=cos(psi);
    stau=sin(psi);
    A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of N of first gap residue
    A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
    A=A*A31; % local frame at N of first gap residue
    
    % loop generation up to Cys-126
    for k = 2:5,
        allowed = false;
        while ~allowed,
            poi = round(rand*Ramachandran.me{rescodes(k)}+0.5);
            rphi=Ramachandran.ephi{rescodes(k)}(poi);
            rpsi=Ramachandran.epsi{rescodes(k)}(poi);
            allowed = Ramachandran_test(rphi,rpsi,sequence(k:k+1),lovell);
        end;
        phivec(k) = 180*rphi/pi;
        psivec(k) = 180*rpsi/pi;
        backbone(4*k-3,:)=acoor'; % coordinates of N
        acoor=acoor+A*b1vec'; % coordinates of C_alpha
        backbone(4*k-2,:)=acoor';
        ctau=cos(rphi);
        stau=sin(rphi);
        A12=[-cfi1,-sfi1,0;sfi1*ctau,-cfi1*ctau,-stau;sfi1*stau,-cfi1*stau,ctau];
        A=A*A12;
        acoor=acoor+A*b2vec'; % coordinates of C
        backbone(4*k-1,:)=acoor';
        ctau=cos(rpsi);
        stau=sin(rpsi);
        A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
        A=A*A23;
        acoor=acoor+A*b3vec'; % coordinates of next N
        A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
        A=A*A31;    
    end;
    backbone(4*k+1,:)=acoor';
    for k = 2:5,
        backbone = add_O(k,backbone);
    end;
    % test for spin label at 126 localization constraint
    % check for label-label constraints

    Ncoor = backbone(4*k-3,:);
    CAcoor = backbone(4*k-2,:);
    Ccoor = backbone(4*k-1,:);

    x=Ncoor-CAcoor; % x axis is along C_alpha-N bond
    x=x/norm(x);    % unit vector along x
    yp=Ccoor-CAcoor; % y axis is in the plane spanned by x axis and C-Ca bond
    yp=yp/norm(yp);
    z=cross(x,yp); % z axis is perpendicular on xy plane
    z=z/norm(z);
    y=cross(z,x); % real (corrected) y axis 
    dircos=[x;y;z];
    Rp=dircos; % rotation matrix for conversion to standard frame
    NO = NO_std*Rp + CAcoor;
    dist = norm(NO - loc126);
    if dist > constraint_threshold,
        failed = true;
        nofit = nofit + 1;
        coor = [];
        acoor = acoor00;
        A = A00;
        continue
    end;
    
    % loop generation after Cys-126
    for k = 6:ngap+2,
        allowed = false;
        while ~allowed,
            poi = round(rand*Ramachandran.me{rescodes(k)}+0.5);
            rphi=Ramachandran.ephi{rescodes(k)}(poi);
            rpsi=Ramachandran.epsi{rescodes(k)}(poi);
            allowed = Ramachandran_test(rphi,rpsi,sequence(k:k+1),lovell);
        end;
        phivec(k) = 180*rphi/pi;
        psivec(k) = 180*rpsi/pi;
        backbone(4*k-3,:)=acoor'; % coordinates of N
        acoor=acoor+A*b1vec'; % coordinates of C_alpha
        backbone(4*k-2,:)=acoor';
        ctau=cos(rphi);
        stau=sin(rphi);
        A12=[-cfi1,-sfi1,0;sfi1*ctau,-cfi1*ctau,-stau;sfi1*stau,-cfi1*stau,ctau];
        A=A*A12;
        acoor=acoor+A*b2vec'; % coordinates of C
        backbone(4*k-1,:)=acoor';
        ctau=cos(rpsi);
        stau=sin(rpsi);
        A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
        A=A*A23;
        acoor=acoor+A*b3vec'; % coordinates of next N
        A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
        A=A*A31;    
    end;
    backbone(4*k+1,:)=acoor';
    for k = 6:ngap+2,
        backbone = add_O(k,backbone);
    end;

    failed = false;
    
    % self-clash test
    k = 1;
    clash = false;
    while ~clash && k <= ngap,
        min_dist = get_min_pair_dist(k,backbone);
        if min_dist < clash_threshold,
            clash = true;
        end;
        k = k + 1;
    end;

    if clash,
        failed = true;
        clashes = clashes + 1;
        errcode = 2;
        coor = [];
        acoor = acoor00;
        A = A00;
    else
        % check for clashes with protein
        pair_dist = get_all_pair_dist(backbone(9:4*ngap+1,:),prot_coor);
        min_dist = min(min(pair_dist));
        if min_dist < clash_threshold,
            failed = true;
            clashes = clashes + 1;
            errcode = 6;
            coor = [];
            acoor = acoor00;
            A = A00;
        end;
    end;
end;
if ~failed,
    errcode = 0;
    coor = backbone;
end;


return

% The following code is for test purposes

% redetermine dihedrals
rphivec = zeros(1,ngap+2);
rpsivec = rphivec;
romvec = rphivec;
for k = 2:ngap+2,
    [phi,psi,omega] = dihedrals(k,backbone);
    rphivec(k) = 180*phi/pi;
    rpsivec(k) = 180*psi/pi;
    romvec(k) = 180*omega/pi;
end;


figure(7); clf;
plot(1:ngap,phivec(2:ngap+1),'k.');
hold on
plot(1:ngap,rphivec(2:ngap+1),'bo');

figure(8); clf;
plot(1:ngap,psivec(2:ngap+1),'k.');
hold on
plot(1:ngap,rpsivec(2:ngap+1),'bo');

figure(9); clf;
hold on
plot(1:ngap,romvec(2:ngap+1),'bo');

disp('Test');


function [phi,psi,omega] = dihedrals(k,backbone)

[m,~] = size(backbone);
CA_back = backbone(4*k-6,:);
C_back = backbone(4*k-5,:);
N = backbone(4*k-3,:);
CA = backbone(4*k-2,:);
C = backbone(4*k-1,:);
omega = dihedral(CA_back,C_back,N,CA);
phi = dihedral(C_back,N,CA,C);
if 4*k+1 > m,
    psi = [];
else
    N_next = backbone(4*k+1,:);
    psi = dihedral(N,CA,C,N_next);
end;

function backbone = add_O(k,backbone)
% carbonyl O position in standard frame (CA-N is x axis, CA-C in xy plane,
% CA is origin)

CA = backbone(4*k-2,:);
C = backbone(4*k-1,:);
Nn = backbone(4*k+1,:);

rO = 1.22;

v1 = CA-C;
v2 = Nn-C;
v1n = v1/norm(v1);
v2n = v2/norm(v2);
z = cross(v1n,v2n);
z = z/norm(z);
Rp = rotmatn(2*pi/3,z);
backbone(4*k,:) = C + rO*v1n*Rp;

function transmat = rotmatn(theta,n)

transmat = zeros(3);
c=cos(theta);
s=sin(theta(1));
t=1-c;
n=n/norm(n);
nx=n(1);
ny=n(2);
nz=n(3);
transmat(1,1)=t*nx^2+c;
transmat(1,2)=t*nx*ny-s*nz;
transmat(1,3)=t*nx*nz+s*ny;
transmat(2,1)=t*nx*ny+s*nz;
transmat(2,2)=t*ny^2+c;
transmat(2,3)=t*ny*nz-s*nx;
transmat(3,1)=t*nx*nz-s*ny;
transmat(3,2)=t*ny*nz+s*nx;
transmat(3,3)=t*nz^2+c;

function min_dist = get_min_pair_dist(k,backbone)

a = backbone(4*k-3:4*k,:);
b = backbone(4*k+5:end,:);

pair_dist = get_all_pair_dist(a,b);
min_dist = min(min(pair_dist));

function pair_dist = get_all_pair_dist(a,b)

[m1,~] = size(a); % get sizes of the coordinates arrays
[m2,~] = size(b);

a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));

function allowed = Ramachandran_test(phi,psi,cn,lovell)
% cn single letter codes for this and next residue

if strcmpi(cn(2),'P'),
    ram_allowed = lovell.prepro;
elseif strcmpi(cn(1),'P')
    ram_allowed = lovell.pro;
elseif strcmpi(cn(1),'G')
    ram_allowed = lovell.gly;
else
    ram_allowed = lovell.general;
end;

phipoi = 181 + round(180*phi/pi);
if phipoi < 1, phipoi = 1; end;
if phipoi >361, phipoi = 361; end;
psipoi = 181 + round(180*psi/pi);
if psipoi < 1, psipoi = 1; end;
if psipoi >361, psipoi = 361; end;

allowed = ram_allowed(psipoi,phipoi);