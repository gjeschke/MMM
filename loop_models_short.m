function [coor,errcode] = loop_models_short(sequence, anchorN, anchorC, anchorNp, anchorCn, maxattempts)
% backbone based on: H. Sugeta, T. Miyazawa, Biopolymers, 1967, 5, 673-679.
%
% (c) G. Jeschke, 2007-2013

global Ramachandran
global residue_defs

silent = true;

errcode = 0;

Cterm2 = sequence(end); % needed to recognize pre-proline Ramachandran case 
sequence = sequence(1:end-1);

if strcmpi(Cterm2,'P'),
    lovell = 'Lovell_Pre_Proline';
elseif strcmpi(sequence(end),'P')
    lovell = 'Lovell_Proline';
elseif strcmpi(sequence,'G')
    lovell = 'Lovell_Glycine';
else
    lovell = 'Lovell_general';
end;

load(lovell);

accstd = 0.5; % standard deviation for Monte Carlo Metropolis
accrad = 1; % maximum accepted distance from target Calpha

clash_threshold = 2.8; % (was 2.5) minimum distance between atoms in non-consecutive residues


if ~exist('maxattempts','var'),
    maxattempts = 1000; % attempts to generate second half loop
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

ngap = length(sequence)-2;
rescodes = zeros(1,ngap+2);

for k = 1:ngap+2,
    rescodes(k) = strfind(residue_defs.single_letter_code,sequence(k));
end;

backbone = zeros(4*(ngap+3),3);
backbone(1:4,:) = anchorN;
backbone(end-7:end-4,:) = anchorC;
backbone(end-3:end,:) = anchorCn;

phivec = zeros(1,ngap);
psivec = zeros(1,ngap);

nout = floor(ngap/2);

rejections = zeros(1,ngap+2);
drvec = rejections;
dRvec = rejections;

pop = Ramachandran.pop{1};
[mr,~] = size(pop);

% phi = pi*Ramachandran.phi/180;
% psi = pi*Ramachandran.psi/180;


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

accepted = false;

failed = 0;

while ~accepted && failed < maxattempts,
    indy = rand(1,length(Ramachandran.ephi{rescodes(1)})); % scramble order
    [~,poi] = sort(indy); % since the first fitting phi will be selected
    ephi = Ramachandran.ephi{rescodes(1)}(poi);
    epsi = Ramachandran.epsi{rescodes(1)}(poi);
    [~,poi] = min(abs(ephi-phi)); % find close match of N-terminal anchor phi
    psi = epsi(poi); % set a corresponding psi value

    ctau=cos(psi);
    stau=sin(psi);
    A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of N of first gap residue
    A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
    A=A*A31; % local frame at N of first gap residue

    acoor00 = acoor;
    A00 = A;
    % loop generation
    for k = 2:ngap+1,
        poi = round(rand*Ramachandran.me{rescodes(k)}+0.5);
        rphi=Ramachandran.ephi{rescodes(k)}(poi);
        rpsi=Ramachandran.epsi{rescodes(k)}(poi);
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
    if norm(acoor'-anchorC(1,:)) > accrad,
        A = A00;
        acoor = acoor00;
        failed = failed +1;
    end;
end;

if failed >= maxattempts,
    if ~silent,
        fprintf(2,'Maximum attempts without finding closed loop solution\n');
    end;
    errcode = 1;
    coor = [];
    return
else
    if ~silent,
        fprintf(1,'%i attempts required to generate loop model.\n',failed + 1);
    end;
end;

% redetermine dihedrals
rphivec = zeros(1,ngap+2);
rpsivec = rphivec;
romvec = rphivec;
for k = 2:ngap+1,
    [phi,psi,omega] = dihedrals(k,backbone);
    rphivec(k) = 180*phi/pi;
    rpsivec(k) = 180*psi/pi;
    romvec(k) = 180*omega/pi;
end;
backbone(4*(ngap+2)-3,:) = acoor';
acoor=acoor+A*b1vec'; % coordinates of C_alpha
backbone(4*(ngap+2)-2,:) = acoor';

% stretch backbone to make C-terminal anchor CA overlap
R0 = anchorC(2,:) - acoor'; % difference vector
steps = 3*ngap + 1; 
dR0 = R0/steps;
for k = 1:steps,
    res = 2 + floor(k/3); % number of residue to be corrected, starts with first gap residue 2
    poi = 4*(res-1) + 1 + mod(k,3);
    backbone(poi,:) = backbone(poi,:) + k*dR0;
end;

for k = 2:ngap+1,
    backbone = add_O(k,backbone);
end;

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
    if ~silent,
        fprintf(2,'Self clash of loop at %4.2f A.\n',min_dist);
    end;
    errcode = 2;
    coor = [];
    return
end;
coor = backbone;

[phi,psi,omega] = dihedrals(ngap+2,backbone);

phipoi = 181 + round(180*phi/pi);
if phipoi < 1, phipoi = 1; end;
if phipoi >361, phipoi = 361; end;
psipoi = 181 + round(180*psi/pi);
if psipoi < 1, psipoi = 1; end;
if psipoi >361, psipoi = 361; end;

allowed = ram_allowed(psipoi,phipoi);

if ~allowed,
    errcode = 3;
    if ~silent,
        fprintf(2,'C-terminal anchor outside allowed Ramachandran region.\n');
    end;
    return
end;

return

% The following code is for test purposes

% redetermine dihedrals
cphivec = zeros(1,ngap+2);
cpsivec = cphivec;
comvec = cphivec;
for k = 2:ngap+1,
    [phi,psi,omega] = dihedrals(k,backbone);
    cphivec(k) = 180*phi/pi;
    cpsivec(k) = 180*psi/pi;
    comvec(k) = 180*omega/pi;
end;

figure(7); clf;
plot(1:ngap,phivec(2:ngap+1),'k.');
hold on
plot(1:ngap,rphivec(2:ngap+1),'bo');
plot(1:ngap,cphivec(2:ngap+1),'rx');

figure(8); clf;
plot(1:ngap,psivec(2:ngap+1),'k.');
hold on
plot(1:ngap,rpsivec(2:ngap+1),'bo');
plot(1:ngap,cpsivec(2:ngap+1),'rx');

figure(9); clf;
hold on
plot(1:ngap,romvec(2:ngap+1),'bo');
plot(1:ngap,comvec(2:ngap+1),'rx');

function [phi,psi,omega] = dihedrals(k,backbone)

[m,n] = size(backbone);
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
