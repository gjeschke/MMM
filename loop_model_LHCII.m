function [coor,p_model,errcode,acc_res] = loop_model_LHCII(residues, sequence, anchorC, anchorCn, prot_coor, p_thresh, r34, r59, trimer_z, trimer_r_sr, ht_r_sr)
% backbone based on: H. Sugeta, T. Miyazawa, Biopolymers, 1967, 5, 673-679.
%
% residues      residue numbers corresponding to the given sequence
% sequence      sequence in single-letter code, must include C-terminal
%               anchor residue
% anchor_C      coordinates of C-terminal anchor residue
% anchor_Cn     coordinates of next residue after C-terminal anchor residue
%               (required for computation of Ramachandran angle psi of
%               C-terminal anchor residue)
% prot_coor     set of protein coordinates for clash test
% p_thresh      probability threshold below which models are rejected
% r34           mean spin label coordinates of residue 34 (LHCIII specific)
% r59           mean spin label coordinates of residue 59 (LHCIII specific)
% trimer_z      water accessibility constraints
% trimer_r_sr   interprotomer distance constraints
% ht_r_sr       intraprotomer distance constraints from heterogeneous
%               trimers
%
% coor          loop coordinates
% p_model       constraint fulfillment probability
% errcode       error code for modeling failure
%               1 clash with protein
%               2 self-clash
%               3 non-fulfillment of accessibility constraint
%               4 probability below threshold
%               5 constraint definition not consistent with model
% acc_res       residue, at which accessibility constraint was violated
%
% (c) G. Jeschke, 2007-2014

global Ramachandran
global residue_defs

constrained = true; % allows for switching off the constraints 
acc_constraints = true; % consider accessibility constraints
tri_constraints = true; % consider trimer constraints for singly labeled protomers
het_tri_constraints = true; % consider constraints from heterogeneous trimers

errcode = 0;
p_model = 1;
acc_res = 0;

% Cterm2 = sequence(end); % needed to recognize pre-proline Ramachandran case 
sequence = sequence(1:end-1);

% if strcmpi(Cterm2,'P'),
%     lovell = 'Lovell_Pre_Proline';
% elseif strcmpi(sequence(end),'P')
%     lovell = 'Lovell_Proline';
% elseif strcmpi(sequence(end),'G')
%     lovell = 'Lovell_Glycine';
% else
%     lovell = 'Lovell_general';
% end;
% 
% load(lovell);

accstd = 0.5; % standard deviation for Monte Carlo Metropolis
accrad = 5.0; % 3.0 maximum accepted distance from target Calpha

clash_threshold = 2.5; % 2.8 (was 2.5) minimum distance between atoms in non-consecutive residues
clash_threshold_lp = 2.5; % 2.8 (was 2.5) minimum distance between loop and protein residues (except terminal residues)

NO_std = [-1.1496,-2.4028,2.8125]; % NO midpoint coordinates in standard frame

dextended = 3.8; % extended mean Calpha-Calpha distance
dhelix = 1.5; % helical (contracted) mean Calpha-Calpha distance
dmean = (dextended + dhelix)/2; % mean Calpha-Calpha distance step

resattempts = 100; % attempts for single-residue step

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
omega = -pi; % sign does not matter, but taken negative for consistency
so = sin(omega); co = cos(omega);

ngap = length(sequence)-1;
rescodes = zeros(1,length(sequence));

for k = 1:length(sequence),
    rescodes(k) = strfind(residue_defs.single_letter_code,sequence(k));
end;

backbone = zeros(4*(ngap+1),3); % ### corrected
backbone(end-3:end,:) = anchorC; % ### corrected

phivec = zeros(1,ngap);
psivec = zeros(1,ngap);


rejections = zeros(1,ngap+2);
drvec = rejections;
dRvec = rejections;

pop = Ramachandran.pop{1};
[mr,~] = size(pop);

% phi = pi*Ramachandran.phi/180;
% psi = pi*Ramachandran.psi/180;


% bootstrapping, phi of C terminal anchor not yet defined
N_forth = anchorCn(1,:);
N = anchorC(1,:);
CA = anchorC(2,:);
C = anchorC(3,:);
psi = dihedral_fast(N,CA,C,N_forth);

% indy = rand(1,length(Ramachandran.epsi{rescodes(end)})); % scramble order
% [~,poi] = sort(indy); % since the first fitting phi will be selected
% ephi = Ramachandran.ephi{rescodes(end)}(poi);
% epsi = Ramachandran.epsi{rescodes(end)}(poi);
poi = round(rand*Ramachandran.me{rescodes(end)}+0.5);
ephi=Ramachandran.ephi{rescodes(end)}(poi);
epsi=Ramachandran.epsi{rescodes(end)}(poi);

[~,poi] = min(abs(epsi-psi)); % find close match of C-terminal anchor psi
phi = ephi(poi); % set a corresponding psi value

% local frame of Calpha of the N-terminal anchor residue
x = anchorC(1,:)-anchorC(2,:);
bondx = norm(x);
x = x/norm(x);
yp = anchorC(2,:)-anchorC(3,:);
yp = yp/norm(yp);
z = cross_rowvec(x,yp);
z = z/norm(z);
y = cross_rowvec(z,x);
y = y/norm(y);
A = [x;y;z];
% transformation matrix into that local frame
A = A';
% coordinates of Calpha of the N-terminal anchor residue
acoor = anchorC(2,:);
acoor = acoor';
acoor = acoor + A*[bondx;0;0]; % coordinates of N of the C-terminal anchor residue

ctau=cos(-phi); % sign changed for 'backward' extension
stau=sin(-phi); % sign changed for 'backward' extension
A23=[-cfi3,-sfi3,0;sfi3*ctau,-cfi3*ctau,-stau;sfi3*stau,-cfi3*stau,ctau];
A=A*A23;
acoor=acoor+A*b3vec'; % coordinates of C of first missing residue
A31=[-cfi2,-sfi2,0;sfi2*co,-cfi2*co,-so;sfi2*so,-cfi2*so,co];
A=A*A31; % local frame at C of first missing residue

pcum = 1;

% loop generation
for k = length(sequence)-1:-1:1,
    resnum = residues(k);
    poi = round(rand*Ramachandran.me{rescodes(k)}+0.5);
    rphi=Ramachandran.ephi{rescodes(k)}(poi);
    rpsi=Ramachandran.epsi{rescodes(k)}(poi);
    phivec(k) = 180*rphi/pi;
    psivec(k) = 180*rpsi/pi;
    backbone(4*k-1,:)=acoor'; % coordinates of C
    acoor=acoor+A*b2vec'; % coordinates of C_alpha
    [mi,poi] = min(abs(trimer_z(:,1) - resnum));
    if mi == 0,
        p = accessibility_constraint(acoor,trimer_z(poi,2),trimer_z(poi,3));
        if p < 0.5,
            coor = [];
            p_model = 0;
            errcode = 3;
            acc_res = resnum;
            return
        end;
    end;
    backbone(4*k-2,:)=acoor';
    ctau=cos(-rpsi); % sign changed for 'backward' prediction
    stau=sin(-rpsi); % sign changed for 'backward' prediction
    A12=[-cfi1,-sfi1,0;sfi1*ctau,-cfi1*ctau,-stau;sfi1*stau,-cfi1*stau,ctau];
    A=A*A12;
    acoor=acoor+A*b1vec'; % coordinates of N
    backbone(4*k-3,:)=acoor';
    [mi,poi] = min(abs(trimer_r_sr(:,1) - resnum));
    if mi == 0,
        Ncoor = backbone(4*k-3,:);
        CAcoor = backbone(4*k-2,:);
        Ccoor = backbone(4*k-1,:);
        x=Ncoor-CAcoor; % x axis is along C_alpha-N bond
        x=x/norm(x);    % unit vector along x
        yp=Ccoor-CAcoor; % y axis is in the plane spanned by x axis and C-Ca bond
        yp=yp/norm(yp);
        z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
        z=z/norm(z);
        y=cross_rowvec(z,x); % real (corrected) y axis 
        dircos=[x;y;z];
        Rp=dircos; % rotation matrix for conversion to standard frame
        NO = NO_std*Rp + CAcoor;
        p = trimer_probability(NO,trimer_r_sr(poi,2),trimer_r_sr(poi,3));
        pcum = pcum * p;
        if pcum < 0.5,
            coor = [];
            p_model = 0;
            errcode = 4;
            acc_res = resnum;
            return
        end;
    end;
    ctau=cos(-rphi); % sign changed for 'backward' prediction
    stau=sin(-rphi); % sign changed for 'backward' prediction
    A23=[-cfi3,-sfi3,0;sfi3*ctau,-cfi3*ctau,-stau;sfi3*stau,-cfi3*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of next N
    A31=[-cfi2,-sfi2,0;sfi2*co,-cfi2*co,-so;sfi2*so,-cfi2*so,co];
    A=A*A31;    
end;


if constrained,
    
    if acc_constraints,
        % check for accessibility constraints

        [m,~] = size(trimer_z);
        for k = 1:m,
            res = trimer_z(k,1);
            [~,poi] = min(abs(residues - res));
            coor = backbone(4*poi-2,:);
            p = accessibility_constraint(coor,trimer_z(k,2),trimer_z(k,3));
            if p < 0.5, % actually it is zero or 1
                errcode = 3;
                p_model = 0;
                coor = [];
                acc_res = res;
                return
            end;
        end;
    end;
    
    if tri_constraints,
        % check for intertrimer constraints

        [m,~] = size(trimer_r_sr);
        for kc = 1:m,
            res = trimer_r_sr(kc,1);
            [~,k] = min(abs(residues - res));
            Ncoor = backbone(4*k-3,:);
            CAcoor = backbone(4*k-2,:);
            Ccoor = backbone(4*k-1,:);

            x=Ncoor-CAcoor; % x axis is along C_alpha-N bond
            x=x/norm(x);    % unit vector along x
            yp=Ccoor-CAcoor; % y axis is in the plane spanned by x axis and C-Ca bond
            yp=yp/norm(yp);
            z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
            z=z/norm(z);
            y=cross_rowvec(z,x); % real (corrected) y axis 
            dircos=[x;y;z];
            Rp=dircos; % rotation matrix for conversion to standard frame
            NO = NO_std*Rp + CAcoor;
            p = trimer_probability(NO,trimer_r_sr(kc,2),trimer_r_sr(kc,3));
            p_model = p_model * p;
        end;

        if p_model < p_thresh,
            errcode = 4;
            coor = [];
            return
        end;
    end;
    
    if het_tri_constraints,
        % check for intratrimer constraints

        [m,~] = size(ht_r_sr);
        for kc = 1:m,
            res = ht_r_sr(kc,1);
            [~,k] = min(abs(residues - res));
            Ncoor = backbone(4*k-3,:);
            CAcoor = backbone(4*k-2,:);
            Ccoor = backbone(4*k-1,:);

            x=Ncoor-CAcoor; % x axis is along C_alpha-N bond
            x=x/norm(x);    % unit vector along x
            yp=Ccoor-CAcoor; % y axis is in the plane spanned by x axis and C-Ca bond
            yp=yp/norm(yp);
            z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
            z=z/norm(z);
            y=cross_rowvec(z,x); % real (corrected) y axis 
            dircos=[x;y;z];
            Rp=dircos; % rotation matrix for conversion to standard frame
            NO = NO_std*Rp + CAcoor;
            switch ht_r_sr(kc,2)
                case 34
                    ref_res = r34;
                case 59
                    ref_res = r59;
                otherwise
                    errcode = 5;
                    coor = [];
                    p_model = 0;
                    return
            end
            p = het_trimer_probability(NO,ref_res,ht_r_sr(kc,3),ht_r_sr(kc,4));
            p_model = p_model * p;
        end;

        if p_model < p_thresh,
            errcode = 4;
            coor = [];
            return
        end;
    end;
end;

% if the model survived up to this point, it has sufficient probability of
% fulfilling all constraints (or flag 'constrained' was set to false)

% generate backbone oxygen coordinates
for k = 1:length(sequence)-1,
    backbone = add_O(k,backbone);
end;

% check for clashes with protein
pair_dist = get_all_pair_dist(backbone(1:4*(length(sequence)-2),:),prot_coor);
min_dist = min(min(pair_dist));
if min_dist < clash_threshold_lp,
    errcode = 1;
    p_model = 0;
    coor = [];
    return
end;

% self-clash test
k = 1;
clash = false;
while ~clash && k <= length(sequence)-2,
    min_dist = get_min_pair_dist(k,backbone);
    if min_dist < clash_threshold,
        clash = true;
    end;
    k = k + 1;
end;

if clash,
    errcode = 2;
    coor = [];
    p_model = 0;
    return
end;

coor = backbone;

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

[m,~] = size(backbone);
CA_back = backbone(4*k-6,:);
C_back = backbone(4*k-5,:);
N = backbone(4*k-3,:);
CA = backbone(4*k-2,:);
C = backbone(4*k-1,:);
omega = dihedral_fast(CA_back,C_back,N,CA);
phi = dihedral_fast(C_back,N,CA,C);
if 4*k+1 > m,
    psi = [];
else
    N_next = backbone(4*k+1,:);
    psi = dihedral_fast(N,CA,C,N_next);
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
z = cross_rowvec(v1n,v2n);
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

% the following constraint test functions assume coordinates in Å units and
% constraints in nm units

function p = accessibility_constraint(coor,lb,ub)
% Given the estimated Calpha(!) coordinate, it is checked whether the
% model fulfils (p = 1) or does not fulfil (p = 0) the water accessibility
% constraint

z = abs(coor(3)/10);
p = 1;
if z < lb || z > ub,
    p = 0;
end;
               
function p = trimer_probability(coor,mr,sr)
% Given the estimated spin label coordinate in a single protomer, the
% distance to the neighboring protomer is computed assuming C3 symmetry
% with z being the symmetry axis,
% the probability of finding this distance is computed with mr and sr being
% the mean value and width of a Gaussian probability density distribution 
        
xl = coor(1)/10;
yl = coor(2)/10;
r = sqrt(3*(xl^2+yl^2)); % see p. 6 Lab book 2014/1
p = exp(-((r-mr)/sr)^2);

function p = het_trimer_probability(coor1,coor2,mr,sr)
% Given the estimated spin label coordinates, the  label-to-label
% distance is computed 
% the probability of finding this distance is computed with mr and sr being
% the mean value and width of a Gaussian probability density distribution 
        
r = norm(coor1-coor2)/10; 
p = exp(-((r-mr)/sr)^2);

