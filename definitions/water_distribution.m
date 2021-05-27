function water_distribution
% computes the distribution functions of molecular groups of DOPC and DPPC
% bilayers according to N. Kucerka et al., Biophys. J. 2008, 95, 2356-2367,
% in particular, the water distributions are needed for ESEEM D2O
% accessibility predictions
%
% G. Jeschke, 2010

% values from Table 3 of Kucerka et al., third column
DPPC.VL=1229;
DPPC.VHL=331;
DPPC.RCG=0.41;
DPPC.RPCN=0.28;
DPPC.r=1.93;
DPPC.r12=0;
DPPC.DB=39.0;
DPPC.DHH=38.0;
DPPC.DC2=28.4;
DPPC.DH1=4.7;
DPPC.A=63.1;
DPPC.zCG=14.7;
DPPC.sigCG=2.11;
DPPC.zPCN=19.7;
DPPC.sigPCN=2.62;
DPPC.zCholCH3=21.6;
DPPC.sigCholCH3=2.98;
DPPC.zCH=0;
DPPC.sigCH=1;
DPPC.sigHC=2.47; % they forgot to square the value, as is apparent from Figure 5, we fix that in line 156
DPPC.zCH3=0;
DPPC.sigCH3=2.73;

% number of CH, CH2, and CH3 groups in lipid tails
DPPC.nCH2=28;
DPPC.nCH=0;
DPPC.nCH3=2;

% values from Table 3 of Kucerka et al., fourth column
DOPC.VL=1303;
DOPC.VHL=331;
DOPC.RCG=0.42;
DOPC.RPCN=0.26;
DOPC.r=1.96;
DOPC.r12=0.79;
DOPC.DB=38.7;
DOPC.DHH=36.7;
DOPC.DC2=28.8;
DOPC.DH1=3.9;
DOPC.A=67.4;
DOPC.zCG=14.8;
DOPC.sigCG=2.05;
DOPC.zPCN=19.1;
DOPC.sigPCN=2.41;
DOPC.zCholCH3=20.6;
DOPC.sigCholCH3=2.98;
DOPC.zCH=9.60;
DOPC.sigCH=3.05;
DOPC.sigHC=2.48;
DOPC.zCH3=0;
DOPC.sigCH3=3.09;

% number of CH, CH2, and CH3 groups in lipid tails
DOPC.nCH2=28;
DOPC.nCH=4;
DOPC.nCH3=2;

z=linspace(0,35,351);
lip=DPPC;

figure(1); clf;
hold on;


VCG=lip.RCG*lip.VHL;
cCG=VCG/(lip.A*lip.sigCG);
PCG=Pgauss(z,cCG,lip.zCG,lip.sigCG);

plot(z,PCG,'m');

VPCN=lip.RPCN*lip.VHL;
cPCN=VPCN/(lip.A*lip.sigPCN);
PPCN=Pgauss(z,cPCN,lip.zPCN,lip.sigPCN);

plot(z,PPCN,'r');

VCholCH3=lip.VHL-VCG-VPCN;
cCholCH3=VCholCH3/(lip.A*lip.sigCholCH3);

PCholCH3=Pgauss(z,cCholCH3,lip.zCholCH3,lip.sigCholCH3);

plot(z,PCholCH3,'r:');

PHC=Perf(z,lip.DC2/2,lip.sigHC^2);

plot(z,PHC,'k');

VCH2=(lip.VL-lip.VHL)/(lip.nCH2+lip.r12*lip.nCH+lip.r*lip.nCH3);
VCH3=lip.r*VCH2;
VCH=lip.r12*VCH2;

cCH=lip.nCH*VCH/(lip.A*lip.sigCH);
PCH=Pgauss(z,cCH,lip.zCH,lip.sigCH);
plot(z,PCH,'y');

cCH3=lip.nCH3*VCH3/(lip.A*lip.sigCH3);
PCH3=Pgauss(z,cCH3,lip.zCH3,lip.sigCH3);
plot(z,PCH3,'g');

PCH2=PHC-PCH3-PCH;
plot(z,PCH2,'b');

PW=1-PCG-PPCN-PCholCH3-PHC;
PW(PW<0)=0;

plot(z,PW,'c');

axis([0,35,0,1.05]);

DPPC.z=z;
DPPC.PW=PW;
DPPC.PCH2=PCH2;
DPPC.PCH3=PCH3;
DPPC.PCH=PCH;
DPPC.PHC=PHC;
DPPC.PCholCH3=PCholCH3;
DPPC.PCG=PCG;
DPPC.PPCN=PPCN;

% plot(z,PHC+PCG,'b:');

lip=DOPC;

figure(2); clf;
hold on;


VCG=lip.RCG*lip.VHL;
cCG=VCG/(lip.A*lip.sigCG);
PCG=Pgauss(z,cCG,lip.zCG,lip.sigCG);

plot(z,PCG,'m');

VPCN=lip.RPCN*lip.VHL;
cPCN=VPCN/(lip.A*lip.sigPCN);
PPCN=Pgauss(z,cPCN,lip.zPCN,lip.sigPCN);

plot(z,PPCN,'r');

VCholCH3=lip.VHL-VCG-VPCN;
cCholCH3=VCholCH3/(lip.A*lip.sigCholCH3);

PCholCH3=Pgauss(z,cCholCH3,lip.zCholCH3,lip.sigCholCH3);

plot(z,PCholCH3,'r:');

PHC=Perf(z,lip.DC2/2,lip.sigHC^2);

plot(z,PHC,'k');

VCH2=(lip.VL-lip.VHL)/(lip.nCH2+lip.r12*lip.nCH+lip.r*lip.nCH3);
VCH3=lip.r*VCH2;
VCH=lip.r12*VCH2;

cCH=lip.nCH*VCH/(lip.A*lip.sigCH);
PCH=Pgauss(z,cCH,lip.zCH,lip.sigCH);
plot(z,PCH,'y');

cCH3=lip.nCH3*VCH3/(lip.A*lip.sigCH3);
PCH3=Pgauss(z,cCH3,lip.zCH3,lip.sigCH3);
plot(z,PCH3,'g');

PCH2=PHC-PCH3-PCH;
plot(z,PCH2,'b');

PW=1-PCG-PPCN-PCholCH3-PHC;
PW(PW<0)=0;
plot(z,PW,'c');


axis([0,35,0,1.05]);

DOPC.z=z;
DOPC.PW=PW;
DOPC.PCH2=PCH2;
DOPC.PCH3=PCH3;
DOPC.PCH=PCH;
DOPC.PHC=PHC;
DOPC.PCholCH3=PCholCH3;
DOPC.PCG=PCG;
DOPC.PPCN=PPCN;

save bilayer_definitions DPPC DOPC

function Pi=Pgauss(z,ci,zi,sigi)
% double Gaussian distribution, Eq. (2) of Kucerka et al.
Pi=ci*(exp(-(z+zi).^2/(2*sigi^2))+exp(-(z-zi).^2/(2*sigi^2)))/sqrt(2*pi);

function PHC=Perf(z,zi,sigi)

arg1=(z+zi)/sqrt(2*sigi);
arg2=(z-zi)/sqrt(2*sigi);
PHC=(erf(arg1)-erf(arg2))/2;
