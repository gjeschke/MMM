% Author : Gerhard Nieuwenhuiys.
% Date of last revision : 2009/01/25.
% The function molecules2.m allows for the calculation of the molecular
% mass of compounds consisting of any number of elements. The element must
% be typed as a recognizable chemical symbol or its corresponding
% atomic number followed by the number of atoms in the molecule of that
% element. Chemical symbols can only be used if the file : atoms.mat is
% activated in the workspace. This is done by simply typing >> load atoms
% in the workspace. All inputs must be separated by commas.
% The function sym2an can also be used for atomic number inputs.
% See also the related functions : sym2an , pertable
% .........................................................................
% Example 1:To calculate the molemass of phenol (C6H6O), type:
% >> molecules2(C,6,H,6,O,1)  ....(after loading atoms.mat !)
% 
% ans =
% 
%   94.1128
%
% Example 2:To calculate the molemass of coppersulphate (CuSO4), type:
% >> molecules2(sym2an('Cu'),1,sym2an('S'),1,8,4)  ....(no need for
%                                                      atoms.mat here)
% ans =
% 
%   159.6036
% .........................................................................
% List of Elements and Atomic Numbers:
% 1(H)        2(He)       3(Li)       4(Be)       5(B)
% 6(C)        7(N)        8(O)        9(F)        10(Ne)
% 11(Na)      12(Mg)      13(Al)      14(Si)      15(P)
% 16(S)       17(Cl)      18(Ar)      19(K)       20(Ca)
% 21(Sc)      22(Ti)      23(V)       24(Cr)      25(Mn)
% 26(Fe)      27(Co)      28(Ni)      29(Cu)      30(Zn)
% 31(Ga)      32(Ge)      33(As)      34(Se)      35(Br)
% 36(Kr)      37(Rb)      38(Sr)      39(Yt)      40(Zr)
% 41(Nb)      42(Mo)      43(Tc)      44(Ru)      45(Rh)
% 46(Pd)      47(Ag)      48(Cd)      49(In)      50(Sn)
% 51(Sb)      52(Te)      53(I)       54(Xe)      55(Cs)
% 56(Ba)      57(La)      58(Ce)      59(Pr)      60(Nd)
% 61(Pm)      62(Sm)      63(Eu)      64(Gd)      65(Tb)
% 66(Dy)      67(Ho)      68(Er)      69(Tm)      70(Yt)
% 71(Lu)      72(Hf)      73(Ta)      74(W)       75(Re)
% 76(Os)      77(Ir)      78(Pt)      79(Au)      80(Hg)
% 81(Tl)      82(Pb)      83(Bi)      84(Po)      85(At)
% 86(Rn)      87(Fr)      88(Ra)      89(Ac)      90(Th)
% 91(Pa)      92(U)       93(Np)      94(Pu)      95(Am)
% 96(Cm)      97(Bk)      98(Cf)      99(Es)      100(Fm)
% .........................................................................
function molmass=molecules2(varargin);
x=cell2mat(varargin);
n=1:0.5*length(x);
x1=x(2*n-1);
x2=x(2*n);
RAMass=[1.0079,4.0026,6.941,9.01218,10.81,12.011,14.0067,15.9994,18.9984,20.179,...
    22.98977,24.305,26.98154,28.086,30.97376,32.06,35.453,39.948,39.098,40.08,...
    44.9559,47.9,50.9414,51.996,54.938,55.847,58.9332,58.7,63.546,65.38,...
    69.72,72.59,74.9216,78.96,79.904,83.8,85.4678,87.62,88.9059,91.22,...
    92.9064,95.94,98.9062,101.07,102.9055,106.4,107.868,112.4,114.82,118.69,...
    121.75,127.6,126.9045,131.3,132.9054,137.34,138.9055,140.12,140.9077,144.24,...
    147,150.4,151.96,157.25,158.9254,162.5,164.9304,167.26,168.9342,173.04,...
    174.97,178.49,180.9479,183.85,186.207,190.2,192.22,195.09,196.9665,200.59,...
    204.37,207.2,208.9804,210,210,222,223,226.0254,227,232.0381,...
    231.0359,238.029,237.0482,244,243,247,247,251,254,257,258,259,262,261,262,...
    266,264,277,268];
molmass=sum(RAMass(x1).*x2);