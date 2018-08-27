function model_LHCII_N_terminus
% Computation is based on PDB structure 2BHW of pea LHCII
% the structure was transformed into a frame where the z axis is the C3
% symmetry axis and membrane normal
% a bilayer model was fitted and the structure further transformed so that
% z = 0 is the membrane center
% residues 10-13 were removed as they feature poor electron density and
% the intratrimer distance constraint for residue 12 matches poorly, while
% those ones for residues 10 and 11 do not match at all
%
% starting from anchoring residue 14, loop models are generated to residue
% 3 (N-terminal extension), on the way probabilities pj are computed for 
% constraint j to be fulfilled, total model probability p_model is the 
% product of all pj
%
% models are rejected if p_model is below a threshold p_thresh
%
% models with sufficient probability are tabulated along with their
% probability

global general
global Ramachandran

addpath(genpath(pwd));

rng('shuffle'); % initialize random number generator to be able to obtain different ensembles in subsequent runs

load([general.Ramachandran 'Ramachandran']);

sequence = 'SATTKKVASSGS'; % last residue must be the achor residue, here Ser-14
residues = 3:14; % residue numbers corresponding to the sequence

p_thresh = 0.01; % probability threshold for accepting a model

% predicted mean label positions for residues 34 and 59 (see lab bokk
% 2014/1, p. 6
r34 = [13.0950   -0.2600  -31.9150];
r59 = [-0.5050   15.2450  -26.9300];

% water accessibility constraints for singly labeled trimers, see
% laboratory Notebook 2014/1, p. 3-5
% format: residue, minimum z coordinate (nm), max. z coordinate (nm)
% z axis is the membrane normal, origin is in the membrane center
% full water accessibility is coded by a 10 nm upper limit
%
% loop model is rejected if lower or upper bound is violated
trimer_z = [3,2.23,10;...
            4,2.13,2.75;...
            7,2.18,2.95;...
            9,2.14,2.76;...
            10,2.23,10.0;...
            11,2.08,2.69;...
            12,2.25,10.0];
        
% distance constraints (Gaussian distribution) for singly labeled trimers
% format: residue, mean distance <r> (nm), Gaussian width sr (nm) 
% p = e^[(r - <r>)/sr]^2
%
% loop model is assigned a probability corresponding to the value of the
% Gaussian distribution at the distance in the model
trimer_r_sr = [3,3.65,1.28;...
          4,3.72,1.21;...
          7,4.00,1.18;...
          9,4.31,1.34;...
          10,4.64,1.46;...
          11,4.57,1.17;...
          12,4.96,1.36];
      
% distance constraints (Gaussian) distribution for doubly-labeled protomers
% in heterogenoeus trimers
% format: residue 1, residue 2, mean distance <r> (nm), width sr (nm)
% p = e^[(r - <r>)/sr]^2
%
% loop model is assigned a probability corresponding to the value of the
% Gaussian distribution at the distance in the model
het_trimer_r_sr = [3,34,2.66,1.58;...
                   3,59,1.50,1.73;...
                   7,34,2.54,1.41;...
                   7,59,1.64,1.86;...
                   11,34,2.50,1.22;...
                   11,59,1.70,2.13];

               ht_r_sr = [3,34,2.66,1.60;...
                   3,59,1.50,1.66;...
                   7,34,2.67,1.62;...
                   7,59,2.75,1.90;...
                   11,34,2.45,1.56;...
                   11,59,1.50,1.73];

[coor,p_model,errcode] = loop_model_LHCII(residues, sequence, anchorC, anchorCn, prot_coor, p_thresh, r34, r59, trimer_z, trimer_r_sr, ht_r_sr);
               
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