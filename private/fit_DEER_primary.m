function [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax,distr,dfile,options)
% [texp,vexp,deer,bckg,param] = fit_DEER_primary(rax,distr,deer_file,options)
%
% Finds best background and modulation depth fit for given primary DEER
% data and given distance distribution
%
% rax       distance axis of the distance distribution
% distr     distance distribution
% dfile     file name of the DEER data file (EleXsys format, no extension)
% options   fit options, structure with fields
%           .depth  if empty, modulation depth is freely fitted, default
%                   if single value, this modulation depth is used
%                   if array(1,2), fit range for modulation depth
%           .dim    if empty, background dimension is 3, default
%                   if single value, this dimension is used
%                   if array(1,2), fit range for background dimension
%           .kdec   if empty, decay constant is freely fitted, default
%                   if single value, this decay constant is used
%                   if array(1,2), fit range for background dimension
%           .cutoff fraction of the time axis, where the data is cut off
%                   between 0.5 (enforced) and 1, defaults to 1
%
% texp      time axis of zero-time corrected DEER data
% vexp      real part of zero-time and phase-corrected experimental data
% deer      fitted DEER data
% bckg      fitted background
% param     fit parameters, structure with fields
%           .zt     zero time
%           .depth  modulation depth
%           .dim    background dimension
%           .kdec   background decay constant
%           .rmsd   fit rmsd
%
% G. Jeschke, 22.11.2018

texp = [];
vexp = [];
deer = [];
bckg = [];
param.zt = [];
param.depth = [];
param.dim = [];
param.kdec = [];
param.rmsd = [];

if ~exist('options','var') || isempty(options)
    options.zt = [];
end
if ~isfield(options,'depth')
    options.depth = [];
end
if ~isfield(options,'dim') || isempty(options.dim)
    options.dim = 3;
end
if ~isfield(options,'kdec')
    options.kdec = [];
end
if ~isfield(options,'cutoff') || isempty(options.cutoff)
    options.cutoff = 1;
end
if options.cutoff < 0.5
    options.cutoff = 0.5;
end
if options.cutoff > 1
    options.cutoff = 1;
end

[t0,~,data]=get_elexsys_MMM(dfile);

if isempty(t0) || isempty(data)
    add_msg_board(sprintf('ERROR: DEER data file %s could not be read.',dfile));
    return
end

load('pake_base40_MMM.mat');
Pake_r = r;
Pake_t = t;
Pake_kernel = base - ones(size(base));
scale = max(t0)/(1000*max(Pake_t));
rscale = scale^(1/3);
pcf=get_std_distr_MMM(rax/rscale,distr,Pake_r);
ff=get_formfactor_MMM(pcf,Pake_kernel,Pake_t);

[texp,vexp,deer,bckg,param] = fit_background(t0,data,scale*Pake_t,ff,options);

function [texp,vexp,deer,bckg,param] = fit_background(t,data,Pake_t,ff,options)

param.zt = [];
param.depth = [];
param.dim = [];
param.kdec = [];
param.rmsd = [];

[texp,vexp,zt] = pre_process_MMM(t,data); 

% cut off end artefact, if requested
pe = round(options.cutoff*length(texp));
texp = texp(1:pe)/1000;
vexp = vexp(1:pe);

param.zt = zt;

ff1=interp1(Pake_t,ff,texp,'pchip');
ff=ff1/max(ff1);

if length(options.kdec) == 1
    v = [options.kdec options.depth(1) options.dim(1)];
    [rmsd,deer,bckg,vexp] = rms_sim(v,texp,vexp,ff);
    param.kdec = v(1);
    param.depth = v(2);
    param.dim = v(3);
    param.rmsd = rmsd;
    return
end

constraints = zeros(3,2);
if length(options.kdec) == 2
    constraints(1,:) = options.kdec;
    if isnan(constraints(1,2))
        constraints(1,2) = 1e6;
    end
else
    constraints(1,1) = 0;
    constraints(1,2) = 1e6;
end

poly = polyfit(texp,log(vexp),1); % linear fit of logarithm
v(1) = -poly(1);

if length(options.depth) ~= 1
    v(2) = 1 - vexp(end);
end

if length(options.depth) == 2
    constraints(2,:) = options.depth;
elseif length(options.depth) == 1
    constraints(2,1) = options.depth;
    constraints(2,2) = options.depth;
else
    constraints(2,1) = 0.0001;
    constraints(2,2) = 0.9999;
end

if length(options.dim) == 2
    v(3) = mean(options.dim);
    constraints(3,:) = options.dim;
else
    constraints(3,1) = options.dim;
    constraints(3,2) = options.dim;
end

v1 = fminsearch(@rms_constrained_bckg,v,[],texp,vexp,ff,constraints);
param.kdec = v1(1);
v = zeros(1,3);
v(1) = param.kdec;
if length(v1) >= 2
    param.depth = v1(2);
    v(2) = param.depth;
else
    param.depth = options.depth;
    v(2) = param.depth;
end
if length(v1) >= 3
    param.dim = v1(3);
else
    param.dim = options.dim;
    v(3) = param.dim;
end
[rmsd,deer,bckg,vexp] = rms_sim(v,texp,vexp,ff);
param.rmsd = rmsd;

function [rms,sim,bckg,vexp] = rms_sim(v,texp,vexp,ff)

bckg = decay_stretched_MMM([v(1) v(3)],texp);
deconvoluted = v(2)*ff + (1-v(2))*ones(size(ff));
sim = deconvoluted.*bckg;
bckg = (1-v(2))*ones(size(ff)).*bckg;

sc = sum(sim.*sim)/sum(sim.*vexp);
vexp = sc*vexp;
diff=sim-vexp;
rms=sqrt(sum(diff.^2)/sum(vexp.^2));