function diagnostics = rigiflex_assembler(restraints,options,handles)
% Assemble rigid-body models and flexible-section models to complete models
% 
% function diagnostics = rigiflex_assembler(restraints,options,handles)
%
% restraints    restraint definition for the complete model
% options       running options
%               .trials     number of assembled models to be tested
%               .snum_rb    structure index of the rigid-body arrangements
%               .snum_flex  array of structure indices for all flexible
%                           sections of all rigid-body arrangements
% handles       handles to the RigiFlex GUI, can be missing for modelling
%               on a remote serverhandles.text_time_left

global model

diagnostics.success = 0;

% if handles of the GUI figure are supplied, the engine will update the MMM
% GUI during the computation, else this is considered as a run on a remote
% server that does not have access to the GUI
if exist('handles','var')
    interactive = true;
else
    interactive = false;
end;

[rigid_models,flex_sections] = size(options.snum_flex);

% determine how many models are available for the flexible sections
flex_models = zeros(rigid_models,flex_sections);
for kr = 1:rigid_models
    for kf = 1:flex_sections
        flex_models(kr,kf) = length(model.structures{options.snum_flex(kr,kf)}(1).residues);
    end
end

trials = options.trials;
if prod(prod(flex_models)) < trials
    trials = prod(prod(flex_models));
end;

poir = 0;
poif = zeros(rigid_models,flex_sections);
fnums = zeros(1,flex_sections);
increments = primes(1000000);
% select combinations that are likely to differ and to use all available
% models
for trial = 1:trials
    poir = poir + 1;
    if poir > rigid_models
        poir = 1;
    end
    for fs = 1:flex_sections
        pf = rem(poif(poir,fs)+increments(fs),flex_models(poir,fs));
        if pf == 0, pf = flex_models(poir,fs); end;
        fnums(fs) = pf;
        poif(poir,fs) = pf;
    end
end;