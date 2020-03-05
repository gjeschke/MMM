function [included,populations,bsl] = fit_ensemble(restraints,model_list,options)
% [included,populations] = fit_ensemble(restraints,model_list,options)
%
% Fits the populations of model conformations so that they best fit a set
% of DEER and small-angle scattering restraints
%
% restraints    unprocessed RigiFlex restraints, to get them, use:
%               rd_restraints_rigiflex(restraint_file,unprocessed)
% model_list    list of file names for model conformations, should not be
%               longer than about 100
% options       options for restraint evaluation, see
%               evaluate_model_restraints.m
% 
% included      indices of models with significant population (>1% of the
%               maximum population), corresponds to model_list
% populations   populations of the included models
% bsl           fitted SANS/SAXS curve baseline offsets
%
% G. Jeschke, 28.6.2019


% ### Code for testing the caller without actually fitting ###
% included = 1:length(model_list);
% coeff1 = rand(1,length(model_list)).^20;
% coeff = coeff1/max(coeff1);
% included = included(coeff >= options.threshold);
% coeff = coeff(coeff >= options.threshold);
% populations = coeff/sum(coeff);
% add_msg_board(sprintf('%i populations are at least 1%% of the maximum population\n',length(coeff))); 
% return

bsl = [];

DEER_scores = zeros(1,length(model_list));
SAS_scores = DEER_scores;

initialized = false;
ensemble_size = length(model_list);

if ~isfield(options,'interactive') || isempty(options.interactive)
    options.interactive = true;
end
if options.interactive
    options.text_status.String = 'Evaluation';
    options.text_info_status.String = sprintf('of %i',ensemble_size);
end

included = zeros(1,ensemble_size);
populations = [];
len_DEER = 0;
for km = 1:length(model_list)
    if options.interactive
        options.text_iteration_counter.String = sprintf('%i',km);
    end
    drawnow
    fname = sprintf('%s_restraints.mat',model_list{km});
    if options.interactive
        fprintf(1,'%s\n',model_list{km});
    end
    if exist(fname,'file')
        s = load(fname);
        model_restraints = s.restraints;
    else
        model_restraints = evaluate_model_restraints(model_list{km},restraints,options);
    end
    [overlaps,DEER_fits] = get_DEER(model_restraints,restraints);
%     if isempty(DEER_fits)
%         disp('Aber hallo!');
%     end
    [chi2_vec,SAS_fits,valid] = get_SAS(model_restraints,restraints);
    if isempty(overlaps) && ~valid
        continue
    end
    if isempty(DEER_fits) && ~isempty(restraints.DEER(1).r)
        continue
    end
    n_DEER = length(overlaps);
    DEER_scores(km) = 1 - prod(overlaps)^(1/n_DEER);
    SAS_scores(km) = sum(chi2_vec);
    if ~initialized
        initialized = true;
        len_DEER = length(DEER_fits);
        mpoi = 1;
        all_DEER_fits = cell(1,len_DEER);
        for k = 1:len_DEER
            dat0 = DEER_fits{k};
            [m,~] = size(dat0);
            data = zeros(m,2+ensemble_size);
            data(:,1:3) = dat0(:,1:3);
            all_DEER_fits{k} = data;
        end
        all_SAS_fits = cell(1,length(SAS_fits));
        for k = 1:length(SAS_fits)
            dat0 = SAS_fits{k};
            [m,~] = size(dat0);
            data = zeros(m,3+ensemble_size);
            data(:,1:2) = dat0(:,1:2);
            data(:,3) = dat0(:,4);
            data(:,4) = dat0(:,3);
            all_SAS_fits{k} = data;
        end
    else
        mpoi = mpoi + 1;
        for k = 1:len_DEER
            dat0 = DEER_fits{k};
            data = all_DEER_fits{k};
            data(:,2+mpoi) = dat0(:,3);
            all_DEER_fits{k} = data;
        end
        for k = 1:length(SAS_fits)
            dat0 = SAS_fits{k};
            data = all_SAS_fits{k};
            data(:,3+mpoi) = dat0(:,3);
            all_SAS_fits{k} = data;
        end
    end
    included(mpoi) = km;
end

if ~initialized
    add_msg_board('ERROR: For none of the models all DEER and SAS restraints could be evaluated.');
    included = [];
    populations = [];
    return
end

% Reduce data arrays to the number of models for which all DEER and SAS
% restraints could be evaluated
included = included(1:mpoi);
included0 = included;

for k = 1:len_DEER
    data = all_DEER_fits{k};
    data = data(:,1:2+mpoi);
    all_DEER_fits{k} = data;
end

for k = 1:length(SAS_fits)
    data = all_SAS_fits{k};
    data = data(:,1:3+mpoi);
    all_SAS_fits{k} = data;
end

if options.interactive
    options.text_status.String = 'Iteration';
    options.text_info_status.String = '· 1000';
end

ensemble_size = mpoi;

if len_DEER > 0

    clear fit_multi_DEER % initialize iteration counter

    fanonym_DEER = @(v_opt)fit_multi_DEER(v_opt,all_DEER_fits,options);

    l_DEER = zeros(1,ensemble_size);
    u_DEER = ones(1,ensemble_size);
    v0 = ones(1,ensemble_size)/ensemble_size;

    fit_options = optimoptions('patternsearch',...
        'MaxFunctionEvaluations',10000*length(v0),'MaxIterations',500*length(v0),...
        'StepTolerance',options.threshold/10);

    tic,
    if isfield(restraints,'no_fitting') && restraints.no_fitting 
        fom_DEER = fit_multi_DEER(v0,all_DEER_fits,options);
        v = v0;
        exitflag = -1;
        fit_output.iterations= 0;
        fit_output.funccount = 0;
        fit_output.meshsize = NaN;
        fit_output.maxconstraint = NaN;
    else
        [v,fom_DEER,exitflag,fit_output] = patternsearch(fanonym_DEER,v0,[],[],[],[],l_DEER,u_DEER,[],fit_options);
    end
    DEER_time = toc;
    th = floor(DEER_time/3600);
    DEER_time = DEER_time - 3600*th;
    tmin = floor(DEER_time/60);
    ts = round(DEER_time - 60*tmin);

    if isfield(restraints,'no_fitting') && restraints.no_fitting
        add_msg_board('DEER restraints were not fitted.');
    else
        add_msg_board(sprintf('DEER restraints were fitted in %i h %i min %i s with overlap deficiency of %5.3f\n',th,tmin,ts,fom_DEER));
    end
    if options.interactive
        options.text_DEER_fom0.String = sprintf('%5.3f',fom_DEER);
    end

    log_fid = fopen(options.logfile,'at');
    if log_fid ~= -1 % there is currently no error reporting here, if the logfile cannot be written
        fprintf(log_fid,'DEER restraints fit took %i h %i min %i s',th,tmin,ts);
        switch exitflag
            case 0
                fprintf(log_fid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
            case 1
                fprintf(log_fid,'Convergence by mesh size.\n');
            case 2
                fprintf(log_fid,'Convergence by change in population.\n');
            case 3
                fprintf(log_fid,'Convergence by DEER overlap deficiency precision.\n');
            case 4
                fprintf(log_fid,'Convergence by machine precision.\n');
            case -1
                fprintf(log_fid,'ERROR: Optimization terminated by output or plot function.\n');
            case -2
                fprintf(log_fid,'ERROR: No feasible solution found.\n');
        end
        fprintf(log_fid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
        fprintf(log_fid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
        if options.interactive
            if exitflag >=0
                esize = sum(v >= options.threshold*max(v));
                options.text_DEER_size.String = sprintf('%i',esize);
                options.text_DEER_fom0.String = sprintf('%5.3f',fom_DEER);
                fprintf(log_fid,'DEER-only ensemble has %i conformers and overlap deficiency %5.3f\n\n',esize,fom_DEER);
            else
                options.text_DEER_size.String = sprintf('%i',NaN);
                options.text_DEER_fom0.String = sprintf('%5.3f',NaN);
            end
        end
        fprintf(log_fid,'\n');
    end
    fclose(log_fid);
    normalize.DEER = fom_DEER;
    coeff1 = v;
    coeff = coeff1/max(coeff1);
    included = included(coeff >= options.threshold);
    coeff = coeff(coeff >= options.threshold);
    populations = coeff/sum(coeff);
    add_msg_board(sprintf('%i populations are at least 1%% of the maximum population for DEER only\n',length(coeff)));
end

if ~isempty(SAS_fits)
    included = included0;
    clear fit_multi_SAS % initialize iteration counter

    fanonym_SAS = @(v_opt)fit_multi_SAS(v_opt,all_SAS_fits,options);

    l_SAS = [zeros(1,ensemble_size) -0.01*ones(1,length(all_SAS_fits))];
    u_SAS = [ones(1,ensemble_size) 0.01*ones(1,length(all_SAS_fits))];
    v0 = [ones(1,ensemble_size)/ensemble_size zeros(1,length(all_SAS_fits))];
    
    fit_options = optimoptions('patternsearch',...
        'MaxFunctionEvaluations',10000*length(v0),'MaxIterations',500*length(v0),...
        'StepTolerance',options.threshold/10);

    tic,
    if isfield(restraints,'no_fitting') && restraints.no_fitting 
        fom_SAS = fit_multi_SAS(v0,all_SAS_fits,options);
        v = v0;
        exitflag = -1;
        fit_output.iterations= 0;
        fit_output.funccount = 0;
        fit_output.meshsize = NaN;
        fit_output.maxconstraint = NaN;
    else
        [v,fom_SAS,exitflag,fit_output] = patternsearch(fanonym_SAS,v0,[],[],[],[],l_SAS,u_SAS,[],fit_options);
    end
    SAS_time = toc;
    th = floor(SAS_time/3600);
    SAS_time = SAS_time - 3600*th;
    tmin = floor(SAS_time/60);
    ts = round(SAS_time - 60*tmin);

    bsl = v(end-length(all_SAS_fits)+1:end);

     if isfield(restraints,'no_fitting') && restraints.no_fitting 
         add_msg_board('SAS restraints were not fitted.');
     else
         add_msg_board(sprintf('SAS restraints were fitted in %i h %i min %i s with sum of chi^2 of %5.2f',th,tmin,ts,fom_SAS));
     end
    log_fid = fopen(options.logfile,'at');
    if log_fid ~= -1 % there is currently no error reporting here, if the logfile cannot be written
        fprintf(log_fid,'SAXS/SANS curves fit took %i h %i min %i s\n',th,tmin,ts);
        switch exitflag
            case 0
                fprintf(log_fid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
            case 1
                fprintf(log_fid,'Convergence by mesh size.\n');
            case 2
                fprintf(log_fid,'Convergence by change in population.\n');
            case 3
                fprintf(log_fid,'Convergence by SAS chi^2 sum precision.\n');
            case 4
                fprintf(log_fid,'Convergence by machine precision.\n');
            case -1
                fprintf(log_fid,'ERROR: Optimization terminated by output or plot function.\n');
            case -2
                fprintf(log_fid,'ERROR: No feasible solution found.\n');
        end
        fprintf(log_fid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
        fprintf(log_fid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
        if exitflag >=0
            esize = sum(v >= options.threshold*max(v));
            if options.interactive
                options.text_SAS_size.String = sprintf('%i',esize);
                options.text_SAS_fom0.String = sprintf('%5.3f',fom_SAS);
            end
            fprintf(log_fid,'SAS-only ensemble has %i conformers and chi^2 sum %5.3f\n\n',esize,fom_SAS);
            fprintf(log_fid,'Fitted SAS-curve baseline offsets:\n');
            bpoi = 0;
            if isfield(model_restraints,'SANS') && ~isempty(model_restraints.SANS)
                for ks = 1:length(model_restraints.SANS)
                    bpoi = bpoi + 1;
                    fprintf(log_fid,'%s: %6.4f\n',model_restraints.SANS(ks).data,bsl(bpoi));
                end
            end
            if isfield(model_restraints,'SAXS') && ~isempty(model_restraints.SAXS)
                for ks = 1:length(model_restraints.SAXS)
                    bpoi = bpoi + 1;
                    fprintf(log_fid,'%s: %6.4f\n',model_restraints.SAXS(ks).data,bsl(bpoi));
                end
            end
        elseif options.interactive
            options.text_SAS_size.String = sprintf('%i',NaN);
            options.text_SAS_fom0.String = sprintf('%5.3f',NaN);
        end
        fprintf(log_fid,'\n');
    end
    fclose(log_fid);
    coeff1 = v(1:end-length(all_SAS_fits));
    bsl = v(end-length(all_SAS_fits)+1:end);
    coeff = coeff1/max(coeff1);
    included = included(coeff >= options.threshold);
    coeff = coeff(coeff >= options.threshold);
    populations = coeff/sum(coeff);
    add_msg_board(sprintf('%i populations are at least 1%% of the maximum population for SAS fitting\n',length(coeff))); 
    normalize.SAS = fom_SAS;
end

if len_DEER > 0 && ~isempty(SAS_fits)
    included = included0;
    clear fit_multi_SAS_DEER % initialize iteration counter

    fanonym_integrative = @(v_opt)fit_multi_SAS_DEER(v_opt,all_SAS_fits,all_DEER_fits,normalize,options);

    l = [zeros(1,ensemble_size) -0.01*ones(1,length(all_SAS_fits))];
    u = [ones(1,ensemble_size) 0.01*ones(1,length(all_SAS_fits))];
    v0 = [ones(1,ensemble_size)/ensemble_size zeros(1,length(all_SAS_fits))];
    fit_options = optimoptions('patternsearch',...
        'MaxFunctionEvaluations',50000*length(v0),'MaxIterations',2500*length(v0),...
        'StepTolerance',options.threshold/10);

    tic,
    if isfield(restraints,'no_fitting') && restraints.no_fitting 
        fom = fit_multi_SAS_DEER(v0,all_SAS_fits,all_DEER_fits,normalize,options);
        v = v0;
        exitflag = -1;
        fit_output.iterations= 0;
        fit_output.funccount = 0;
        fit_output.meshsize = NaN;
        fit_output.maxconstraint = NaN;
    else
        [v,fom,exitflag,fit_output] = patternsearch(fanonym_integrative,v0,[],[],[],[],l,u,[],fit_options);
    end
    integrative_time = toc;



    th = floor(integrative_time/3600);
    integrative_time = integrative_time - 3600*th;
    tmin = floor(integrative_time/60);
    ts = round(integrative_time - 60*tmin);

    if isfield(restraints,'no_fitting') && restraints.no_fitting 
        add_msg_board('No integrative fitting performed.');
    else
        add_msg_board(sprintf('Integrative fit in %i h %i min %i s with loss of %5.3f',th,tmin,ts,fom));
    end

    coeff1 = v(1:end-length(all_SAS_fits));
    bsl = v(end-length(all_SAS_fits)+1:end);

    fom_DEER = sim_multi_DEER(coeff1,all_DEER_fits);
    [fom_SAS,fits_SAS] = sim_multi_SAS(v,all_SAS_fits);

    if options.interactive
        for ksas = 1:length(fits_SAS)
            figure(3000+ksas); clf; hold on;
            plot(fits_SAS(ksas).s,fits_SAS(ksas).curve,'k');
            plot(fits_SAS(ksas).s,fits_SAS(ksas).sim,'r');
            title(sprintf('chi2 = %5.2f',fits_SAS(ksas).chi2));
        end
    end

    coeff = coeff1/max(coeff1);
    included = included(coeff >= options.threshold);
    coeff = coeff(coeff >= options.threshold);
    populations = coeff/sum(coeff);
    add_msg_board(sprintf('%i populations are at least 1%% of the maximum population\n',length(coeff))); 

    log_fid = fopen(options.logfile,'at');
    if log_fid ~= -1 % there is currently no error reporting here, if the logfile cannot be written
        fprintf(log_fid,'Integrative ensemble fit took %i h %i min %i s\n',th,tmin,ts);
        switch exitflag
            case 0
                fprintf(log_fid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
            case 1
                fprintf(log_fid,'Convergence by mesh size.\n');
            case 2
                fprintf(log_fid,'Convergence by change in population.\n');
            case 3
                fprintf(log_fid,'Convergence by loss precision.\n');
            case 4
                fprintf(log_fid,'Convergence by machine precision.\n');
            case -1
                fprintf(log_fid,'ERROR: Optimization terminated by output or plot function.\n');
            case -2
                fprintf(log_fid,'ERROR: No feasible solution found.\n');
        end
        fprintf(log_fid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
        fprintf(log_fid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
        if exitflag >=0
            esize = sum(v >= options.threshold*max(v));
            if options.interactive
                options.text_ensemble_size.String = sprintf('%i',esize);
                options.text_loss.String = sprintf('%5.3f',fom);
                options.text_DEER_fom.String = sprintf('%5.3f',fom_DEER);
                options.text_SAS_fom.String = sprintf('%5.3f',fom_SAS);
            end
            fprintf(log_fid,'Integrative ensemble has %i conformers and loss %5.3f\n\n',length(included),fom);
            fprintf(log_fid,'Fitted SAS-curve baseline offsets:\n');
            bpoi = 0;
            if isfield(model_restraints,'SANS') && ~isempty(model_restraints.SANS)
                for ks = 1:length(model_restraints.SANS)
                    bpoi = bpoi + 1;
                    fprintf(log_fid,'%s: %6.4f\n',model_restraints.SANS(ks).data,bsl(bpoi));
                end
            end
            if isfield(model_restraints,'SAXS') && ~isempty(model_restraints.SAXS)
                for ks = 1:length(model_restraints.SAXS)
                    bpoi = bpoi + 1;
                    fprintf(log_fid,'%s: %6.4f\n',model_restraints.SAXS(ks).data,bsl(bpoi));
                end
            end
        elseif options.interactive
            options.text_ensemble_size.String = sprintf('%i',NaN);
            options.text_loss.String = sprintf('%5.3f',NaN);
            options.text_DEER_fom.String = sprintf('%5.3f',NaN);
            options.text_SAS_fom.String = sprintf('%5.3f',NaN);
        end
        if exitflag >=0
        end
        fprintf(log_fid,'\n');
    end
    fclose(log_fid);
end



function [overlaps,fits] = get_DEER(restraints,restraints0)



% make DEER distance distributions

core = length(restraints.DEER);
flex = 0;

if isfield(restraints,'pflex') && ~isempty(restraints.pflex) && isfield(restraints.pflex,'DEER')...
        && isfield(restraints0,'pflex') && ~isempty(restraints0.pflex) && isfield(restraints0.pflex,'DEER')
    for kl = 1:length(restraints.pflex)
        flex = flex + length(restraints.pflex(kl).DEER);
    end
end

overlaps = zeros(1,core+flex);
fits = cell(1,core+flex);



n_DEER = 0;
score_DEER = 1;
DEER_valid = true;

for k = 1:min(length(restraints.DEER),length(restraints0.DEER))
    if ~isempty(restraints0.DEER(k).r) && restraints0.DEER(k).r ~=0 &&  restraints0.DEER(k).sigr ~=0  % this is a restraint to be tested
        rax = restraints.DEER(1).rax;
        distr = restraints.DEER(k).distr;
        n_DEER = n_DEER + 1;
        if ~isempty(distr)
            distr = distr/sum(distr);
            data = zeros(length(rax),3);
            data(:,1) = rax.';
            argr = (restraints.DEER(k).r-rax)/(sqrt(2)*restraints.DEER(k).sigr);
            distr_sim = exp(-argr.^2);
            distr_sim = distr_sim/sum(distr_sim);
            if isfield(restraints.DEER(k),'distr_exp') && ~isempty(restraints.DEER(k).distr_exp)
                data(:,2) = restraints.DEER(k).distr_exp.';
                overlaps(n_DEER) = sum(min([restraints.DEER(k).distr_exp;distr]));
            else
                data(:,2) = distr_sim.';
                overlaps(n_DEER) = sum(min([distr_sim;distr]));
            end
            data(:,3) = distr.';
            fits{n_DEER} = data;            
            score_DEER = score_DEER * overlaps(n_DEER);
        else
            DEER_valid = false;
            break
        end
    end
end

if ~DEER_valid
    overlaps = [];
    fits = {};
    return
end


pflex = core;

if isfield(restraints,'pflex') && ~isempty(restraints.pflex) && isfield(restraints.pflex,'DEER')...
        && isfield(restraints0,'pflex') && ~isempty(restraints0.pflex) && isfield(restraints0.pflex,'DEER')
    for kl = 1:length(restraints.pflex)
        for k = 1:length(restraints.pflex(kl).DEER)
            rax = restraints.pflex(kl).DEER(k).rax;
            pflex = pflex + 1;
            if restraints.pflex(kl).DEER(k).r ~=0 &&  restraints.pflex(kl).DEER(k).sigr ~=0 % valid restraint
                n_DEER = n_DEER + 1;
                distr = restraints.pflex(kl).DEER(k).distr;
                if ~isempty(distr)
                    distr = distr/sum(distr);
                    data = zeros(length(rax),3);
                    data(:,1) = rax.';
                    argr = (restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*restraints.pflex(kl).DEER(k).sigr);
                    distr_sim = exp(-argr.^2);
                    distr_sim = distr_sim/sum(distr_sim);
                    if isfield(restraints.pflex(kl).DEER(k),'distr_exp') && ~isempty(restraints.pflex(kl).DEER(k).distr_exp)
                        data(:,2) = restraints.pflex(kl).DEER(k).distr_exp.';
                        overlaps(n_DEER) = sum(min([restraints.pflex(kl).DEER(k).distr_exp;distr]));
                    else
                        data(:,2) = distr_sim.';
                        overlaps(n_DEER) = sum(min([distr_sim;distr]));
                    end
                    data(:,3) = distr.';
                    fits{n_DEER} = data;
                    overlaps(n_DEER) = sum(min([distr_sim;distr]));
                    score_DEER = score_DEER * overlaps(n_DEER);
                else
                    DEER_valid = false;
                    break
                end
            end
        end
    end
end
% score_DEER = 1 - score_DEER^(1/n_DEER);
% fprintf(1,'Total DEER score : %5.3f\n',score_DEER);

if ~DEER_valid
    overlaps = [];
    fits = {};
    return
end

overlaps = real(overlaps(1:n_DEER));
fits = fits(1:n_DEER);

function [chi2_vec,fits,valid] = get_SAS(restraints,restraints0)

chi2_vec = zeros(1,10);
fits = cell(1,10);
valid = true;

n_SANS = 0;
if isfield(restraints,'SANS') && ~isempty(restraints.SANS) && isfield(restraints0,'SANS') && ~isempty(restraints0.SANS)
    n_SANS = length(restraints.SANS);
    for ks = 1:n_SANS
        fits{ks} = restraints.SANS(ks).fit;
        chi2_vec(ks) = restraints.SANS(ks).chi2;
    end
end

n_SAXS = 0;
if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS) && isfield(restraints0,'SAXS') && ~isempty(restraints0.SAXS)
    n_SAXS = length(restraints.SAXS);
    for ks = 1:n_SAXS
        fits{n_SANS+ks} = restraints.SAXS(ks).fit;
        chi2_vec(n_SANS+ks) = restraints.SAXS(ks).chi2;
    end
end

chi2_vec = chi2_vec(1:n_SANS+n_SAXS);
fits = fits(1:n_SANS+n_SAXS);


