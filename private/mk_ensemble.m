function [score_DEER,sum_chi2] = mk_ensemble(ensemble_list,populations,restraints,options)
% mk_ensemble(ensemble_list,populations,restraints,options)
%
% Write ensemble file and make restraint fulfillment plots
%
% ensemble_list     list of conformer PDB files making up the ensemble
% populations       vector of populations, defaults to uniform populations
% restraints        restraint definition
% options           .script     flag indicating whether an MMM display
%                               script is to be written, default: true
%                   .write_pdb  flag indicating whether ensemble file is to
%                               be written out, default: true
%                   .ordering   ensemble ordering, defaults to sequence in
%                               file list
%                   .individual flag indicating whether distance
%                               distributions of individual conformers are
%                               plotted, defaults to false
%                   .groups     groups of ensemble members, for which a
%                               sum distance distribution is plotted,
%                               array of struct with fields
%                               .conformers  conformer numbers
%                               .color       RGB color triple
%                               defaults to empty
%
% G. Jeschke, 7.1.2020



global model


if ~exist('populations','var') || isempty(populations)
    populations = ones(1,length(ensemble_list));
    populations = populations/sum(populations);
end

if ~exist('options','var') || isempty(options) || ~isfield(options,'script') || isempty(options.script)
    options.script = true;
end

if ~isfield(options,'write_pdb') || isempty(options.write_pdb)
    options.write_pdb = true;
end

if ~isfield(options,'ordering') || isempty(options.ordering)
    options.ordering = 1:length(ensemble_list);
end


if ~isfield(options,'plot') || isempty(options.plot)
    options.plot = true;
end

populations = populations(options.ordering);

if ~isfield(options,'individual') || isempty(options.individual)
    options.individual = false;
end

if ~isfield(options,'groups') || isempty(options.groups)
    options.groups = [];
end

DEER.initialize = true;

if options.script
    fid = fopen(options.script_file,'wt');
    fprintf(fid,'show [%s] ribbon\n',options.pdb_ID);
    fprintf(fid,'color [%s](B){:} grey\n',options.pdb_ID);
    fprintf(fid,'colorscheme [%s](A){:} sequence\n',options.pdb_ID);

    for km = 1:length(populations)
        alpha = populations(km)/max(populations);
        fprintf(fid,'transparency [%s](:){%i} %5.3f\n',options.pdb_ID,km,alpha);
    end
    fclose(fid);
end

for km = 1:length(ensemble_list)
    fname = sprintf('%s_restraints.mat',ensemble_list{options.ordering(km)});
    if exist(fname,'file')
        s = load(fname);
        model_restraints = s.restraints;
    else
        model_restraints = evaluate_model_restraints(ensemble_list{km},restraints,options);
    end
    % fprintf(1,'%s: %6.3f\n',fname,populations(km));
    DEER = display_DEER_single(km,model_restraints,populations,DEER,options);
    if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
        for ks = 1:length(restraints.SANS)
            if km == 1
                restraints.SANS(ks).fit = model_restraints.SANS(ks).fit;
                restraints.SANS(ks).fit(:,3) = populations(1)*model_restraints.SANS(ks).fit(:,3);
            else
                restraints.SANS(ks).fit(:,3) = restraints.SANS(ks).fit(:,3) +...
                    populations(km)*model_restraints.SANS(ks).fit(:,3);
            end
        end
    end
    if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS)
        for ks = 1:length(restraints.SAXS)
            if km == 1
                restraints.SAXS(ks).fit = model_restraints.SAXS(ks).fit;
                restraints.SAXS(ks).fit(:,3) = populations(1)*model_restraints.SAXS(ks).fit(:,3);
            else
                restraints.SAXS(ks).fit(:,3) = restraints.SAXS(ks).fit(:,3) +...
                    populations(km)*model_restraints.SAXS(ks).fit(:,3);
            end
        end
    end
end

sum_chi2 = 0;
bpoi = 0;
if isfield(restraints,'SANS') && ~isempty(restraints.SANS)
    for ks = 1:length(restraints.SANS)
        bpoi = bpoi + 1;
        [fit,chi2] = scale_offset_SAS_fit(restraints.SANS(ks).fit);
        sum_chi2 = sum_chi2 + chi2;
        if options.plot
            figure(10000+bpoi); clf;
            plot(fit(:,1),fit(:,2),'k');
            hold on
            plot(fit(:,1),fit(:,3),'Color',[0.6,0,0]);
            title(sprintf('SANS %s: chi^2 %5.2f\n',restraints.SANS(ks).data,chi2));
            xlabel('Scattering vector s');
            ylabel('Scattering amplitude [n.a.]');
        end
    end
end
if isfield(restraints,'SAXS') && ~isempty(restraints.SAXS)
    for ks = 1:length(restraints.SAXS)
        bpoi = bpoi + 1;
        [fit,chi2] = scale_offset_SAS_fit(restraints.SAXS(ks).fit);
        sum_chi2 = sum_chi2 + chi2;
        if options.plot
            figure(10000+bpoi); clf;
            plot(fit(:,1),fit(:,2),'k');
            hold on
            plot(fit(:,1),fit(:,3),'Color',[0.6,0,0]);
            title(sprintf('SAXS %s: chi^2 %5.2f\n',restraints.SAXS(ks).data,chi2));
            xlabel('Scattering vector s');
            ylabel('Scattering amplitude [n.a.]');
        end
    end
end

score_DEER = 1;
n_DEER = 0;

for k = 1:DEER.core+DEER.flex
    DEER.all_distr(k,:) = DEER.all_distr(k,:)/sum(DEER.all_distr(k,:));
    if options.plot
        figure(k); hold on;
    end
    for kgrp = 1:length(options.groups)
        DEER.groups_distr{kgrp}(k,:) = DEER.groups_distr{kgrp}(k,:)/sum(DEER.all_distr(k,:));
        if options.plot
            plot(DEER.rax,DEER.groups_distr{kgrp}(k,:),'Color',options.groups(kgrp).color);
        end
    end
    if options.plot
        plot(DEER.rax,DEER.all_distr(k,:),'Color',[0.6,0,0]);
    end
    axis_vec = DEER.axis_vecs(k,:);
    if 1.05*max(DEER.all_distr(k,:)) > axis_vec(4)
        axis_vec(4) = 1.05*max(DEER.all_distr(k,:));
    end
    overlap_exp = NaN;
    if sum(DEER.all_distr_exp(k,:)) > 10*eps
        overlap_exp = sum(min([DEER.all_distr_exp(k,:);DEER.all_distr(k,:)]));
        ma = max(DEER.all_distr_exp(k,:));
        if 1.05*ma > axis_vec(4)
            axis_vec(4) = 1.05*ma;
        end
    end
    overlap = sum(min([DEER.all_distr_sim(k,:);DEER.all_distr(k,:)]));
    if ~isnan(overlap) && overlap > 10*eps && DEER.all_flags(k)
        score_DEER = score_DEER*overlap;
        n_DEER = n_DEER + 1;
    end
    adr1 = DEER.all_adr1{k};
    adr2 = DEER.all_adr2{k};
    if options.plot
        if DEER.all_flags(k)
            ftype = 'fit';
            title(sprintf('%s-%s. o_{exp}: %5.3f, mod.depth %5.3f',adr1,adr2,overlap_exp,DEER.mod_depths(k)));
            % title(sprintf('%s-%s. o_{res}: %5.3f, o_{exp}: %5.3f, mod.depth %5.3f',adr1,adr2,overlap,overlap_exp,DEER.mod_depths(k)));
            % Restraint plotting in model was disabled
            % colr = get_color(overlap);
            % fprintf(fid,'plot %s.CA %s.CA 6 %5.3f %5.3f %5.3f :\n',adr1,adr2,colr);
        else
            ftype = 'control';
            title(sprintf('%s-%s. o_{exp}: %5.3f, mod.depth %5.3f',adr1,adr2,overlap_exp,DEER.mod_depths(k)));
            % Restraint plotting in model was disabled
            % colr = get_color(overlap_exp);
            % fprintf(fid,'plot %s.CA %s.CA 6 %5.3f %5.3f %5.3f :\n',adr1,adr2,colr);
        end
        axis(axis_vec);
    end
    [texp,vexp,deer,bckg,param] = fit_DEER_primary(DEER.rax/10,DEER.all_distr(k,:),strcat('deer\',DEER.all_fnames{k}),options);
    if ~isempty(texp) && options.plot
        figure(1000+k); clf;
        plot(texp,vexp,'k');
        hold on;
        plot(options.cutoff*[max(texp),max(texp)],[min(vexp),max(vexp)],'b');
        plot(texp,deer,'Color',[0.75,0,0]);
        plot(texp,bckg,'Color',[0,0.6,0]);
        title(sprintf('%s: %s-%s. rmsd: %6.4f',ftype,adr1,adr2,param.rmsd));
    end
end

score_DEER = 1 - score_DEER^(1/n_DEER);

add_msg_board(sprintf('Mean DEER overlap deficiency: %6.3f',score_DEER));

add_msg_board(sprintf('Sum of all SANS/SAXS chi^2: %5.2f',sum_chi2));

if options.write_pdb
    % Make ensemble PDB file
    [msg,snum0] = add_pdb(ensemble_list{1});
    if msg.error
        add_msg_board(sprintf('ERROR: First PDB file %s could not be read (%s).',ensemble_list{1},msg.text));
        return
    end
    for kf = 2:length(ensemble_list)
        copy_structure(snum0,'+mod',[],kf,snum0); % make a copy of the first model
        [msg,snumc] = add_pdb(ensemble_list{kf});
        if msg.error
            add_msg_board(sprintf('ERROR: PDB file #%i (%s) could not be read (%s).',kf,ensemble_list{kf},msg.text));
            return
        end
        for kc = 1:length(model.structures{snum0})
            model.structures{snum0}(kc).atoms{kf} = model.structures{snumc}(kc).atoms{1};
            model.structures{snum0}(kc).residues{kf} =model.structures{snumc}(kc).residues{1};
            model.structures{snum0}(kc).xyz{kf} = model.structures{snumc}(kc).xyz{1};
            model.structures{snum0}(kc).Bfactor{kf} =  model.structures{snumc}(kc).Bfactor{1};
            model.structures{snum0}(kc).Btensor{kf} = model.structures{snumc}(kc).Btensor{1};
        end
    end
    model.selected = cell(1,1);
    model.selected{1} = snum0;
    message = wr_pdb_selected(options.pdb_file,options.pdb_ID);

    if message.error
        add_msg_board(strcat('ERROR (wr_pdb_selected): ',message.text));
    else
        add_msg_board(sprintf('Structure file %s written',options.pdb_file));
    end
end



function DEER = display_DEER_single(km,model_restraints,populations,DEER,options)

core = length(model_restraints.DEER);
rax = model_restraints.DEER(1).rax;
flex = 0;

if DEER.initialize
    DEER.all_distr = zeros(core,length(rax));
    for k = 1:length(options.groups)
        DEER.groups_distr{k} = zeros(core,length(rax));
    end
    DEER.all_distr_sim = zeros(core,length(rax));
    DEER.all_distr_exp = zeros(core,length(rax));
    DEER.axis_vecs = zeros(core,4);
    DEER.all_fnames = cell(1,core);
    DEER.all_adr1 = cell(1,core);
    DEER.all_adr2 = cell(1,core);
    DEER.all_flags = zeros(1,core);
    DEER.mod_depths = zeros(1,core);
    DEER.core = core;
    DEER.flex = flex;
    DEER.rax = rax;
end


if isfield(model_restraints,'pflex') && ~isempty(model_restraints.pflex) && isfield(model_restraints.pflex,'DEER')
    for kl = 1:length(model_restraints.pflex)
        flex = flex + length(model_restraints.pflex(kl).DEER);
    end
    if DEER.initialize
        DEER.all_distr = zeros(core+flex,length(rax));
        for k = 1:length(options.groups)
            DEER.groups_distr{k} = zeros(core+flex,length(rax));
        end
        DEER.all_distr_sim = zeros(core+flex,length(rax));
        DEER.all_distr_exp = zeros(core+flex,length(rax));
        DEER.axis_vecs = zeros(core+flex,4);
        DEER.all_fnames = cell(1,core+flex);
        DEER.all_adr1 = cell(1,core+flex);
        DEER.all_adr2 = cell(1,core+flex);
        DEER.all_flags = zeros(1,core+flex);
        DEER.mod_depths = zeros(1,core+flex);
        DEER.core = core;
        DEER.flex = flex;
        DEER.rax = rax;
    end
end

for k = 1:length(model_restraints.DEER)
    DEER.all_adr1{k} = model_restraints.DEER(k).adr1;
    DEER.all_adr2{k} = model_restraints.DEER(k).adr2;
    distr = model_restraints.DEER(k).distr;
    if options.plot
        figure(k);
        if km == 1
            clf; hold on;
            title(sprintf('%s-%s',model_restraints.DEER(k).adr1,model_restraints.DEER(k).adr2));
        end
    end
    axis_vec = [0,120,-1e-6,1e-6];
    if isfield(model_restraints.DEER(k),'file') && ~isempty(model_restraints.DEER(k).file)
        dfname = strcat(model_restraints.DEER(k).file,'_distr.dat');
        deer_basname = strcat('DeerLab\',model_restraints.DEER(k).file);
        if km == 1
            [axis_vec,md] = plot_exp_dist(deer_basname,rax,options);
            DEER.mod_depths(k) = md;
        end
        dfname = strcat('DeerLab\',dfname);
        Pdata = load(dfname);
        rexp = 10*Pdata(:,1).';
        distr_exp = Pdata(:,2).';
        distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
        distr_exp = distr_exp/sum(distr_exp);
        DEER.all_distr_exp(k,:) = distr_exp;
    end
    if ~isempty(distr)
        argr = (model_restraints.DEER(k).r-rax)/(sqrt(2)*model_restraints.DEER(k).sigr);
        distr_sim = exp(-argr.^2);
        distr_sim = distr_sim/sum(distr_sim);
        DEER.all_distr_sim(k,:) = distr_sim;
        if km == 1 && options.plot
            plot(rax,distr_sim,'Color',[0,0.6,0]);
        end
        if options.individual && options.plot
            plot(rax,populations(km)*distr,'Color',[0.2,0.2,1]);
        end
        DEER.all_distr(k,:) = DEER.all_distr(k,:) + populations(km)*distr;
        for kgrp = 1:length(options.groups)
            if min(abs(options.groups(kgrp).conformers-km)) == 0
                DEER.groups_distr{kgrp}(k,:) = DEER.groups_distr{kgrp}(k,:) + populations(km)*distr;
            end
        end
        if isfield(model_restraints.DEER(k),'file')
            DEER.all_fnames{k} = model_restraints.DEER(k).file;
        end
        if model_restraints.DEER(k).r ~=0 &&  model_restraints.DEER(k).sigr ~=0
            DEER.all_flags(k) = 1;
            if 1.05*max(distr_sim) > axis_vec(4)
                axis_vec(4) = 1.05*max(distr_sim);
            end
        end
    end
    if km == 1
        DEER.axis_vecs(k,:) = axis_vec;
    end
end

pflex = core;
if isfield(model_restraints,'pflex') && ~isempty(model_restraints.pflex) && isfield(model_restraints.pflex,'DEER')
    for kl = 1:length(model_restraints.pflex)
        for k = 1:length(model_restraints.pflex(kl).DEER)
            pflex = pflex + 1;
            if options.plot
                figure(pflex);
                if km == 1
                    clf; hold on;
                    title(sprintf('%s-%s',model_restraints.DEER(k).adr1,model_restraints.DEER(k).adr2));
                end
            end
            distr = model_restraints.pflex(kl).DEER(k).distr;
            axis_vec = [0,120,-1e-6,1e-6];
            if isfield(model_restraints.pflex(kl).DEER(k),'file') && ~isempty(model_restraints.pflex(kl).DEER(k).file)
                dfname = strcat(model_restraints.pflex(kl).DEER(k).file,'_distr.dat');
                deer_basname = strcat('DeerLab\',model_restraints.pflex(kl).DEER(k).file);
                if km == 1
                    [axis_vec,md] = plot_exp_dist(deer_basname,rax,options);
                    DEER.mod_depths(pflex) = md;
                end
                dfname = strcat('DeerLab\',dfname);
                Pdata = load(dfname);
                rexp = 10*Pdata(:,1).';
                distr_exp = Pdata(:,2).';
                distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
                distr_exp = distr_exp/sum(distr_exp);
                DEER.all_distr_exp(pflex,:) = distr_exp;
            end
            if ~isempty(distr)
                if options.plot
                    title(sprintf('%s-%s',model_restraints.pflex(kl).DEER(k).adr1,model_restraints.pflex(kl).DEER(k).adr2));
                end
                argr = (model_restraints.pflex(kl).DEER(k).r-rax)/(sqrt(2)*model_restraints.pflex(kl).DEER(k).sigr);
                distr_sim = exp(-argr.^2);
                distr_sim = distr_sim/sum(distr_sim);
                DEER.all_distr_sim(pflex,:) = distr_sim;
                if km == 1 && options.plot
                    plot(rax,distr_sim,'Color',[0,0.6,0]);
                end
                if options.individual && options.plot
                    plot(rax,populations(km)*distr,'Color',[0.2,0.2,1]);
                end
                DEER.all_distr(pflex,:) = DEER.all_distr(pflex,:) + populations(km)*distr;
                for kgrp = 1:length(options.groups)
                    if min(abs(options.groups(kgrp).conformers-km)) == 0
                        DEER.groups_distr{kgrp}(pflex,:) = DEER.groups_distr{kgrp}(pflex,:) + populations(km)*distr;
                    end
                end
                if isfield(model_restraints.pflex(kl).DEER(k),'file')
                    DEER.all_fnames{pflex} = model_restraints.pflex(kl).DEER(k).file;
                end
                if model_restraints.pflex(kl).DEER(k).r ~=0 &&  model_restraints.pflex(kl).DEER(k).sigr ~=0
                    DEER.all_flags(pflex) = 1;
                    if 1.05*max(distr_sim) > axis_vec(4)
                        axis_vec(4) = 1.05*max(distr_sim);
                    end
                end
                DEER.all_adr1{pflex} = model_restraints.pflex(kl).DEER(k).adr1;
                DEER.all_adr2{pflex} = model_restraints.pflex(kl).DEER(k).adr2;
            end
            if km == 1
                DEER.axis_vecs(pflex,:) = axis_vec;
            end
        end
    end
end

DEER.initialize = false;

function [fit,chi2] = scale_offset_SAS_fit(fit)

curve = fit(:,2);
sim = fit(:,3);
errors = fit(:,4);
sc = sum(curve.*sim)/sum(sim.*sim);

v0 = [sc,0];

[v,chi2] = fminsearch(@chi2_scale_offset,v0,[],curve,sim,errors);
fit(:,3) = v(1)*sim+v(2);

function chi2 = chi2_scale_offset(v,curve,sim,errors)

m = length(curve);
chi2 = sum(((curve-v(1)*sim-v(2))./errors).^2)/(m-1);


function [axis_vec,md] = plot_exp_dist(deer_basname,rax,options)

rmax = 120;

data = load(sprintf('%s_distr.dat',deer_basname));
rexp = 10*data(:,1);
distr_exp = data(:,2);
distr_lb = data(:,3);
distr_ub = data(:,4);
distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
scal = sum(distr_exp);
distr_exp = distr_exp/scal;
distr_lb = interp1(rexp,distr_lb,rax,'pchip',0);
distr_lb = distr_lb/scal;
distr_ub = interp1(rexp,distr_ub,rax,'pchip',0);
distr_ub = distr_ub/scal;
if options.plot
    fill([rax'; flipud(rax')],[distr_ub'; flipud(distr_lb')],0.75*[1,1,1],'LineStyle','none')
    plot(rax,distr_exp,'k','LineWidth',1.5);
    grid on
    box on
end
ma = max(distr_ub);
axis_vec = [0,rmax,-0.1*ma,1.05*ma];
if options.plot
    axis(axis_vec);
    xlabel('distance r [?]');
    ylabel('P(r)');
end

data = load(sprintf('%s_bckg.dat',deer_basname));
md = 1-max(data(:,3));
if options.plot
    title(sprintf('m.d. %4.2f',md));
end
data = load(sprintf('%s_fit.dat',deer_basname));
tmax = data(end,1);
sc = (tmax/2)^(1/3);
r0 = 15;
r1 = 30*sc;
skip = false;
if r1 > rmax
    r1 = rmax;
    skip = true;
end
if options.plot
    plot([r0,r1],[-0.04*ma,-0.04*ma],'Color',[0,0.6,0],'LineWidth',4);
end
if ~skip
    r2 = 40*sc;
    if r2 > rmax
        r2 = rmax;
        skip = true;
    end
    if options.plot
        plot([r1,r2],[-0.04*ma,-0.04*ma],'Color',[0.8,0.8,0],'LineWidth',4);
    end
end
if ~skip
    r3 = 50*sc;
    if r3 > rmax
        r3 = rmax;
        skip = true;
    end
    if options.plot
        plot([r2,r3],[-0.04*ma,-0.04*ma],'Color',[0.8,0.6,0],'LineWidth',4);
    end
end
if ~skip
    r4 = 60*sc;
    if r4 > rmax
        r4 = rmax;
    end
    if options.plot
        plot([r3,r4],[-0.04*ma,-0.04*ma],'Color',[0.6,0,0],'LineWidth',4);
    end
end

% function rgb = get_color(overlap,limit)
% 
% if ~exist('limit','var')
%     limit = 0.8;
% end
% sigdim = 0.4;
% red = [1 0 0];
% green = [0 1 0];
% if overlap > limit
%     overlap = limit;
% end
% dimit = exp(-(overlap-limit/2)^2/(2*sigdim^2));
% fulfill = overlap/limit;
% rgb = dimit*(fulfill*green + (1-fulfill)*red);