global model

trials = 5000000;
maxtime = 4*3600;
SANS_threshold = 2.5;

% reference points
ref(1) = 71;
ref(2) = 80;
ref(3) = 109;
ref(4) = 205;
ref(5) = 235;
ref(6) = 240;
ref(7) = 352;
ref(8) = 388;
ref(9) = 475;

if length(model.structures) == 1
    [snum,message]=copy_structure(1,'+PTB');
    fprintf(1,'Result is stored in structure %i with identifier +PTB\n',snum);
end

SANS_data = 'SANS_hPTB_100D2O_11m_clean.dat';

chain_A = model.structures{2}(1).xyz{1};
chain_B = model.structures{2}(2).xyz{1};
chain_C = model.structures{2}(3).xyz{1};
chain_D = model.structures{2}(4).xyz{1};
chain_E = model.structures{2}(5).xyz{1};
chain_F = model.structures{2}(6).xyz{1};
chain_G = model.structures{2}(7).xyz{1};

model.selected{1} = [2 1 1];
model.selected{2} = [2 2 1];
model.selected{3} = [2 3 1];
model.selected{4} = [2 4 1];
model.selected{5} = [2 5 1];
model.selected{6} = [2 6 1];
model.selected{7} = [2 7 1];

% message = wr_pdb_selected('curr_model','+PTB');
% pdbfile = 'curr_model';
% [chi,outname,status,result,fit] = fit_SANS_by_cryson(SANS_data,pdbfile);
% 
% figure(7); clf;
% plot(fit(:,1),fit(:,2),'k');
% hold on;
% plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);

% RNA gap definitions
maxr_per_nt = 7; % maximum length of extended RNA between consecutive P atoms

[msg,SLD3.coor] = get_object('(G)537.P','coor');
if msg.error
    fprintf(2,'Error reading RNA anchor coordinate. %s\n',msg.text);
end;
SLD3.nt = 21; % nucleotide number in construct
[msg,SLE5.coor] = get_object('(B)322.P','coor');
if msg.error
    fprintf(2,'Error reading RNA anchor coordinate. %s\n',msg.text);
end;
SLE5.nt = 39; % nucleotide number in construct
[msg,SLE3.coor] = get_object('(B)330.P','coor');
if msg.error
    fprintf(2,'Error reading RNA anchor coordinate. %s\n',msg.text);
end;
SLE3.nt = 47; % nucleotide number in construct
[msg,LinkEF5.coor] = get_object('(F)548.P','coor');
if msg.error
    fprintf(2,'Error reading RNA anchor coordinate. %s\n',msg.text);
end;
LinkEF5.nt = 57; % nucleotide number in construct
[msg,LinkEF3.coor] = get_object('(F)552.P','coor');
if msg.error
    fprintf(2,'Error reading RNA anchor coordinate. %s\n',msg.text);
end;
LinkEF3.nt = 61; % nucleotide number in construct
[msg,SLF5.coor] = get_object('(D)357.P','coor');
if msg.error
    fprintf(2,'Error reading RNA anchor coordinate. %s\n',msg.text);
end;
SLF5.nt = 74; % nucleotide number in construct

rna_gaps(1).rb1 = 3;
rna_gaps(1).coor1 = SLD3.coor;
rna_gaps(1).rb2 = 1;
rna_gaps(1).coor2 = SLE5.coor;
rna_gaps(1).length = SLE5.nt - SLD3.nt;

rna_gaps(2).rb1 = 1;
rna_gaps(2).coor1 = SLE3.coor;
rna_gaps(2).rb2 = 3;
rna_gaps(2).coor2 = LinkEF5.coor;
rna_gaps(2).length = LinkEF5.nt - SLE3.nt;

rna_gaps(3).rb1 = 3;
rna_gaps(3).coor1 = LinkEF3.coor;
rna_gaps(3).rb2 = 2;
rna_gaps(3).coor2 = SLF5.coor;
rna_gaps(3).length = SLF5.nt - LinkEF3.nt;

pmodel = 0.5;
pthr = exp(-erfinv(pmodel)^2);

forgive = 0.8;
min_approach = 5; % minimal approach of two reference points [Å]
max_extension = 180; % maximum distance between any two reference points [Å]
rad_ref=0.5; % radius for reference point sphere 
clash_threshold = 1.5*forgive;

restraints = load('nonahedron_restraints.dat');

data = load('rigid_body_protein_IAP_new.dat');
labels = zeros(9,3);
for k = 1:9
    labels(k,:) = data(k,2:4);
end;
lab_I = labels(1:3,:);
lab_II = labels(4:6,:);
lab_III = labels(7:9,:);

lab_II = [lab_II; 2.89, 21.28, 6.21]; % appends SLF U 15
lab_III = [lab_III; -44.04, 11.65, -12.64];

restrain_rna.r = 64.4;
restrain_rna.sigr = 5.8;

[msg,rrm1]=get_object('[PTBP](A)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading RRM1 coordinates');
    add_msg_board(msg.text);
    return
end
[msg,sle]=get_object('[PTBP](B)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading SLE coordinates');
    add_msg_board(msg.text);
    return
end
[msg,rrm2]=get_object('[PTBP](C)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading RRM2 coordinates');
    add_msg_board(msg.text);
    return
end
[msg,slf]=get_object('[PTBP](D)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading SLF coordinates');
    add_msg_board(msg.text);
    return
end
[msg,rrm34]=get_object('[PTBP](E)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading RRM34 coordinates');
    add_msg_board(msg.text);
    return
end
[msg,linkef]=get_object('[PTBP](F)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading Linker E/F coordinates');
    add_msg_board(msg.text);
    return
end
[msg,sld]=get_object('[PTBP](G)','xyz_heavy');
if msg.error,
    add_msg_board('ERROR in reading SLD coordinates');
    add_msg_board(msg.text);
    return
end
rigid_I = [rrm1;sle];
[ma_I,~] = size(rigid_I);
rigid_II = [rrm2;slf];
[ma_II,~] = size(rigid_II);
rigid_III = [rrm34;sld;linkef];
[ma_III,~] = size(rigid_III);


sig_fac = 1; % factor for multiplying standard deviation to make lower and upper bounds

dmat0 = zeros(9);
sigmat = zeros(9);
lb = zeros(9);
ub = zeros(9);

% translate restraint table to polyhedron indexing, make distance matrix,
% and lower and upper bound matrices
[m,~] = size(restraints);
for k = 1:m,
    [mi,poi1] = min(abs(restraints(k,1)-ref));
    if mi > 0
        fprintf(2,'Warning: Index 1 (i) not recognized for restraint %i\n',restraints(k,1),k);
    end
    restraints(k,1) = poi1;
    [mi,poi2] = min(abs(restraints(k,2)-ref));
    if mi > 0
        fprintf(2,'Warning: Index 2 (%i) not recognized for restraint %i\n',restraints(k,2),k);
    end
    restraints(k,2) = poi2;
    % conversion to Angstroem
    restraints(k,3:4) = restraints(k,3:4)*10;
    dmat0(poi1,poi2) = restraints(k,3); 
    sigmat(poi1,poi2) = restraints(k,4);
    lb(poi1,poi2) = restraints(k,3) - sig_fac*restraints(k,4);
    ub(poi1,poi2) = restraints(k,3) + sig_fac*restraints(k,4);
    dmat0(poi2,poi1) = dmat0(poi1,poi2);
    lb(poi2,poi1) = lb(poi1,poi2);
    ub(poi2,poi1) = ub(poi1,poi2);
    sigmat(poi2,poi1) = sigmat(poi1,poi2);
end

% augment lower and upper bounds
for k1 = 1:8
    for k2 = k1+1:9
        if lb(k1,k2) < min_approach
            lb(k1,k2) = min_approach;
            lb(k2,k1) = min_approach;
        end;
        if ub(k1,k2) < 0.1 % the unset bounds are zero
            ub(k1,k2) = max_extension;
            ub(k2,k1) = max_extension;
        end
    end
end

[lb,ub,err]=triangle(lb,ub,true);

switch err
    case 0
        fprintf(1,'Successful bound smoothing with experimental constraints.\n');
    case 1
        fprintf(2,'ERROR: Some distance restraints are inconsistent.\n');
    case 2
        fprintf(2,'ERROR: At least one bound matrix is not square.\n');
    case 3
        fprintf(2,'ERROR: Upper and lower bound matrices differ in size.\n');
    otherwise
        fprintf(2,'Unspecified error in bound smoothing.\n');
end;

met_err = 0;
embed_err = 0;
bound_err = 0;
clash_err = 0;
rna_fail = 0;
rna_link_fail = 0;
label_fail = 0;
sans_fail = 0;
success = 0;

nonahedra = cell(1,trials);

if trials <= 50000
    diagnostics = true;
else
    diagnostics = false;
end;
if diagnostics
    rna_dist = zeros(1,trials);
    prob_labels = zeros(1,trials);
    rms_trace = zeros(3,trials);
    rna_links = zeros(3,trials);
end;
chi_SANS = zeros(1,trials);
transmats = cell(1,3);

runtime = 0;
k = 0;
tic,
while runtime<=maxtime && k < trials
    k = k + 1;
    [dmatr,err]=metrize(lb,ub);
    if err==1, % metrization failed (restraints inconsistent), next trial, increment error counter
        met_err=met_err+1;
    else
        coor0=dmat2coor(dmatr); % embed distance matrix to Cartesian space
        if isempty(coor0)
            embed_err = embed_err + 1;
        else
            [coor1,err]=bound_refiner(coor0,lb,ub);
            if err > 0
                bound_err = bound_err + 1;
            else
                [rms1,lab_I_t,transmat_I] = rmsd_superimpose(coor1(1:3,:),lab_I(1:3,:));
                [rms2,lab_II_t,transmat_II] = rmsd_superimpose(coor1(4:6,:),lab_II(1:3,:));
                [rms3,lab_III_t,transmat_III] = rmsd_superimpose(coor1(7:9,:),lab_III(1:3,:));
                transmats{1} = transmat_I;
                transmats{2} = transmat_II;
                transmats{3} = transmat_III;   
                if diagnostics
                    rms_trace(1,k) = rms1;
                    rms_trace(2,k) = rms2;
                    rms_trace(3,k) = rms3;
                end
                t_labels = [lab_I_t;lab_II_t;lab_III_t];
                fulfill = true;
                lab_I_f = affine_coor_set(lab_I,transmat_I);
                lab_II_f = affine_coor_set(lab_II,transmat_II);
                lab_III_f = affine_coor_set(lab_III,transmat_III);
                rna_57_78 = norm(lab_II_f(4,:)-lab_III_f(4,:));
                if diagnostics
                    rna_dist(k) = rna_57_78;
                end;
                prob = prob_Gaussian(rna_57_78,restrain_rna.r,restrain_rna.sigr);
                if prob < pthr
                    rna_fail = rna_fail + 1;
                    fulfill = false;
                end;
                % check how well label-to-label distances are reproduced
                plabels = 1;
                for kr = 1:length(restraints)
                    rsim = norm(t_labels(restraints(kr,1),:) - t_labels(restraints(kr,2),:));
                    plabels = plabels*prob_Gaussian(rsim,restraints(kr,3),restraints(kr,4));
                    if diagnostics
                        prob_labels(k) = plabels;
                    end
                end;
                if fulfill
                    if plabels < pthr^9
                        fulfill = false;
                        label_fail = label_fail + 1;
                    end;
                end;
                % check for RNA linker lengths
                if fulfill
                    for kl = 1:length(rna_gaps)
                        a3 = affine_coor_set(rna_gaps(kl).coor1,transmats{rna_gaps(kl).rb1});
                        a5 = affine_coor_set(rna_gaps(kl).coor2,transmats{rna_gaps(kl).rb2});
                        rlink = norm(a5-a3);
                        if diagnostics
                            rna_links(kl,k) = rlink;
                        end
                        if rlink > maxr_per_nt*rna_gaps(kl).length
                            fulfill = false;
                        end
                    end
                    if ~fulfill,
                        rna_link_fail = rna_link_fail + 1;
                    end;
                end
                if fulfill
                    rigid_Ib = affine_coor_set(rigid_I,transmat_I);
                    rigid_IIb = affine_coor_set(rigid_II,transmat_II);
                    mdist = min_dist(rigid_Ib,rigid_IIb);
                    if mdist < clash_threshold
                        clash_err = clash_err + 1;
                    else
                        rigid_IIIb = affine_coor_set(rigid_III,transmat_III);
                        mdist = min_dist(rigid_Ib,rigid_IIIb);
                        if mdist < clash_threshold
                            clash_err = clash_err + 1;
                        else
                            mdist = min_dist(rigid_IIb,rigid_IIIb);
                            if mdist < clash_threshold
                                clash_err = clash_err + 1;
                            else
                                success = success + 1;
                                nonahedra{success} = [transmat_I;transmat_II;transmat_III];
                                t_chain_A = affine_coor_set(chain_A,transmat_I);
                                t_chain_B = affine_coor_set(chain_B,transmat_I);
                                t_chain_C = affine_coor_set(chain_C,transmat_II);
                                t_chain_D = affine_coor_set(chain_D,transmat_II);
                                t_chain_E = affine_coor_set(chain_E,transmat_III);
                                t_chain_F = affine_coor_set(chain_F,transmat_III);
                                t_chain_G = affine_coor_set(chain_G,transmat_III);
                                [~,~,tmstd] = rmsd_superimpose(chain_C,t_chain_C);
                                t2_chain_A = affine_coor_set(t_chain_A,tmstd);
                                t2_chain_B = affine_coor_set(t_chain_B,tmstd);
                                t2_chain_C = affine_coor_set(t_chain_C,tmstd);
                                t2_chain_D = affine_coor_set(t_chain_D,tmstd);
                                t2_chain_E = affine_coor_set(t_chain_E,tmstd);
                                t2_chain_F = affine_coor_set(t_chain_F,tmstd);
                                t2_chain_G = affine_coor_set(t_chain_G,tmstd);
                                model.structures{2}(1).xyz{success} = t2_chain_A;
                                model.structures{2}(2).xyz{success} = t2_chain_B;
                                model.structures{2}(3).xyz{success} = t2_chain_C;
                                model.structures{2}(4).xyz{success} = t2_chain_D;
                                model.structures{2}(5).xyz{success} = t2_chain_E;
                                model.structures{2}(6).xyz{success} = t2_chain_F;
                                model.structures{2}(7).xyz{success} = t2_chain_G;
                                if success > 1,
                                    for kc = 1:7,
                                        model.structures{2}(kc).residues{success} = model.structures{2}(kc).residues{1};
                                        model.structures{2}(kc).Bfactor{success} = model.structures{2}(kc).Bfactor{1};
                                        model.structures{2}(kc).Btensor{success} = model.structures{2}(kc).Btensor{1};                                        
                                    end;
                                end;
                                % select only the RRMs of the current model
                                model = rmfield(model,'selected');
                                model.selected{1} = [2 1 success];
                                model.selected{2} = [2 3 success];
                                model.selected{3} = [2 5 success];
                                pdbfile = sprintf('t%i',k);
                                message = wr_pdb_selected(pdbfile,'+PTB');
                                [chi,outname,status,result,fit] = fit_SANS_by_cryson(SANS_data,pdbfile);
                                if ~isempty(chi)
                                    chi_SANS(success) = chi;
                                    if chi < SANS_threshold
                                        figure(10+success); clf;
                                        plot(fit(:,1),fit(:,2),'k');
                                        hold on;
                                        plot(fit(:,1),fit(:,3),'Color',[0.75,0,0]);
                                        title(sprintf('SANS fit for model %i',success));
                                        drawnow
                                    else
                                        fprintf(1,'SANS fit above threshold (chi = %5.2f)\n',chi);
                                        success = success - 1;
                                        sans_fail = sans_fail + 1;
                                    end;
                                else
                                    fprintf(2,'Warning: SANS fitting failed in trial %i:\n',k);
                                    fprintf(2,'%s',result);
                                    success = success - 1;
                                    sans_fail = sans_fail + 1;
                                end;
                            end
                        end
                    end
                end;
            end;
        end;
    end;
    runtime=toc;
end;
toc,

trials = k;

chi_SANS = chi_SANS(1:success);
figure(100); clf;
plot(chi_SANS,'o');

if success > 0
    model = rmfield(model,'selected');
    spoi = 0;
    for kc = 1:7,
        for km = 1:success,
            spoi = spoi + 1;
            model.selected{spoi} = [2 kc km];
        end;
    end;
    message = wr_pdb_selected('PTBP1_model','+PTB');
end

left_trials = trials;
fprintf(1,'%5.2f%% of %i trials were successful. %i models were obtained.\n',100*success/trials,trials,success);
fprintf(1,'%5.2f%% of all trials failed in metrization.\n',100*met_err/trials);
left_trials = left_trials - met_err;
fprintf(1,'%5.2f%% of all trials (%i) failed in embedding (%5.2f%% of still active trials).\n',100*embed_err/trials,embed_err,100*embed_err/left_trials);
left_trials = left_trials - embed_err;
fprintf(1,'%5.2f%% of all trials (%i) failed on bound refinement (%5.2f%% of still active trials).\n',100*bound_err/trials,bound_err,100*bound_err/left_trials);
left_trials = left_trials - bound_err;
fprintf(1,'%5.2f%% of all trials (%i) failed on RNA restraint (%5.2f%% of still active trials).\n',100*rna_fail/trials,rna_fail,100*rna_fail/left_trials);
left_trials = left_trials - rna_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed on label restraints (%5.2f%% of still active trials).\n',100*label_fail/trials,label_fail,100*label_fail/left_trials);
left_trials = left_trials - label_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed on RNA linker length (%5.2f%% of still active trials).\n',100*rna_link_fail/trials,rna_link_fail,100*rna_link_fail/left_trials);
left_trials = left_trials - rna_link_fail;
fprintf(1,'%5.2f%% of all trials (%i) failed in clash tests (%5.2f%% of still active trials).\n',100*clash_err/trials,clash_err,100*clash_err/left_trials);
left_trials = left_trials - clash_err;
fprintf(1,'%5.2f%% of all trials (%i) failed due to poor SANS fit (%5.2f%% of still active trials).\n',100*sans_fail/trials,sans_fail,100*sans_fail/left_trials);
left_trials = left_trials - sans_fail;
if left_trials ~= success
    fprintf(2,'Trial dissipation. Expected success: %i. Found success %i.\n',left_trials,success);
end

if diagnostics
    figure(777); clf;
    plot(rna_dist,'k');
    hold on;
    plot([1,trials],[restrain_rna.r,restrain_rna.r],'g');
    plot([1,trials],[restrain_rna.r-restrain_rna.sigr,restrain_rna.r-restrain_rna.sigr],'r:');
    plot([1,trials],[restrain_rna.r+restrain_rna.sigr,restrain_rna.r+restrain_rna.sigr],'r:');
    
    figure(333); clf;
    plot(rms_trace(1,:));
    hold on;
    plot(rms_trace(2,:));
    plot(rms_trace(3,:));
    
    figure(111); clf;
    plot(prob_labels,'k');
    hold on;
    plot([1,trials],[pthr^9,pthr^9],'r:');
    
    cols = [0.6,0,0;0,0.4,0;0,0,0.7];
    figure(555); clf;
    hold on;
    for kl = 1:length(rna_gaps)
        plot(rna_links(kl,:),'.','Color',cols(kl,:));
        plot([1,trials],[maxr_per_nt*rna_gaps(kl).length,maxr_per_nt*rna_gaps(kl).length],':','Color',cols(kl,:));
    end;
end;

