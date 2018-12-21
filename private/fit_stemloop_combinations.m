function [solutions,trafo,libs] = fit_stemloop_combinations(restraints,snum0,snum,repname)

global model

forgive = 0.8;
clash_threshold = 1.5*forgive;


link_nt = zeros(1,length(restraints.RNA.bind)-1);
for k = 2:length(restraints.RNA.bind)
    link_nt(k-1) = restraints.RNA.bind(k).nta-restraints.RNA.bind(k-1).nte;
end
thr_length_per_nt = 6;
max_length_per_nt = 5.5;
max_link_length = link_nt*thr_length_per_nt; 

nmod = length(model.structures{snum}(1).residues);
trafo = cell(nmod,3);
template = cell(1,3);
target = cell(1,3);

stag = mk_address(snum);
stag0 = mk_address(snum0);

libs = cell(1,length(restraints.stemlibs));
chains = char(1,3);
sls = char(1,3);
for k = 1:length(restraints.stemlibs)
    lib = load(restraints.stemlibs{k}.name);
    libs{k} = lib.library;
    [~,xyz] = get_object(sprintf('%s%s{1}',stag0,restraints.stemlibs{k}.rrm),'xyz_heavy');
    template{k} = xyz;
    chains(k) = restraints.stemlibs{k}.rrm(2);
    sls(k) = restraints.stemlibs{k}.chaintag(2);
end

report = fopen(repname,'wt');
% for kl = 1:3
%     figure(10+kl); clf;
%     hold on;
%     plot([1,nmod],[max_link_length(kl),max_link_length(kl)],'k');
% end
tic,

solutions = zeros(10000,5);
valid_rbas = zeros(1,nmod);
spoi = 0;
for k = 1:nmod % 1:nmod loop over all models
    % determine the transformations
    for kc = 1:3
        adr = sprintf('%s(%c){%i}',stag,chains(kc),k);
        [~,xyz_new] = get_object(adr,'xyz_heavy');
        target{kc} = xyz_new;
        [~,~,transmat] = rmsd_superimpose(xyz_new,template{kc});
        trafo{k,kc} = transmat;
    end
    % find non-clashing RNA stemloops
    valid_decoys = zeros(3,100);
    vpoi = zeros(1,3);
    for ksl = 1:3
        library = libs{ksl};
        for kdecoy = 1:length(library.chains)
            xyz_sl = library.chains{kdecoy}.xyz{1};
            transmat = trafo{k,ksl};
            xyz_sl = affine_trafo_coor(xyz_sl,transmat);
            full_cost = 0;
            all_costs = zeros(1,3);
            for kc = 1:3
                xyz_rrm = target{kc};
                cost = clash_cost(xyz_rrm,xyz_sl,clash_threshold);
                all_costs(kc) = cost;
                if kc ~= ksl % exclude clashes with the own RRM
                    full_cost = full_cost + cost;
                end
            end
            if full_cost < 5*eps
                vpoi(ksl) = vpoi(ksl) + 1;
                valid_decoys(ksl,vpoi(ksl)) = kdecoy;
                % fprintf(report,'Model %i, stemloop %c(%i) has clash cost %6.2f (%6.2f, %6.2f, %6.2f)\n',k,sls(ksl),kdecoy,full_cost,all_costs);
            end
        end
    end
    fprintf(report,'Model %i, %i SL%c models, %i SL%c models, and %i SL%c models are non-clashing.\n',k,vpoi(1),sls(1),vpoi(2),sls(2),vpoi(3),sls(3));
    % check link restraints
    transmat1 = trafo{k,1};
    transmat2 = trafo{k,2};
    transmat3 = trafo{k,3};
    library = libs{1};
    anchor_12 = library.linksites(1).coor;
    anchor_12 = affine_trafo_coor(anchor_12,transmat1);
    library = libs{2};
    anchor_21 = library.linksites(1).coor;
    anchor_21 = affine_trafo_coor(anchor_21,transmat2);
    anchor_23 = library.linksites(2).coor;
    anchor_23 = affine_trafo_coor(anchor_23,transmat2);
    library = libs{3};
    anchor_43 = library.linksites(1).coor;
    anchor_43 = affine_trafo_coor(anchor_43,transmat3);
    [~,a32c] = get_object(sprintf('%s(F){%i}342.C5''',stag,k),'coor');
    [~,a34c] = get_object(sprintf('%s(F){%i}344.C5''',stag,k),'coor');
    r12max = link_nt(1)*thr_length_per_nt;
    r23max = link_nt(2)*thr_length_per_nt;
    r34max = link_nt(3)*thr_length_per_nt;
    min_stretch = 1e6;
    best_r = [1000,1000,1000];
    best_combi = [0,0,0];
    % loop over all non-clashing decoys
    for k1 = 1:vpoi(1)
        a12c = anchor_12(valid_decoys(1,k1),:);
%         fprintf(report,'Anchor 12(%i): %4.1f, %4.1f, %4.1f Å\n',valid_decoys(1,k1),a12c);
        for k2 = 1:vpoi(2)
            a21c = anchor_21(valid_decoys(2,k2),:);
            r12 = norm(a12c-a21c);
%             figure(11);
%             plot(k,r12,'r.');
            a23c = anchor_23(valid_decoys(2,k2),:);
            r23 = norm(a23c-a32c);
%             if k1 == 1
%                 fprintf(report,'Anchor 21(%i): %4.1f, %4.1f, %4.1f Å\n',valid_decoys(2,k2),a21c);
%                 fprintf(report,'Anchor 23(%i): %4.1f, %4.1f, %4.1f Å\n',valid_decoys(2,k2),a23c);
%             end
%             if k1 == 1
%                 figure(12);
%                 plot(k,r23,'g.');
%             end
            for k3 = 1:vpoi(3)
                a43c = anchor_43(valid_decoys(3,k3),:);
                r34 = norm(a34c-a43c);
%                 if k1 == 1 && k2 == 1
%                     figure(13);
%                     plot(k,r34,'b.');
%                 end
%                 if k1 == 1 && k2 == 1
%                     fprintf(report,'Anchor 43(%i): %4.1f, %4.1f, %4.1f Å\n',valid_decoys(3,k3),a43c);
%                 end
                stretch = 0;
                lnt12 = r12/link_nt(1);
                if lnt12 > max_length_per_nt
                    stretch = stretch + (lnt12-max_length_per_nt)^2;
                end
                lnt23 = r23/link_nt(2);
                if lnt23 > max_length_per_nt
                    stretch = stretch + (lnt23-max_length_per_nt)^2;
                end
                lnt34 = r34/link_nt(3);
                if lnt34 > max_length_per_nt
                    stretch = stretch + (lnt34-max_length_per_nt)^2;
                end
                if stretch < min_stretch
                    min_stretch = stretch;
                    best_r = [r12,r23,r34];
                    best_combi = [valid_decoys(1,k1),valid_decoys(2,k2),valid_decoys(3,k3)];
                end
                if r12 <= r12max && r23 <= r23max && r34 <= r34max
                    fprintf(report,'M(%i)SL(%i,%i,%i): %5.1f, %5.1f, %5.1f Å\n',k,...
                        valid_decoys(1,k1),valid_decoys(2,k2),valid_decoys(3,k3),...
                        r12,r23,r34);
                    spoi = spoi + 1;
                    solutions(spoi,1) = k;
                    solutions(spoi,2) = valid_decoys(1,k1);
                    solutions(spoi,3) = valid_decoys(2,k2);
                    solutions(spoi,4) = valid_decoys(3,k3);
                    valid_rbas(k) = 1;
                    solutions(spoi,5) = stretch;
                end
            end
        end
    end
    fprintf(report,'Best combination:\n, M(%i)SL(%i,%i,%i): %5.1f, %5.1f, %5.1f Å\n',k,best_combi,best_r);
end
toc
solutions = solutions(1:spoi,:);
save test_SL solutions trafo 
fprintf(report,'\n%i solutions were found in %i valid RBAs.\n',spoi,sum(valid_rbas));
fclose(report);
add_msg_board(sprintf('%i solutions were found in %i valid RBAs.\n',spoi,sum(valid_rbas)));
