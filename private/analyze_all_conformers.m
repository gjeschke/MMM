function analyze_all_conformers(listfile)

model_list = get_file_list(listfile);
scores = zeros(1,length(model_list));
sas_chi2 = zeros(1,length(model_list));
for k = 1:length(model_list)
    [~,~,score_DDR,chi2_vec] = analyze_conformer(model_list{k});
    scores(k) = score_DDR;
    sas_chi2(k) = sum(chi2_vec);
end

fprintf(1,'Mean DDR overlap deficiency        : %5.3f\n',mean(scores));
fprintf(1,'Std. dev. of DDR overlap deficiency: %5.3f\n',std(scores));
fprintf(1,'Mean SAS chi^2                     : %5.3f\n',mean(sas_chi2));
fprintf(1,'Std. dev. of SAS chi^2             : %5.3f\n',std(sas_chi2));

