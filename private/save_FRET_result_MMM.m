function handles=save_FRET_result_MMM(handles,bas_name,path)
%
% Saves a log file, listing the parameters of data analysis, the distance distribution 
% and, if any, the experimental transfer efficiency and simulated transfer efficiencies,
% and with RMSD if selected
% 
%
% Content and file name convention:
% - filenames are derived from bas_name
% - the method string is appended to the name of the data set
% - up to six ASCII output files are generated (depending on available data), 
%   which are distinguished by another
%   appendix to the file name:
%   _res    main parameters of the analysis and moment analysis of the
%           distribution
%   _distr  distance distribution, two columns, distance axis and
%           the simulated distribution
%   _AVdistr distance distribution from AV simulation
%   _Esim   simulated FRET transfer efficiencies with diffusion model if selected
%   _rmsd   if RMSD checkbox is selected also the RMSD data is saved in the 
%           same format as the Esim data
% 
% - files _distr, _Esim, and _rmsd have the extension '.dat', file _res has
%   the extension '.txt'
% - time axes are in microseconds, distance axes are in nanometers
%

    if exist('onCleanup','class')==8  % this is done only for Matlab version 
                                      % 2008a and later, where onCleanup exists
        c = onCleanup(@myCleanup);
    end

    fname_distr=[path bas_name '_distr.dat'];
    fname_AVdistr=[path bas_name '_AVdistr.dat'];
    fname_Esim=[path bas_name '_Esim.dat'];
    fname_rmsd=[path bas_name '_rmsd.dat'];
    fname_res=[path bas_name '_res.txt'];

    rmsd_flag=get(handles.checkbox_plot_rmsd,'Value');
    diff_flag=get(handles.checkbox_diffusion_model,'Value');
    AV_flag=get(handles.checkbox_display_AV,'Value');
    r=handles.rsim;
    simdistr=handles.dsim;
    
    AV_rax = [];
    AV_distr = [];
    if AV_flag
        AV_rax = handles.AV_rax;
        AV_distr = handles.AV_distr;
    end
    
%     texp=handles.texp_fit;
%     vexp=handles.vexp_fit;
% %     dipevo=handles.cluster;
%     sim=handles.vsim;
%     ffsim=handles.ff;
%     bckg=handles.bckg;

    
%     distr=handles.dexp_fit;
    

    % Normalize simulated dipolar evolution function, dipolar spectrum, and distance distribution
%     if ~isempty(bckg) && ~isempty(texp) && ~isempty(vexp) && ~ff_flag,
%         data2=[texp real(vexp) bckg sim];
%         save(fname_bckg,'data2','-ascii');
%     end;
%     if ~isempty(ffsim),
%         if ~isempty(vexp) && ~isempty(bckg),
%             cluster=real(vexp)./bckg;
%             cluster=cluster/max(cluster);
%             data3=[texp cluster ffsim];
%         else
%             data3=[handles.tsim ffsim];
%         end;
%         save(fname_fit,'data3','-ascii');
%     end;
    
    if ~isempty(simdistr)
        data1=[r' simdistr'];
        save(fname_distr,'data1','-ascii');
    end
    if AV_flag && ~isempty(AV_distr)
        data2=[AV_rax' AV_distr'];
        save(fname_AVdistr,'data2','-ascii');
    end

%     rms=handles.rmsd;

    % Determine analysis method
    deltaR0 = handles.F_radius.*[1-handles.error_F_radius/100 1+handles.error_F_radius/100]; %R0.*[0.93 1.07];
    nEffGrid = handles.R0_grid_points;
    if nEffGrid >= 2
        R0vec = linspace(deltaR0(1), deltaR0(2), nEffGrid);
    else
        R0vec = handles.F_radius;
    end
    
    if diff_flag %&& ~isempty(handles.FRET_efficiency_matrix)
        
        diffvec = linspace(handles.diff_const_lower_bound, handles.diff_const_upper_bound, handles.diffusion_grid_size);
        data3 = [NaN R0vec; diffvec' handles.FRET_efficiency_matrix'];
        
        str = '# First row: R0 values [nm]\n# First colunm: D values [nm^2/ns]';
        save(fname_Esim, 'str', '-ascii');
        save(fname_Esim, 'data3','-ascii','-append');
        
        if rmsd_flag
            data4 = [NaN R0vec; diffvec' handles.FRET_RMSD_matrix'];
            str = ['# First row: R0 values [nm]\n# First column: D values [nm^2/ns]\n# Data: RMSD values from exp. E_FRET ' num2str(handles.exp_FRET_efficiency)];
            save(fname_rmsd, 'str','-ascii');
            save(fname_rmsd, 'data4','-ascii','-append');
        end
        
    else
        data5 = [R0vec' handles.FRET_efficiency_static'];
        str = '# First column: R0 values [nm]\n# Second colunm: Static FRET efficiencies';
        save(fname_Esim, 'str','-ascii');
        save(fname_Esim, 'data5','-ascii','-append');
        
        if rmsd_flag
            data6 = [R0vec' handles.FRET_RMSD_static'];
            str = ['# First row: R0 values [nm]\n# Second column: RMSD values from exp. E_FRET ' num2str(handles.exp_FRET_efficiency)];
            save(fname_rmsd, 'str','-ascii');
            save(fname_rmsd, 'data6','-ascii','-append');
        end
    end
      
    
    if diff_flag
        format='FRET simulation with diffusional averaging calculation';
    else
        format='FRET static simulation only';
    end
   
    wfile=fopen(fname_res,'w+');
    fprintf(wfile,'%s%s%s\n','>>> MMM FRET simulation data set: ',bas_name,' <<<');
    fprintf(wfile,'%s%s%s\n\n','- ',format,' -');
    
    if handles.exp_FRET_efficiency ~= 0
        fprintf(wfile,'%s%4.3f%s\n\n','Experimental Transfer Efficiency <Eexp>: ',handles.exp_FRET_efficiency, ' ');
    else
        fprintf(wfile,'%s%2.1f%s\n\n','Simulated Transfer Efficiency <Eexp>: ',handles.exp_FRET_efficiency, ' (default)');
    end
    
    fprintf(wfile,'%s\n','FRET simulation parameters:');
    fprintf(wfile,'%s%4.2f%s\n','Foerster Radius (R0): ',handles.F_radius, ' nm');
    fprintf(wfile,'%s%4.1f%s\n','Error of R0: ',handles.error_F_radius, ' %');
    fprintf(wfile,'%s%d%s\n','R0 grid points: ',handles.R0_grid_points, ' ');
    fprintf(wfile,'%s%s%s\n\n','Simulated Transfer Efficiency <Esim>: ',handles.text_sim_FRET_efficiency.String, ' ');
    
    if diff_flag
        fprintf(wfile,'%s\n','Diffusion model parameters:');
        fprintf(wfile,'%s%5.3f%s%5.3f%s\n','Diffusion constant: ',handles.diff_const_lower_bound, ' - ', handles.diff_const_upper_bound, ' nm^2/ns');
        fprintf(wfile,'%s%5.2f%s\n','Excited donor lifetime: ',handles.exci_donor_lifetime, ' ns');
        fprintf(wfile,'%s%d%s\n','Diffusion grid size: ',handles.diffusion_grid_size, ' ');
        fprintf(wfile,'%s%d%s\n','Diffusion matrix size: ',handles.diff_matrix_size, ' ');
    end
    
    if rmsd_flag
        fprintf(wfile,'%s%s\n','Minimum RMSD found: ',handles.text_minimum_rmsd.String);
        fprintf(wfile,'%s%s%s\n','with R0: ',handles.text_best_R0.String, ' nm');
        if diff_flag
            fprintf(wfile,'%s%s%s\n','Diffusional averaging with D: ',handles.text_fitted_diffusion_constant.String, ' nm^2/ns');
        else
            fprintf(wfile,'%s\n','using Static averaging.');
        end
        fprintf(wfile,'%s%s\n\n','Resulting <Esim>: ',handles.text_sim_FRET_efficiency.String);
    end
    
    fprintf(wfile,'Suggested minimum Foerster Radius R0: %3.1f nm\n',handles.R0_min);
    fprintf(wfile,'Suggested optimum Foerster Radius R0: %3.1f nm\n',handles.R0_opt);
    fclose(wfile);

    function myCleanup

        if ~isempty(wfile==fopen('all'))
            fclose(wfile);
        end

    end
end