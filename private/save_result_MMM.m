function handles=save_result_MMM(handles,bas_name,path)
%
% Saves a log file, listing the parameters of data analysis 
% the original data and background data, if any,
% the experimental (if any) and simulated dipolar form factors,
% and the distance distribution
%
% the format is similar to DeerAnalysis output format
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
%   _bckg   fit of primary data, three columns, time axis, primary
%           experimental data (real part), background function, simulated
%           primary data, this file exists only if there were experimental
%           data and only if primary data (as opposed to the form factor)
%           were fitted
%   _fit    fit of the dipolar evolution function, up to three columns,
%           time axis, experimental form factor, simulated
%           form factor, if only two columns exist, the
%           experimental form factor is missing
%   _distr  distance distribution, up to three columns, distance axis and
%           experimental distribution and simulated distribution, if only
%           two columns exist, the experimental distribution is missing
% - files _distr, _bckg, and _fit have the extension '.dat', file _res has
%   the extension '.txt'
% - time axes are in microseconds, distance axes are in nanometers
%

    if exist('onCleanup','class')==8, % this is done only for Matlab version 
                                      % 2008a and later, where onCleanup exists
        c = onCleanup(@myCleanup);
    end;

    fname_distr=[path bas_name '_distr.dat'];
    fname_distr_anyrotamers=[path bas_name '_anyrotamers_distr.dat'];
    fname_bckg=[path bas_name '_bckg.dat'];
    fname_fit=[path bas_name '_fit.dat'];
    fname_res=[path bas_name '_res.txt'];

    ff_flag=get(handles.checkbox_form_factor,'Value');
    texp=handles.texp_fit;
    vexp=handles.vexp_fit;
    dipevo=handles.cluster;
    sim=handles.vsim;
    ffsim=handles.ff;
    bckg=handles.bckg;
    r=handles.rsim;
    distr=handles.dexp_fit;
    simdistr=handles.dsim;

    % Normalize simulated dipolar evolution function, dipolar spectrum, and distance distribution
    if ~isempty(bckg) && ~isempty(texp) && ~isempty(vexp) && ~ff_flag,
        data2=[texp real(vexp) bckg sim];
        save(fname_bckg,'data2','-ascii');
    end;
    if ~isempty(ffsim),
        if ~isempty(vexp) && ~isempty(bckg),
            cluster=real(vexp)./bckg;
            cluster=cluster/max(cluster);
            data3=[texp cluster ffsim];
        else
            data3=[handles.tsim ffsim];
        end;
        save(fname_fit,'data3','-ascii');
    end;
    if ~isempty(simdistr),
        if ~isempty(distr),
            data1=[r' distr' simdistr'];
        else
            data1=[r' simdistr'];
        end;
        save(fname_distr,'data1','-ascii');
    end;
    if ~isempty(handles.tweak_distr)
        anyrotamers = interp1(handles.tweak_rax,handles.tweak_distr,r,'pchip',0);
        if ~isempty(distr),
            data4=[r' distr' anyrotamers'];
        else
            data4=[r' anyrotamers'];
        end
        save(fname_distr_anyrotamers,'data4','-ascii');
    end

    rms=handles.rmsd;

    % Determine analysis method
    fitted=1;
    format='Rotamer simulation and experimental data';
    if isempty(texp) || isempty(vexp) || isempty(bckg),
        format='Rotamer simulation only';
        fitted=0;
    end;

    contents = get(handles.popupmenu_bckg,'String');
    bckg_mode=contents{get(handles.popupmenu_bckg,'Value')};
    if strcmp(bckg_mode,'fractal'),
        bckg_mode=sprintf('%s with exponent %s',bckg_mode,get(handles.edit_bckg,'String'));
    end;


    wfile=fopen(fname_res,'w+');
    fprintf(wfile,'%s%s%s\n','>>> MMM DEER rotamer simulation of data set: ',bas_name,' <<<');
    fprintf(wfile,'%s%s%s\n\n','- ',format,' -');
    if fitted,
        fprintf(wfile,'%s\n','### Description of data set ###');
        fprintf(wfile,'%s\n%s\n','Source file:',handles.source_file);
        fprintf(wfile,'\n%s\n','### Pre-processing ###');
        fprintf(wfile,'%s%d%s\n','Time shift: ',handles.zero_time,' ns');
        fprintf(wfile,'%s\n','No phase correction by MMM. Only real part considered.');
        if ~ff_flag,
            fprintf(wfile,'Background correction: %s\n',bckg_mode);
            fprintf(wfile,'%s%6.3f\n','Background decay constant: ',handles.bckg_k);
        end;
        fprintf(wfile,'%s%5.3f\n','Fitted modulation depth: ',handles.mod_depth);
        fprintf(wfile,'%s%9.6f\n','r.m.s. error of fit: ',rms);
    end;
    if ~isempty(handles.ff_multi),
        fprintf(wfile,'Multi-spin effects were considered for total modulation depth of %5.3f.\n',handles.exp_depth);
    elseif length(handles.labels)>2,
        fprintf(wfile,'### WARNING ### More than 2 spins, but multi-spin effects were not considered.\n');
    end;
    fprintf(wfile,'Suggested minimum of the total dipolar evolution time: %3.1f microseconds\n',handles.tmax_min);
    fprintf(wfile,'Suggested optimum of the total dipolar evolution time: %3.1f microseconds\n',handles.tmax_opt);
    fclose(wfile);

    function myCleanup

        if ~isempty(wfile==fopen('all')),
            fclose(wfile);
        end;

    end
end