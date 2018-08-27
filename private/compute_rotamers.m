function msg=compute_rotamers(indices,label,temperature,no_context)
% Compute rotamers for spin labeling using a rotamer library
%
% rotamers for the spin label with short name label are computed for all 
% residues specified by array indices using the rotamer library for the 
% temperature whose inverse is closest to the selected temperature
%
% returns a message record that reports on missing indices, missing rotamer
% library, or wrong temperature argument
%
% temperature   'cryo' (corresponds to 175 K)
%               'ambient' (298 K)
%               a string translating to a number (temperature in K)
% label         'IA1','IAP','IA-PROXYL','MTSL','MTSSL','R1A' (case-insensitive)
% no_context    optional flag that requests context-free computation,
%               defaults to false
%
% Results are stored in the field 
% model.sites{j} for the j-th scan performed on this model
% structure:
% model.sites{j}(k)             information for the k-th chain model
% model.sites{j}(k).residue(l)  rotamer analysis for l-th residue, with
% ...residue(l).label           label type
% ...residue(l).indices         indices of the labeled residue
% ...residue(l).NOpos           xyz coordinates and populations for all
%                               label positions
% ...residue(l).rotamers(n)     information for the n-th leading rotamer
%   ...rotamers(n).pop          populations of the leading rotamers making
%                               up 99.5% of the population, sorted by
%                               decreasing population
%   ...rotamers(n).coor         coordinates of the leading rotamers
%
% G. Jeschke, 2009

global hMain
global model
global rotamer_libraries

partition_function_threshold=0.05;

msg.error=0;
msg.text='No error.';

if ~exist('no_context','var')
    no_context = false;
end;

label=upper(label);
switch label
    case {'R1A','MTSSL'}
        label='MTSL';
    case {'IA1','IAP'}
        label='IA-PROXYL';
end;
    
[temperature,rem]=strtok(temperature); % get rid of unit, if provided

libflag=false;
switch temperature
    case 'cryo'
        T=175;
    case 'ambient'
        T=298;
    case 'library'
        T = 298;
        libflag=true;
    otherwise
        T=str2double(temperature);
end;

if isnan(T),
    msg.error=2;
    msg.text='Invalid argument for rotamer library temperature.';
    return
end;

if isempty(indices),
    msg.error=1;
    msg.text='Nothing to label.';
    return;
end;

if libflag,
    if ~exist(label,'file'),
        msg.error=13;
        msg.text='Library does not exist.';
        return;
    end;
    found=true;
    library=label;
    load(library);
    T = rot_lib.calibration.T;
else
    libs=length(rotamer_libraries);
    found=0;
    for k=1:libs,
        if strcmpi(rotamer_libraries(k).label,strtrim(label)) || strcmpi(rotamer_libraries(k).tc,strtrim(label)),
            Tvec=rotamer_libraries(k).T;
            found=1;
            Tinv=1/T;
            Tinv_vec=ones(size(Tvec))./Tvec;
            [mi,id]=min(abs(Tinv_vec-Tinv));
            library=id2tag(id,rotamer_libraries(k).files);
            break
        end;
    end;
end;

if ~found,
    msg.error=3;
    msg.text=sprintf('No rotamer library for label %s.',label);
    return;
end;

load(library);

scanned=0;

[mmm,nn]=size(indices);
cindices=zeros(mmm,4);
poi=0;
for k=1:mmm,
    aind=indices(k,:);
    aind=aind(aind>0);
    if length(aind)==4,
        poi=poi+1;
        cindices(poi,:)=aind;
    end;
end;
cindices=cindices(1:poi,:);

starting_time=clock;
if poi>=1,
    add_msg_board(sprintf('Scanning %i selected residues.',mmm));
    scanned=scanned+mmm;
    calc_positions=rotamer_populations(cindices,rot_lib,T,0,-1,'','',no_context);
    sites(k).residue=calc_positions;
    for kk=1:length(calc_positions),
        text=sprintf('rotamers computed: %s using library %s at a temperature of %4.0f K',calc_positions(kk).label,library,T);
        add_annotation(calc_positions(kk).indices,'Spin',text,{'rotamers computed'});
        numr=length(calc_positions(kk).rotamers);
        rmsd=NOpos_rmsd(calc_positions(kk).NOpos);
        rmsdz=NOpos_rmsdz(calc_positions(kk).NOpos);
        Z=calc_positions(kk).partition_function;
        add_annotation(calc_positions(kk).indices,'Spin',sprintf('%i significant rotamers with partition function %7.5f',numr,Z),{'rotamers computed'});
        if Z<partition_function_threshold,
            add_annotation(calc_positions(kk).indices,'Spin','Tight site. Labeling may fail here.',{'rotamers computed'});
        end;
        text=sprintf('NO position r.m.s.d. %4.2f nm',rmsd);
        if hMain.z_analysis,
            text=sprintf('%s; NO z position r.m.s.d. %4.2f nm',text,rmsdz);
        end;
        add_annotation(calc_positions(kk).indices,'Spin',text,{'rotamers computed'});                
    end;
else
    add_msg_board('Warning: No residues were selected for rotamer computation.');
end;

if scanned>0,
    scan=1;
    if isfield(model,'sites'),
        scan=length(model.sites)+1;
    end;
    model.sites{scan}=sites;
    needed_time=etime(clock,starting_time);
    hour=floor(needed_time/3600);
    minute=floor((needed_time-3600*hour)/60);
    sec=floor(needed_time-3600*hour-60*minute);
    add_msg_board(sprintf('Site scan completed in %i h, %i min, %i s',hour,minute,sec));
end;


function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function rmsd=NOpos_rmsdz(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
zmean=sum(NOall(:,3).*pop);
dz=(NOall(:,3)-zmean);
nNO=length(dz);
rmsd=sqrt(0.005/3+nNO*sum(dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm
