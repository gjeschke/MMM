function msg=attach_spin_label(indices,label,temperature)
% function msg=attach_spin_label(indices,label,temperature)
%
% attach spin label with short name label to all residues specified by
% array indices using the rotamer library for the temperature whose
% inverse is closest to the selected temperature
%
% returns a message record that reports on missing indices, missing rotamer
% library, or wrong temperature argument
%
% temperature   'cryo' (corresponds to 175 K)
%               'ambient' (298 K)
%               a string translating to a number (temperature in K)
% label         'IA1','IAP','IA-PROXYL','MTSL','MTSSL','R1A' (case-insensitive)

global rotamer_libraries

msg.error=0;
msg.text='No error.';

label=upper(label);
switch label
    case {'R1A','MTSSL'}
        label='MTSL';
    case {'IA1','IAP'}
        label='IA-PROXYL';
end;
    
[temperature,rem]=strtok(temperature); % get rid of unit, if provided

switch temperature
    case 'cryo'
        T=175;
    case 'ambient'
        T=298;
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


libs=length(rotamer_libraries);
found=0;
for k=1:libs,
    if strcmpi(rotamer_libraries(k).label,strtrim(label)),
        Tvec=rotamer_libraries(k).T;
        found=1;
        Tinv=1/T;
        Tinv_vec=ones(size(Tvec))./Tvec;
        [mi,id]=min(abs(Tinv_vec-Tinv));
        library=id2tag(id,rotamer_libraries(k).files);
        break
    end;
end;

if ~found,
    msg.error=3;
    msg.text=sprintf('No rotamer library for label %s.',label);
    return;
end;

[m,n]=size(indices);

load(library);

calc_positions=rotamer_populations(indices,rot_lib,T,true);

set(gcf,'Pointer','watch');

n=length(calc_positions);
if n>=1,
    for k=1:n,
        set_residue(calc_positions(k).indices,'label',{calc_positions(k).label,calc_positions(k).rotamers,calc_positions(k).attach_frame,calc_positions(k).exclude});
        text=sprintf('label attached: %s using library %s at a temperature of %4.0f K',calc_positions(k).label,library,T);
        add_annotation(calc_positions(k).indices,'Spin',text,{'spin label attached'});
    end;
end;

set(gcf,'Pointer','arrow');
