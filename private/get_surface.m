function [argout,msg]=get_surface(address,property)
% Get a property of a surface specified by its address
% the surface address must begin with $, otherwise an empty object array is
% returned
% an empty object array is also returned, if a surface with the suppliedd
% address does not exist
% surface addresses have the format $type:tag
% where type is the surface type (currently density or SAS) and tag the
% surface tag
%
% G. Jeschke, 2009

global model

argout=[];
msg.error=0;
msg.text='No error.';

if ~isfield(model,'surfaces') || isempty(model.surfaces), 
    msg.error=1;
    msg.text='ERROR: No surfaces defined in this model.';
    return; 
end; % return immediately if there are no surfaces

if ~strcmp(address(1),'$') || length(address)<2, % return immediately if address is invalid
    msg.error=2;
    msg.text='ERROR: Surface address must begin with $.';
    return
else
    mytag=address(2:end);
    found=0;
    for k=1:length(model.surfaces),
        if strcmp(mytag,model.surfaces(k).tag),
            found=1;
            id=k;
            break
        end;
    end;
end;

if ~found,
    msg.error=3;
    msg.text=sprintf('ERROR: Surface with address "%s" does not exist.',address);
    return
end;

switch property
    case 'all'
        argout=model.surfaces(id);
    case 'graphics'
        argout=model.surfaces(id).gobjects;
    case 'color'
        argout=model.surfaces(id).color;
    case 'level'
        argout=model.surfaces(id).level;
    case 'transparency'
        argout=model.surfaces(id).transparency;
    case 'type'
        [type,tag]=strtok(model.surfaces(id).tag,':');
        argout=type;
    otherwise
        msg.error=4;
        msg.text=sprintf('ERROR: Specified property "%s" does not exist.',property);
end;
