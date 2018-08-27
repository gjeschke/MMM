function [structure,chain,modtag,residue,atom,location,type] = dissect_address(address)
% function [structure,chain,modtag,residue,atom,location] = dissect_address(address)
%
% Dissects an MMM address into its parts, missing parts are output as empty
% strings
%
% G. Jeschke, 2012

structure='';
chain='';
modtag='';
residue='';
atom='';
location='';
type='';

stres=1;

st=strfind(address,'[');
en=strfind(address,']');
if ~isempty(st) && ~isempty(en),
    structure=address(st:en);
    stres=en+1;
end;

st=strfind(address,'(');
en=strfind(address,')');
if ~isempty(st) && ~isempty(en),
    chain=address(st:en);
    stres=en+1;
end;

st=strfind(address,'{');
en=strfind(address,'}');
if ~isempty(st) && ~isempty(en),
    modtag=address(st:en);
    stres=en+1;
end;

enres=strfind(address,'.');
enatom=strfind(address,':');
enloc=strfind(address,';');

if isempty(enres) && isempty(enatom) && isempty(enloc),
    if length(address)>=stres,
        residue=address(stres:end);
        return;
    end;
end;

if ~isempty(enres) && enres>stres,
    residue=address(stres:enres-1);
end;

if ~isempty(enres) && length(address)>enres,
    if isempty(enatom),
        if isempty(enloc),
            atom=address(enres+1:end);
            return
        elseif enloc-enres>1
            atom=address(enres+1:enloc-1);
        end;
    elseif enatom-enres>1
        atom=address(enres+1:enatom-1);
    end;
end;

if length(address)>enatom,
    if isempty(enloc),
        location=address(enatom+1:end);
        return;
    elseif enloc-enatom>1
        location=address(enatom+1:enloc-1);
    end;
end;

if length(address)>enloc,
    type=address(enloc+1:end);
end;