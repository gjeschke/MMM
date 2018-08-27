function msg = replace_residues(chimera,template,alignment)
% function replace_residues(chimera,template,alignment)
%
% Replaces a range of residues in chimera by the corresponding residues of
% template, correspondence of residue numbers is specified by table
% alignment
% only coordinates are replaced atom by atom
% if there is no atom-by-atom correspondence, an error is reported
% corresponding atoms must have the same number of locations, otherwise an
% error is reported
%
% the chain chimera should not have any displayed graphics elements,
% otherwise graphics is inconsistent after replacement
%
% chimera   indices of a chain
% template  indices of a chain
% alignment array [res,2] of aligned residues in the chain, they must be of
%           the same residue type
%
% msg       message structure reporting on success
%           .error  error number, 0 for no error
%           .text   error text
%
% G. Jeschke, 14.9.2017

global model

msg.error = 0;
msg.text = '';

[res,~] = size(alignment);

chimera_adr = mk_address(chimera);
template_adr = mk_address(template);

for kr = 1:res
    [indc,msg0] =resolve_address(sprintf('%s%i',chimera_adr,alignment(kr,1)));
    if msg0.error
        msg.error = 1;
        msg.text = sprintf('Residue %i in chimera could not be resolved (%s)',alignment(kr,1),msg0.text);
        return
    end
    [msg0,catoms] = get_residue(indc,'children');
    if msg0.error
        msg.error = 2;
        msg.text = sprintf('For residue %i in chimera atoms could not be retrieved (%s)',alignment(kr,2),msg0.text);
        return
    end
    [nat,~] = size(catoms); % all atom indices for atoms to be replaced
    for ka = 1:nat
        [~,clocations] = get_atom(catoms(ka,:),'children');
        [nlc,~] = size(clocations); 
        [stag,~,~,~,~,atag] = mk_address_parts(catoms(ka,:));
        if isempty(stag)
            msg.error = 3;
            msg.text = sprintf('Atom resolution error in chimera');
            return
        end
        [tempat,msg0] = resolve_address(sprintf('%s%i.%s',template_adr,alignment(kr,2),atag));
        if msg0.error
            msg.error = 4;
            msg.text = sprintf('Atom %s in residue %i in template could not be resolved (%s)',atag,alignment(kr,2),msg0.text);
            return
        end
        [~,tlocations] = get_atom(tempat,'children');
        [nlt,~] = size(tlocations);
        if nlt ~= nlc
            msg.error = 5;
            msg.text = sprintf('Atom %s in residue %i in template has mismatch in number of locations',atag,alignment(kr,2));
            return
        end
        for kl = 1:nlc
            [msg0,xyz] = get_location(tlocations(kl,:),'xyz');
            if msg0.error
                msg.error = 6;
                msg.text = sprintf('Failed to retrieve coordinates for atom %s in residue %i in template (%s)',atag,alignment(kr,2),msg0.text);
                return
            end
            passit{1} = xyz;
            msg0 = set_location(clocations(kl,:),'xyz',passit);
            if msg0.error
                msg.error = 7;
                msg.text = sprintf('Failed to set coordinates for atom %s in residue %i in chimera (%s)',atag,alignment(kr,1),msg0.text);
                return
            end
        end
    end
end

