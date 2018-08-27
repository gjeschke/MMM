function [coor,msg] = get_CYANA_atom_coor(adr,residue)
% function [coor,msg] = get_CYANA_atom_coor(adr,residue)
%
% Returns an atom coordinate conforming to CYANA restraint format, given a
% pseudo-MMM address, the pseudo-MMM address can be eithe a well-formed MMM
% atom address or can refer to a pseudo-atom as specified in the CYANA
% residue library file
%
% if the address does not contain a dot, an empty coordinate is returned
%
% adr       MMM (pseudo-)atom address
% residue   residue type, three-letter code
%
% coor      coordinate of the atom or pseudoatom, if there are several
%           locations, the mean coordinate is returned
% msg       error message in the format msg.error, msg.text 
%
% G. Jeschke, 17.8.2016

coor = [];
msg.error = 0;
msg.text = '';

point = strfind(adr,'.');
if isempty(point)
    msg.error = 1;
    msg.text = 'Invalid atom address';
    return
else
    restag = adr(1:point);
    atag = adr(point+1:end);
end

if atag(1) == 'Q' % this is a pseudo-atom specification
    switch upper(residue)
        case 'ALA'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'ARG'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QG'
                    [msg,coor1] = get_object(strcat(restag,'HG2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD'
                    [msg,coor1] = get_object(strcat(restag,'HD2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QH1'
                    [msg,coor1] = get_object(strcat(restag,'HH11'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HH12'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QH2'
                    [msg,coor1] = get_object(strcat(restag,'HH21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HH22'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'ASN'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD2'
                    [msg,coor1] = get_object(strcat(restag,'HD21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD22'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'ASP','CYS','CYSS','HIS','SER','TRP'}
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'GLN'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QG'
                    [msg,coor1] = get_object(strcat(restag,'HG2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QE1'
                    [msg,coor1] = get_object(strcat(restag,'HE21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HE22'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'GLU'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QG'
                    [msg,coor1] = get_object(strcat(restag,'HG2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'GLY'
            switch atag
                case 'QA'
                    [msg,coor1] = get_object(strcat(restag,'HA2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HA3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'ILE'
            switch atag
                case 'QG2'
                    [msg,coor1] = get_object(strcat(restag,'HG21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG22'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HG23'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                case 'QG1'
                    [msg,coor1] = get_object(strcat(restag,'HG12'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG13'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD1'
                    [msg,coor1] = get_object(strcat(restag,'HD11'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD12'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HD13'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'LEU'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD1'
                    [msg,coor1] = get_object(strcat(restag,'HD11'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD12'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HD13'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                case 'QD2'
                    [msg,coor1] = get_object(strcat(restag,'HD21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD22'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HD23'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                 case 'QQD'
                    [msg,coor1] = get_object(strcat(restag,'HD11'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD12'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HD13'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor4] = get_object(strcat(restag,'HD21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor5] = get_object(strcat(restag,'HD22'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor6] = get_object(strcat(restag,'HD23'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3+coor4+coor5+coor6)/6;
               otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'LYS'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QG'
                    [msg,coor1] = get_object(strcat(restag,'HG2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD'
                    [msg,coor1] = get_object(strcat(restag,'HD2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QE'
                    [msg,coor1] = get_object(strcat(restag,'HE2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HE3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QZ'
                    [msg,coor1] = get_object(strcat(restag,'HZ1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HZ2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HZ3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end
        case 'MET'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QG'
                    [msg,coor1] = get_object(strcat(restag,'HG2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QE'
                    [msg,coor1] = get_object(strcat(restag,'HE1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HE2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HE3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'PHE','TYR'}
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD'
                    [msg,coor1] = get_object(strcat(restag,'HD1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD2'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QE'
                    [msg,coor1] = get_object(strcat(restag,'HE1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HE2'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QR'
                    [msg,coor1] = get_object(strcat(restag,'HD1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HE1'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor4] = get_object(strcat(restag,'HE2'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3+coor4)/4;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'PRO'
            switch atag
                case 'QB'
                    [msg,coor1] = get_object(strcat(restag,'HB2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HB3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QG'
                    [msg,coor1] = get_object(strcat(restag,'HG2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'QD'
                    [msg,coor1] = get_object(strcat(restag,'HD2'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HD3'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'THR'
            switch atag
                case 'QG2'
                    [msg,coor1] = get_object(strcat(restag,'HG21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG22'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HG23'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case 'VAL'
            switch atag
                case 'QG1'
                    [msg,coor1] = get_object(strcat(restag,'HG11'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG12'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HG13'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/2;
                case 'QG2'
                    [msg,coor1] = get_object(strcat(restag,'HG21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG22'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HG23'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                case 'QQG'
                    [msg,coor1] = get_object(strcat(restag,'HG11'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'HG12'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'HG13'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor4] = get_object(strcat(restag,'HG21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor5] = get_object(strcat(restag,'HG22'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor6] = get_object(strcat(restag,'HG23'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3+coor4+coor5+coor6)/6;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'ADE',' DA'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q2'''
                    [msg,coor1] = get_object(strcat(restag,'H2'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H2'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q6'
                    [msg,coor1] = get_object(strcat(restag,'H61'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H62'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'RADE','  A'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q6'
                    [msg,coor1] = get_object(strcat(restag,'H61'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H62'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'CYT',' DC'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q2'''
                    [msg,coor1] = get_object(strcat(restag,'H2'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H2'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q4'
                    [msg,coor1] = get_object(strcat(restag,'H41'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H42'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'RCYT','  C'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q4'
                    [msg,coor1] = get_object(strcat(restag,'H41'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H42'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'GUA',' DG'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q2'''
                    [msg,coor1] = get_object(strcat(restag,'H2'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H2'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q2'
                    [msg,coor1] = get_object(strcat(restag,'H21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H22'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'RGUA','  G'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q2'
                    [msg,coor1] = get_object(strcat(restag,'H21'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H22'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'THY',' DT'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q2'''
                    [msg,coor1] = get_object(strcat(restag,'H2'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H2'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                case 'Q7'
                    [msg,coor1] = get_object(strcat(restag,'H71'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H72'),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor3] = get_object(strcat(restag,'H73'),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2+coor3)/3;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        case {'URA','  U'}
            switch atag
                case 'Q5'''
                    [msg,coor1] = get_object(strcat(restag,'H5'''),'coor');
                    if msg.error
                        return
                    end
                    [msg,coor2] = get_object(strcat(restag,'H5'''''),'coor');
                    if msg.error
                        return
                    end
                    coor = (coor1+coor2)/2;
                otherwise
                    msg.error = 3;
                    msg.text = sprintf('Pseudo-atom %s in residue type %s is unknown.',atag,residue);
                    return
            end                    
        otherwise
            msg.error = 2;
            msg.text = sprintf('Pseudo-atom in unknown residue type %s.',residue);
            return
    end
else
    [msg,coor] = get_object(adr,'coor');
end;


