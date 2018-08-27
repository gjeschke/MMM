function character=rotamer_character(chivec)
% function char_table=rotamer_character(chi1,chi2,chi3,chi4,chi5)
%   returns a rotamer character string for given dihedrals
% 
% chivec        vector of dihedral angles (degree)
%
% character     character string
%               characters are:
%               o   near 0 degree
%               p   near 60 degree
%               u   near 120 degree
%               t   near 180/-180 degree
%               l   near -120 degree
%               m   near -60 degree
%               list, e.g. 'm,m,t,t,p'
%
% G. Jeschke, 2012

generic_canonicals=[0,60,120,180,-120,-60];
generic_can_char=':o:p:u:t:l:m:';


char_string='';
for k=1:length(chivec),
    ang=chivec(k);
    chi_table=generic_canonicals;
    for kk=1:length(chi_table),
        while chi_table(kk)>180, chi_table(kk)=chi_table(kk)-360; end;
        while chi_table(kk)<-180, chi_table(kk)=chi_table(kk)+360; end; 
    end;
    while ang>180, ang=ang-360; end;
    while ang<-180, ang=ang+360; end; 
    [diff1,id1]=min(abs(chi_table-ang));
    [diff2,id2]=min(abs(chi_table-ang+360));
    [diff3,id3]=min(abs(chi_table-ang-360));
    ids=[id1 id2 id3];
    diffs=[diff1 diff2 diff3];
    [diff,dpoi]=min(abs(diffs));
    id=ids(dpoi);
    if diff>30,
        fprintf(2,'Warning: Unsafe rotamer character for dihedral %i\n',k);
        fprintf(2,'Difference is %4.1f degree between found angle %4.1f and canonical angle %4.1f\n',diff,ang,chi_table(id));
    end;
    my_character=id2tag(id,generic_can_char);
    char_string=[char_string my_character ','];
end;
character=char_string(1:end-1); % remove last comma