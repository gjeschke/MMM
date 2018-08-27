function char_table=rotamer_char_table(library)
% function char_table=rotamer_characters(library)
%   returns a rotamer character table for a given rotamer library
% 
% library       name of the rotamer library
%
% char_table    charater table, cell array of n strings, where n is the
%               number of rotamers in the library
%               characters are:
%               o   near 0 degree
%               p   near 60 degree
%               u   near 120 degree
%               t   near 180/-180 degree
%               l   near -120 degree
%               m   near -60 degree
%               list, e.g. 'm,m,t,t,p'
%
% returns empty cell array if library does not exist
%
% G. Jeschke, 2010

load(library);
if ~exist('rot_lib','var'),
    char_table={};
    return;
end;

generic_canonicals=[0,60,120,180,-120,-60];
generic_can_char=':o:p:u:t:l:m:';


% assign characters for all rotamers in library
rotamers=length(rot_lib.library);
char_table=cell(1,rotamers);
for r=1:rotamers,
    rotamer_def=rot_lib.library(r).dihedrals;
    char_string='';
    for k=1:length(rotamer_def),
        ang=rotamer_def(k);
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
            fprintf(2,'Warning: Unsafe rotamer character for dihedral %i in rotamer %i of library %s\n',k,r,library);
            fprintf(2,'Difference is %4.1f degree between found angle %4.1f and canonical angle %4.1f\n',diff,ang,chi_table(id));
        end;
        my_character=id2tag(id,generic_can_char);
        char_string=[char_string my_character ','];
    end;
    char_table{r}=char_string(1:end-1); % remove last comma
end;