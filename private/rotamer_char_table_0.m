function char_table=rotamer_char_table(library)
% function char_table=rotamer_characters(library)
%   returns a rotamer character table for a given rotamer library
% 
% library       name of the rotamer library
%
% char_table    charater table, cell array of n strings, where n is the
%               number of rotamers in the library
%               characters for each individual dihedral angle:
%               2 canonical values:     g0  near 0 degree
%                                       o+  near 90 degree
%                                       g   near 180 degree
%                                       o-  near -90 degree
%               3 canonical values:     g0  near 0 degree
%                                       t+  near 60 degree
%                                       g   near 180 degree
%                                       t-  near -60 degree
%               4 canonical values:     g0  near 0 degree
%                                       t+  near 60 degree
%                                       o+  near 120 degree
%                                       g   near 180 degree
%                                       o-  near -120 degree
%                                       t-  near -60 degree
%               these characters are derived from rot_lib.chi_select
%               individual rotamers are assigned the characters
%               corresponding to the values in rot_lib.chi_select that are
%               closest to the actual value
%               the character string for each rotamer is a comma-separated
%               list, e.g. 't+,t-,o+,t+,t-'
%
% returns empty cell array if library does not exist
%
% G. Jeschke, 2010

load(library);
if ~exist('rot_lib','var'),
    char_table={};
    return;
end;

generic_canonical=[0,60,120,180,-180,-120,-60];
generic_can_char=':g0:t+:o+:g:g:o-:t-:';

% make table of canonical characters
dihedrals=length(rot_lib.chi_select);
canonical_val=zeros(1,dihedrals);
characters=cell(dihedrals,4);
for k=1:dihedrals,
    angles=rot_lib.chi_select(k).angles;
    canonical_val(k)=length(angles);
    switch length(angles) % different characters for 2, 3 or 4 canonical values
        case 2
            can_char=':g0:o+:g:g:o-:';
            canonical=[0,90,180,-180,-90];
        case 3
            can_char=':g0:t+:g:g:t-:';
            canonical=[0,60,180,-180,-60];
        case 4
            can_char=':g0:t+:o+:g:g:o-:t-:';
            canonical=[0,60,120,180,-180,-120,-60];
    end
    can_char=generic_can_char;
    canonical=generic_canonical;
    for kk=1:length(angles),
        ang=angles(kk);
        while ang>180, ang=ang-360; end;
        while ang<-180, ang=ang+360; end;
        [diff,id]=min(abs(canonical-ang));
        if diff>30,
            fprintf(2,'Warning: Unsafe rotamer character for dihedral %i in library %s\n',kk,library);
            fprintf(2,'Difference is %4.1f degree between found angle %4.1f and canonical angle %4.1f\n',diff,ang,canonical(id));
        end;
        characters{k,kk}=id2tag(id,can_char);
    end;
end;

% assign characters for all rotamers in library
rotamers=length(rot_lib.library);
char_table=cell(1,rotamers);
for r=1:rotamers,
    rotamer_def=rot_lib.library(r).dihedrals;
    char_string='';
    for k=1:dihedrals,
        ang=rotamer_def(k);
        chi_table=rot_lib.chi_select(k).angles;
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
        my_character=characters{k,id};
        char_string=[char_string my_character ','];
    end;
    char_table{r}=char_string(1:end-1); % remove last comma
end;