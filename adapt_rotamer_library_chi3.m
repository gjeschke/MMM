function adapt_rotamer_library_chi3(T,chi3_mode,fname)


iname=sprintf('R1A_%iK_090619.mat',T);
load(iname);

char_table=rotamer_char_strings(rot_lib);

[r,~]=size(char_table);
sum0=sum(rot_lib.calibration.pop);
sumpop=0;
cstr='';
for k=1:r,
    cstr=[cstr char_table(k,3)];
    if char_table(k,3) ~= chi3_mode,
        rot_lib.calibration.pop(k)=0;
    end;
    sumpop=sumpop+rot_lib.calibration.pop(k);
end;
fprintf(1,'Population fraction of %s rotamers: %5.2f%%\n',chi3_mode,100*sumpop/sum0);
rot_lib.calibration.pop=sum0*rot_lib.calibration.pop/sumpop;
fprintf(1,'Population test sum is %12.2f\n',sum(rot_lib.calibration.pop));
disp(cstr);
save(fname,'rot_lib');

function char_table=rotamer_char_strings(rot_lib)
% function char_table=rotamer_char_strings(rot_lib)
%   returns a rotamer character table for a given rotamer library
% 
% rot_lib       rotamer library
%
% char_table    charater table, array of n strings, where n is the
%               number of rotamers in the library
%               characters are:
%               o   near 0 degree
%               p   near 60 degree
%               u   near 120 degree
%               t   near 180/-180 degree
%               l   near -120 degree
%               m   near -60 degree
%               string, e.g. 'mmttp'
%
% returns empty cell array if library does not exist
%
% G. Jeschke, 2010

generic_canonicals=[90,-90];
generic_can_char=':p:m:';

rotamer_def=rot_lib.library(1).dihedrals;
rotamers=length(rot_lib.library);
nchi=length(rotamer_def);
char_table=char(rotamers,nchi);

% assign characters for all rotamers in library
for r=1:rotamers,
    rotamer_def=rot_lib.library(r).dihedrals;
    for k=1:nchi,
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
        if diff>75,
            fprintf(2,'Warning: Unsafe rotamer character for dihedral %i in rotamer %i\n',k,r);
            fprintf(2,'Difference is %4.1f degree between found angle %4.1f and canonical angle %4.1f\n',diff,ang,chi_table(id));
        end;
        my_character=id2tag(id,generic_can_char);
        char_table(r,k)=my_character;
    end;
end;