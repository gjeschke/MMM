function analyze_rotamer_library(fname)

load(fname);

char_table=rotamer_char_strings(rot_lib);

% make a list of different chi1-chi2 rotamers
chi12_tags=':';
[r,nchi]=size(char_table);
chi12_id=zeros(1,r);
maxid=0;
for k=1:r,
    tag=char_table(k,1:2);
    id=tag2id(tag,chi12_tags);
    if isempty(id),
        chi12_tags=[chi12_tags tag ':'];
        id=tag2id(tag,chi12_tags);
        maxid=id;
    end;
    chi12_id(k)=id;
end;
% determine total populations of chi1-chi2 rotamer groups 
fprintf(1,'Rotamer group populations for library %s\n',fname);
pops=zeros(1,maxid);
for k=1:maxid,
    sumpop=0;
    for kr=1:r,
        if chi12_id(kr)==k,
            sumpop=sumpop+rot_lib.calibration.pop(kr);
        end;
    end;
    pops(k)=sumpop;
end;
frames=sum(pops);
fprintf(1,'Statistics for %i frames\n',frames);
spops=pops/frames;
for k=1:maxid,
    tag=id2tag(k,chi12_tags);
    fprintf(1,'%s: %5.2f%%\n',tag,100*spops(k));
end;

% Analysis of transition rates
transitions=rot_lib.calibration.transitions;
[m,~]=size(transitions);
chi_trans=zeros(1,nchi);
trans13=0;
trans45=0;
for k=1:m,
    char1=char_table(k,:);
    for kk=1:m,
        char2=char_table(kk,:);
        for kc=1:nchi,
            if char1(kc)~=char2(kc),
                chi_trans(kc)=chi_trans(kc)+transitions(k,kk);
            end;
        end;
        if ~strcmp(char1(1:3),char2(1:3)),
            trans13=trans13+transitions(k,kk);
        end;
        if ~strcmp(char1(4:5),char2(4:5)),
            trans45=trans45++transitions(k,kk);
        end;
    end;
end;

fprintf(1,'Rotamer transition analysis:\n');
for kc=1:nchi,
    fprintf(1,'chi%i: %5.2f transitions/ns\n',kc,chi_trans(kc)/100);
end;
fprintf(1,'chi1-3: %5.2f transitions/ns\n',trans13/100);
fprintf(1,'chi4-5: %5.2f transitions/ns\n',trans45/100);

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

generic_canonicals=[0,60,120,180,-120,-60];
generic_can_char=':o:p:u:t:l:m:';

% generic_canonicals=[60,180,-60];
% generic_can_char=':p:t:m:';

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
        if diff>30,
            fprintf(2,'Warning: Unsafe rotamer character for dihedral %i in rotamer %i\n',k,r);
            fprintf(2,'Difference is %4.1f degree between found angle %4.1f and canonical angle %4.1f\n',diff,ang,chi_table(id));
        end;
        my_character=id2tag(id,generic_can_char);
        char_table(r,k)=my_character;
    end;
end;