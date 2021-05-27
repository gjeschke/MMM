function adapt_rotamer_library(T,fname)

sezer=false;

cull=true;

cull_tags=':mm:tp:tm:';
allowed=3;
cull_not_found=0.01;

not_found=0.15; % population fraction of rotamer types that are assumed 
                % to be not (yet) found in crystal structures
xray_tags=':mm:tp:tm:pl:pp:pm:mt:';
xray_occ=[21,3,3,2,2,1,1];

if sezer,
    xray_tags=':mt:mm:tt:tp:tm:';
    xray_occ=[0.314,0.271,0.061,0.066,0.02];
    not_found=1-sum(xray_occ);
end;

nx=length(xray_occ);
xray_occ=(1-not_found)*xray_occ/sum(xray_occ);

fprintf(1,'Population fractions of observed rotamer groups:\n');
for k=1:nx,
    tag=id2tag(k,xray_tags);
    fprintf(1,'%s: %5.2f%%\n',tag,100*xray_occ(k));
end;

iname=sprintf('R1A_%iK_090619.mat',T);
load(iname);

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
if cull,
    nn=maxid-allowed;
    npop=cull_not_found/nn;
else
    nn=maxid-nx; % number of rotamer groups not observed in xray structures
    npop=not_found/nn; % assumed population fraction for unobserved rotamer groups
end;
fprintf(1,'Population fraction for other rotamer groups: %5.2f%%\n',100*npop);
% determine total populations of chi1-chi2 rotamer groups 
fprintf(1,'Rotamer group populations for library %s\n',iname);
pops=zeros(1,maxid);
spop_allowed=0;
spop_culled=0;
for k=1:maxid,
    sumpop=0;
    for kr=1:r,
        if chi12_id(kr)==k,
            sumpop=sumpop+rot_lib.calibration.pop(kr);
        end;
    end;
    pops(k)=sumpop;
    tag=id2tag(k,chi12_tags);
    id=tag2id(tag,cull_tags);
    if isempty(id),
        spop_culled=spop_culled+sumpop;
    else
        spop_allowed=spop_allowed+sumpop;
    end;
end;
spop_culled=spop_culled/(spop_culled+spop_allowed);
spop_allowed=1-spop_culled;
cull_fac=(1-cull_not_found)/spop_allowed;

frames=sum(pops);
fprintf(1,'Statistics for %i frames\n',frames);
spops=pops/frames;
for k=1:maxid,
    tag=id2tag(k,chi12_tags);
    fprintf(1,'%s: %5.2f%%\n',tag,100*spops(k));
end;

full_sum=0;
for kr=1:r,
    oldpop=spops(chi12_id(kr));
    tag=id2tag(chi12_id(kr),chi12_tags);
    if cull,
        id=tag2id(tag,cull_tags);
    else
        id=tag2id(tag,xray_tags);
    end;
    if isempty(id),
        newpop=npop;
    else
        if cull,
            newpop=cull_fac*oldpop;
        else
            newpop=xray_occ(id);
        end;
    end;
    rot_lib.calibration.pop(kr)=newpop*rot_lib.calibration.pop(kr)/oldpop;
    full_sum=full_sum+rot_lib.calibration.pop(kr);
end;
fprintf(1,'Population test sum is %12.2f\n',full_sum);

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

generic_canonicals=[0,60,120,180,-120,-60];
generic_can_char=':o:p:u:t:l:m:';

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
            fprintf(2,'Warning: Unsafe rotamer character for dihedral %i in rotamer %i of library %s\n',k,r,library);
            fprintf(2,'Difference is %4.1f degree between found angle %4.1f and canonical angle %4.1f\n',diff,ang,chi_table(id));
        end;
        my_character=id2tag(id,generic_can_char);
        char_table(r,k)=my_character;
    end;
end;