function adapt_rotamer_library_chi1_chi3(T,fname)

R=1.985877534e-3; % gas constant in kcal/(mol K)

sezer=true;

%w_tags=':tpm:mmp:tmp:mmm:tmm:tpp:mtp:ttp:ttm:'; % Warshaviak et al., JPC B, 115, 397-405 (2011) 
w_tags=':tpm:mmp:tmu:mmm:tmm:tpp:mtp:ttp:ttm:'; % tmp replaced by tmu
w_E=[0,1.11,1.77,1.91,1.96,2.03,2.47,3.03,3.38]; % Table 3, total energies

nw=length(w_E);
w_assigned=zeros(1,nw);
w_pop=exp(-w_E/(R*T));

not_found=0.05; % total population fraction of rotamer types not treated 
                % by warshaviak et al.
if sezer,
    w_tags=':mtm:mtp:mmm:mmp:ttm:tpm:tmm:'; % tmp replaced by tmu
    w_pop=[0.236,0.078,0.145,0.126,0.061,0.066,0.02]; % Table 10
    nw=length(w_pop);
    not_found=1-sum(w_pop);
end;

w_pop=w_pop/sum(w_pop);

if sezer
    fprintf(1,'Sezer population distribution:\n');
else
    fprintf(1,'Warshaviak population distribution:\n');
end;
for k=1:nw,
    tag=id2tag(k,w_tags);
    fprintf(1,'%s: %6.2f%%\n',tag,100*w_pop(k));
end;


w_pop=(1-not_found)*w_pop/sum(w_pop);

iname=sprintf('R1A_%iK_090619.mat',T);
load(iname);

char_table=rotamer_char_strings(rot_lib);

% make a list of different chi1-chi3 rotamers
chi13_tags=':';
[r,nchi]=size(char_table);
chi13_id=zeros(1,r);
maxid=0;
for k=1:r,
    tag=char_table(k,1:3);
    id=tag2id(tag,chi13_tags);
    if isempty(id),
        chi13_tags=[chi13_tags tag ':'];
        id=tag2id(tag,chi13_tags);
        maxid=id;
    end;
    chi13_id(k)=id;
end;
nn=maxid-nw; % number of rotamer groups not treated by Warshaviak et al.
npop=not_found/nn; % assumed population fraction for unobserved rotamer groups
fprintf(1,'Population fraction for other rotamer groups: %5.2f%%\n',100*npop);
% determine total populations of chi1-chi3 rotamer groups 
fprintf(1,'Rotamer group populations for library %s\n',iname);
pops=zeros(1,maxid);
for k=1:maxid,
    sumpop=0;
    for kr=1:r,
        if chi13_id(kr)==k,
            sumpop=sumpop+rot_lib.calibration.pop(kr);
        end;
    end;
    pops(k)=sumpop;
end;

frames=sum(pops);
fprintf(1,'Statistics for %i frames\n',frames);
spops=pops/frames;
for k=1:maxid,
    tag=id2tag(k,chi13_tags);
    fprintf(1,'%s: %5.2f%%\n',tag,100*spops(k));
end;

full_sum=0;
for kr=1:r,
    oldpop=spops(chi13_id(kr));
    tag=id2tag(chi13_id(kr),chi13_tags);
    id=tag2id(tag,w_tags);
    if isempty(id),
        newpop=npop;
    else
        newpop=w_pop(id);
        w_assigned(id)=1;
    end;
    rot_lib.calibration.pop(kr)=newpop*rot_lib.calibration.pop(kr)/oldpop;
    full_sum=full_sum+rot_lib.calibration.pop(kr);
end;
fprintf(1,'Population test sum is %12.2f\n',full_sum);
for k=1:nw,
    if ~w_assigned(k),
        tag=id2tag(k,w_tags);
        fprintf(2,'Rotamer group %s does not occur in library.\n',tag);
    end;
end;

save(fname,'rot_lib');

fprintf(1,'Populations in adapted library\n');
% determine total populations of chi1-chi3 rotamer groups 
fprintf(1,'Rotamer group populations for library %s\n',fname);
pops=zeros(1,maxid);
for k=1:maxid,
    sumpop=0;
    for kr=1:r,
        if chi13_id(kr)==k,
            sumpop=sumpop+rot_lib.calibration.pop(kr);
        end;
    end;
    pops(k)=sumpop;
end;
frames=sum(pops);
fprintf(1,'Statistics for %i frames\n',frames);
spops=pops/frames;
for k=1:maxid,
    tag=id2tag(k,chi13_tags);
    fprintf(1,'%s: %5.2f%%\n',tag,100*spops(k));
end;


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