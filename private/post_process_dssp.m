function post_process_dssp(dssp,snum)
% function dssp=post_process_dssp(dssp,snum)
%
% Cleanes up stretches of residues to get a secondary structure assignment
% that is better adapted to the viewing habits of ribbon model users and
% stores this assignment in the MMM structure model
%
% dssp  a Matlab structure variable with secondary structure information as
%       obtained with rd_dssp.m
% snum  structure index
%
% G. Jeschke, 2009

global model
global third_party

dssp_ref=true;
id=tag2id('Kabsch:1983_dssp',third_party.tags,[],'|');
if isfield(model,'auto_references'),
    if ~isempty(find(id==model.auto_references, 1)),
        dssp_ref=false;
    end;
else
    model.auto_references=[];
end;
if dssp_ref,
    if ~isfield(model,'references'),
        model.references(1)=third_party.references(id);
    elseif isempty(model.references)
        model=rmfield(model,'references');
        model.references(1)=third_party.references(id);
    else
        model.references(end+1)=third_party.references(id);
    end;
    model.auto_references(end+1)=id;
end;
    
extend_threshold=3; % threshold for extending a secondary structure element 
                    % (in multiples of the standard deviation of the phi
                    % and psi backbone dihedrals in the existing stretch
                    % larger numbers mean more generous granting of helix
                    % or sheet status
residues=length(dssp);
indices=zeros(residues,4);

% make a table of continuous chain segments and store dssp secondary
% structure identifier for residues
segments=zeros(200,2);
cpoi=1;
segments(1,1)=1;
cnum=tag2id(dssp(1).chain,model.chain_tags{snum});
for k=1:residues,
    if strcmp(dssp(k).slc,'!'), % the DSSP chain break identifier
        cpoi=cpoi+1;
        segments(cpoi,1)=k+1;
        segments(cpoi-1,2)=k-1;
        cnum=tag2id(dssp(k+1).chain,model.chain_tags{snum});
    else % this appears to be a residue
        if isempty(cnum);
            cnum=tag2id(dssp(k).chain,model.chain_tags{snum});
        end;
        if ~isempty(cnum),
            rnum=tag2id(dssp(k).tag,model.structures{snum}(cnum).residues{1}.residue_tags); % does not yet treat several chain models
            if ~isempty(rnum),
                model.structures{snum}(cnum).residues{1}.info(rnum).dssp=dssp(k).sec;
                indices(k,:)=[snum cnum 1 rnum];
            end;
        end;
    end;
end;
if ~strcmp(dssp(k).slc,'!'), % final chain break identifier might be missing
    segments(cpoi,2)=k;
end;
segments=segments(1:cpoi,:);

% try to straighten out assignments by using context
clean_sec=zeros(1,residues);
for k=1:cpoi,
    %scan for continuous stretches of residues with same simplified
    %secondary structure
    stretches=zeros(50,2);
    sec=get_sec(segments(k,1),dssp);
    clean_sec(segments(k,1))=sec;
    stretches(1,1)=segments(k,1);
    poi=1;
    for kk=segments(k,1):segments(k,2),
        nsec=get_sec(kk,dssp);
        clean_sec(kk)=nsec;
        if sec~=nsec,
            stretches(poi,2)=kk-1;
            poi=poi+1;
            stretches(poi,1)=kk;
            sec=nsec;
        end;
    end;
    stretches(poi,2)=segments(k,2);
    % extend these stretches, based on similar backbone dihedrals
    extensions=0;
    for kk=1:poi,
        l_overlap=0;
        h_overlap=0;
        if clean_sec(stretches(kk,1))>0
            extend=1;
        else
            extend=0;
        end;
        while ~(l_overlap && h_overlap) && extend,
            extend=0;
            if stretches(kk,1)-1<1,
                l_overlap=1;
            end;
            if kk>1 && stretches(kk-1,2)>=stretches(kk,1),
                l_overlap=1;
            end;
            if stretches(kk,2)+1>residues,
                h_overlap=1;
            end;
            if kk<poi && stretches(kk+1,1)<=stretches(kk,2),
                h_overlap=1;
            end;
            if ~(l_overlap && h_overlap),
                phi=zeros(1,stretches(kk,2)-stretches(kk,1)+1);
                psi=zeros(1,stretches(kk,2)-stretches(kk,1)+1);
                kkk=0;
                for ks=stretches(kk,1):stretches(kk,2),
                    kkk=kkk+1;
                    phi(kkk)=dssp(ks).phi;
                    psi(kkk)=dssp(ks).psi;
                end;
                mphi=mean(phi);
                sphi=std(phi);
                mpsi=mean(phi);
                spsi=std(phi);
                if ~l_overlap,
                    if abs(dssp(stretches(kk,1)-1).phi-mphi)<extend_threshold*sphi,
                        if abs(dssp(stretches(kk,1)-1).psi-mpsi)<extend_threshold*spsi,
                            extend=1;
                            extensions=extensions+1;
                            stretches(kk,1)=stretches(kk,1)-1;
                            clean_sec(stretches(kk,1)-1)=clean_sec(stretches(kk,1));
                        end;
                    end;
                end;
                if ~h_overlap,
                    if abs(dssp(stretches(kk,2)+1).phi-mphi)<extend_threshold*sphi,
                        if abs(dssp(stretches(kk,2)+1).psi-mpsi)<extend_threshold*spsi,
                            extend=1;
                            extensions=extensions+1;
                            stretches(kk,2)=stretches(kk,2)+1;
                            clean_sec(stretches(kk,2)+1)=clean_sec(stretches(kk,2));
                        end;
                    end;
                end;
            end;
        end;
    end;
%     disp(sprintf('%i secondary structure element extensions performed.',extensions));
    % reassign secondary structure
    for kk=1:residues,
        if sum(indices(kk,:))~=0,
            model.structures{indices(kk,1)}(indices(kk,2)).residues{indices(kk,3)}.info(indices(kk,4)).secondary=clean_sec(kk);
        end;
    end;
end;

function sec=get_sec(poi,dssp)
% simplified secondary structure code
switch dssp(poi).sec
    case 'H'
        sec=1;
    case 'E';
        sec=2;
    case 'G';
        sec=1;
    case 'I';
        sec=1;
    otherwise
        sec=0;
end;