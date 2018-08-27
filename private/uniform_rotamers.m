function [rax,distr]=uniform_rotamers(adr1,adr2,handles)
% Creates a distance distribution for all rotamer pairs (approximately) uniformly populated,
% this is the broadest possible distance distribution
sig=0.1;

NOpos1=get_NOpos(adr1,handles);
NOpos2=get_NOpos(adr2,handles);

pop1=ones(size(NOpos1(:,4)));
pop1=pop1/length(pop1);
pop2=ones(size(NOpos2(:,4)));
pop2=pop2/length(pop2);
NOpos1(:,4)=pop1;
NOpos2(:,4)=pop2;

[rax,distr]=get_distribution(NOpos1,NOpos2,sig);

% tdistr=interp1(trax0,tdistr0,rax,'pchip',0);
% % sc=sum(tdistr.*distr)/sum(tdistr.*tdistr);
% % tdistr=sc*tdistr;
% tdistr=tdistr*sum(distr)/sum(tdistr);
% trax=rax;
        
function NOpos=get_NOpos(adr,handles)

global model

if strcmp(adr(1),'#');
    nt=str2double(adr(6:end));
    NOpos=handles.frames{nt};
    return;
end;

NOpos=[];

indices=resolve_address(adr);
if isempty(indices) || length(indices)~=4,
    return;
end;

nscan=length(model.sites);
for k=1:nscan,
    for kc=1:length(model.sites{k}),
        for kk=1:length(model.sites{k}(kc).residue),
            diff=sum(abs(indices-model.sites{k}(kc).residue(kk).indices));
            if diff==0,
                NOpos=model.sites{k}(kc).residue(kk).NOpos;
            end;
        end;
    end;
end;

