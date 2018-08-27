function [trax,tdistr]=tweak_rotamer_populations(adr1,adr2,rax,distr)
% Analyze which rotamer pairs contribute to distance distribution in the
% range specified by handles.range

sig=0.1;

NOpos1=get_NOpos(adr1);
NOpos2=get_NOpos(adr2);

% [rax0,distr0]=get_distribution(NOpos1,NOpos2,sig);

[mrot1,nr]=size(NOpos1);
[mrot2,nr]=size(NOpos2);

dist=zeros(1,mrot1*mrot2);
rmin=1000;
rmax=0;
pop=zeros(1,mrot1*mrot2);
poi=0;
for k=1:mrot1,
    xyz1=NOpos1(k,1:3);
    for kk=1:mrot2,
        poi=poi+1;
        xyz2=NOpos2(kk,1:3);
        dist(poi)=norm(xyz2-xyz1)/10;
        if dist(poi)<rmin,
            rmin=dist(poi);
        end;
        if dist(poi)>rmax,
            rmax=dist(poi);
        end;
        pop(poi)=NOpos1(k,4)*NOpos2(kk,4);
    end;
end;
rmin=min(rax);
rmax=max(rax);
nr=length(rax)-1;

pop0=zeros(size(pop));
popsum=0;
while sqrt(sum((pop0-pop).^2)/length(pop))>1e-12
    pop0=pop/sum(pop);

    % Tweak the rotamer distribution
    pop1=zeros(1,mrot1);
    pop2=zeros(1,mrot2);
    popsum=0;
    poi=0;
    for k=1:mrot1,
        for kk=1:mrot2,
            poi=poi+1;
            rr=dist(poi);
            rpoi=round(1+nr*(rr-rmin)/(rmax-rmin));
            if rpoi<1 || rpoi>nr+1,
                ppop=0;
            else
                ppop=distr(rpoi);
            end;
            pop1(k)=pop1(k)+pop0(poi)*sqrt(ppop);
            pop2(kk)=pop2(kk)+pop0(poi)*sqrt(ppop);
            popsum=popsum+pop(poi)*ppop;
            pop(poi)=ppop;
        end;
    end;
    pop1=pop1/sum(pop1);
    pop2=pop2/sum(pop2);
    pop=pop/sum(pop);
    % disp(sqrt(sum((pop0-pop).^2)/length(pop)));
end;

add_msg_board(sprintf('Selected rotamer pair merit is %6.5f\n',popsum));
        
NOpos1(:,4)=pop1';
NOpos2(:,4)=pop2';

[trax,tdistr]=get_distribution(NOpos1,NOpos2,sig);

% tdistr=interp1(trax0,tdistr0,rax,'pchip',0);
% % sc=sum(tdistr.*distr)/sum(tdistr.*tdistr);
% % tdistr=sc*tdistr;
% tdistr=tdistr*sum(distr)/sum(tdistr);
% trax=rax;
        
function NOpos=get_NOpos(adr)

global model

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
