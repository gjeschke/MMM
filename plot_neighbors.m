function plot_neighbors(pdbid)

[protein,chain_tags,chain_ids]=rd_pdb_coarse(pdbid);

ra=1;
re=10;
dr=0.01;
rax=ra:dr:re;
nr=length(rax);
first=zeros(size(rax));
second=zeros(size(rax));

for k=1:length(protein),
    Ca=protein(k).Ca;
    resnum=protein(k).resnum;
    for kk=1:length(resnum),
        res1=resnum(kk);
        res2=res1+1;
        poi=find(resnum==res2,1);
        if ~isempty(poi),
            rvec=Ca(kk,:)-Ca(poi,:);
            r=norm(rvec);
            poi=1+round((r-ra)/dr);
            if poi>=1 && poi<=nr,
                first(poi)=first(poi)+1;
            else
                fprintf(2,'Direct neighbor distance outside histogram\n');
            end;
        end;
        res3=res1+2;
        poi=find(resnum==res3,1);
        if ~isempty(poi),
            rvec=Ca(kk,:)-Ca(poi,:);
            r=norm(rvec);
            poi=1+round((r-ra)/dr);
            if poi>=1 && poi<=nr,
                second(poi)=second(poi)+1;
            else
                fprintf(2,'Next neighbor distance outside histogram\n');
            end;
        end;
    end;
end;

first=first/sum(first);
second=second/sum(second);

figure(1); clf;
plot(rax,first,'b');
hold on;
plot(rax,second,'r');

mean1=sum(first.*rax);
mean2=sum(second.*rax);

fprintf(1,'Direct neighbor mean distance %4.2f Å\n',mean1);
fprintf(1,'Next neighbor mean distance %4.2f Å\n',mean2);
