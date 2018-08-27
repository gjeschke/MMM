function [protein,chain_tags,chain_ids]=test_min_Ca_dist(pdbid)

[protein,chain_tags,chain_ids]=rd_pdb_coarse(pdbid);

maxn=6;
min_distances=1e6*ones(1,maxn);
max_distances=zeros(1,maxn);
mean_distances=zeros(1,maxn);
count_distances=zeros(1,maxn);
nc=length(protein);
nrt=0;
for k=1:nc,
    nr=length(protein(k).resnum);
    nrt=nrt+nr;
    for kk=1:nr-1,
        for kkk=kk+1:nr,
            dc=abs(protein(k).resnum(kk)-protein(k).resnum(kkk));
            if dc<=maxn,
                r=norm(protein(k).Ca(kk,:)-protein(k).Ca(kkk,:));
                mean_distances(dc)=mean_distances(dc)+r;
                count_distances(dc)=count_distances(dc)+1;
                if min_distances(dc)>r,
                    min_distances(dc)=r;
                end;
                if max_distances(dc)<r,
                    max_distances(dc)=r;
                end;
            end;
        end;
    end;
end;

mean_distances=mean_distances./count_distances;

fprintf(1,'--- Protein structure %s with %i chains and a total of %i residues. ---\n',pdbid,nc,nrt);
for k=1:maxn,
    fprintf(1,'Min. geometric distance for sequence distance %i is %4.2f Å\n',k,min_distances(k));  
end;
for k=1:maxn,
    fprintf(1,'Mean geometric distance for sequence distance %i is %4.2f Å\n',k,mean_distances(k));  
end;
for k=1:maxn,
    fprintf(1,'Max. geometric distance for sequence distance %i is %4.2f Å\n',k,max_distances(k));  
end;
