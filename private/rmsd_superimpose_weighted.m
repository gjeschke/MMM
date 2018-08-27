function [rms,coor2b,transmat,weights,nit]=rmsd_superimpose_weighted(coor1,coor2)

transmat0=zeros(4,4);
[rms,coor2b,transmat]=rmsd_superimpose(coor1,coor2);

nit=0;
while max(max(abs(transmat-transmat0)))>1e-4 && nit<1000,
    nit=nit+1;
    transmat0=transmat;
    diff=coor2b-coor1;
    weights=sum(diff.^2,2);
    weights=ones(size(weights))./weights;
    weights=weights/sum(weights);
    [rms,coor2b,transmat]=rmsd_superimpose(coor1,coor2,weights);
end;