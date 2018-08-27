function rmsd=rms_cluster(bsl,deer0,vexp)

n=length(vexp);
deer=deer0+bsl;
sc=sum(vexp.*deer)/sum(deer.*deer);
diff=sc*deer-vexp;
rmsd=sqrt(sum(diff.*diff)/(n-1));