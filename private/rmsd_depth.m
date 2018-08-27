function rmsd=rmsd_depth(depth,vexp,ff)

if depth<0,
    rmsd=1e6;
    return
end;
if depth>1,
    rmsd=1e6;
    return
end;
sim=ff*depth+ones(size(ff))*(1-depth);
diff=sim-vexp;
rmsd=sqrt(sum(diff.^2)/length(diff));