maxdev = 0;
profile on
for k = 1:1000000
    coor = 30*rand(27,3);
    dmat1 = coor2dmat(coor);
    dmat2 = coor2dmatsq(coor);
    dmat3 = coor2dmatself(coor);
    maxdevc = max(max(abs(dmat1-dmat3)));
    if maxdevc > maxdev
        maxdev = maxdevc;
    end
end
fprintf(1,'Maximum deviation: %g\n',maxdev);
profile viewer
profile off