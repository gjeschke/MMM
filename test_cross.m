function test_cross
mauz = rand(1,3);
miez = rand(1,3);
c1 = cross(mauz,miez);
c2 = cross_rowvec(mauz,miez);
disp(c1-c2);
