function test_magic_fit

angles=2*pi*rand(1,3);
transmat1=affine('Euler',angles);
translate=10*(rand(1,3)-0.5);
transmat2=affine('translation',translate);

transmat=transmat1*transmat2;

coor1=30*(rand(20,3)-0.5);
coor1b=coor1+0.5*(rand(20,3)-0.5);
d0=coor1b-coor1;
d0=sum(d0.^2,2);
d0=ones(size(d0),1)./d0;
d0=d0/sum(d0);

[mm,~]=size(coor1b);
coor2=[coor1b ones(mm,1)];
coor2=transmat*coor2';
coor2=coor2';
coor2=coor2(:,1:3);

[rms,coor2b,transmat_fit]=rmsd_superimpose(coor1,coor2);

[rms,coor2b,transmat_fit,weights,nit]=rmsd_superimpose_weighted(coor1,coor2);

fprintf(1,'R.m.s.d. = %6.4f\n',rms);
fprintf(1,'Max. dev. of coordinates = %8.6f\n',max(max(abs(coor2b-coor1b))));
test_tm=transmat*transmat_fit-eye(4);
fprintf(1,'Max. dev. of transformation matrix = %6.4f\n',max(max(abs(test_tm))));
fprintf(1,'Number of iterations: %i\n',nit);


figure(1); clf;
plot(d0,weights,'.');



