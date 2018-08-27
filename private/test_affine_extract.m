function test_affine_extract

transmat1 = affine('xyz',[27,47,113]*pi/180);
transmat2 = affine('translation',[-3.73,5.67,12.47]);
transmat = transmat2*transmat1;
[trans,euler] = affine2transrot(transmat);
euler = 180*euler/pi;

fprintf(1,'Rextracted Euler angles: %6.2f, %6.2f, %6.2f degree\n',euler);
fprintf(1,'Rextracted translation: %6.2f, %6.2f, %6.2f Å\n',trans);
