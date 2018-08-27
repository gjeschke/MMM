data1=load('MtsslWizard.dat');
data2=load('Pronox_Hagelueken.dat');
data3=load('MMM_Hagelueken.dat');

% Test of consistency
diff12=data1(:,3)-data2(:,3);
fprintf(1,'Experimental difference 12: %4.2f\n',sum(sqrt(diff12.^2/length(diff12))));
diff13=data1(:,3)-data3(:,3);
fprintf(1,'Experimental difference 13: %4.2f\n',sum(sqrt(diff13.^2/length(diff13))));
diff23=data2(:,3)-data3(:,3);
fprintf(1,'Experimental difference 23: %4.2f\n',sum(sqrt(diff23.^2/length(diff23))));

diff1=data1(:,5)-data1(:,3);
diff2=data2(:,5)-data2(:,3);
diff3=data3(:,5)-data3(:,3);
fprintf(1,'MtsslWizard: mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff1),std(diff1));
fprintf(1,'Pronox     : mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff2),std(diff2));
fprintf(1,'MMM        : mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff3),std(diff3));

diff12=data1(:,5)-data2(:,5);
fprintf(1,'MtsslWizard/Pronox mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff12),std(diff12));
diff13=data1(:,5)-data3(:,5);
fprintf(1,'MtsslWizard/MMM    mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff13),std(diff13));
diff23=data2(:,5)-data3(:,5);
fprintf(1,'Pronox/MMM         mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff23),std(diff23));

T4L=1:36;
HCO=37:52;

T4L=1:32;
HCO=33:48;

fprintf(1,'\nT4 Lysozyme\n');
diff1=data1(T4L,5)-data1(T4L,3);
diff2=data2(T4L,5)-data2(T4L,3);
diff3=data3(T4L,5)-data3(T4L,3);
fprintf(1,'MtsslWizard: mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff1),std(diff1));
fprintf(1,'Pronox     : mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff2),std(diff2));
fprintf(1,'MMM        : mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff3),std(diff3));

fprintf(1,'\nHistone core octamer\n');
diff1=data1(HCO,5)-data1(HCO,3);
diff2=data2(HCO,5)-data2(HCO,3);
diff3=data3(HCO,5)-data3(HCO,3);
fprintf(1,'MtsslWizard: mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff1),std(diff1));
fprintf(1,'Pronox     : mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff2),std(diff2));
fprintf(1,'MMM        : mean dev. %4.2f nm, std. dev. %4.2f nm\n',mean(diff3),std(diff3));
