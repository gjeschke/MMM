function [rax,distr]=rd_MtsslWizard(fname)

sig=0.1;
ra=-5;
re=20;
nr=501;

rax=linspace(ra,re,nr);

fid=fopen(fname);
if fid==-1,
    distr=[];
    return;
% else % skip header line
%     fgetl(fid);
end;

distr=zeros(1,nr);

nl=0;
while 1
    tline = fgetl(fid);
    nl=nl+1;
    if ~ischar(tline), break, end
    bin=str2num(tline);
    r=bin(2)/10;
    poi=1+round(nr*(r-ra)/(re-ra));
    distr(poi)=distr(poi)+bin(3);
end;

vari=rax/sig;
vari=exp(-vari.^2); % broadening function
inv_distr=ifft(distr).*ifft(vari);
distr=real(fft(inv_distr));

rax=rax(101:end-100);
%rax=rax(101:401);
distr=distr(201:end);
distr=distr/sum(distr);
