function [rax,dist]=rd_PRONOX(fname)

sig=0.1; % Gaussian broadening of the distance distribution

ra=-5;
re=20;
nr=501;
rax=linspace(ra,re,nr);
dist=zeros(1,length(rax));
vari=rax/sig;
vari=exp(-vari.^2); % broadening function

key='nitroxide - nitroxide distances(A)';
key0='SSlabel.zmat';
keyfav='t1,t2 fav';
keyunfav='t1,t2 unfav';


fid=fopen(fname,'r');

ndat=0;
datamode=false;
confmode=false;
rvec=zeros(1,80000);
wvec=zeros(1,80000);
favconf=zeros(20000,2);
mconf=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    poi=strfind(tline,key0);
    if ~isempty(poi),
        confmode=true;
        fgetl(fid);
        tline = fgetl(fid);
    end;
    poi=strfind(tline,key);
    if ~isempty(poi),
        datamode=true;
        fgetl(fid);
        tline = fgetl(fid);
        favconf=favconf(1:mconf,:);
    end;
    if length(tline)<50,
        datamode=false;
        confmode=false;
    end;
    if datamode,
        lab1=str2double(tline(4:8));
        if favconf(lab1,1)>0,
            w1=0.9;
        elseif favconf(lab1,1)<0,
            w1=0.1;
        else
            w1=1;
        end;
        lab2=str2double(tline(20:24));
        if favconf(lab2,2)>0,
            w2=0.9;
        elseif favconf(lab2,2)<0,
            w2=0.1;
        else
            w2=1;
        end;
        r = str2double(tline(45:end))/10;
        ndat=ndat+1;
        rvec(ndat)=r;
        wvec(ndat)=w1*w2;
        poi=1+round(nr*(r-ra)/(re-ra));
        dist(poi)=dist(poi)+w1*w2;
    end;
    if confmode,
        lab=str2double(tline(1:6));
        conf=str2double(tline(7:12));
        if conf>mconf,
            mconf=conf;
        end;
        poi=strfind(tline,keyfav);
        if ~isempty(poi),
            favconf(conf,lab)=1;
        end;
        poi=strfind(tline,keyunfav);
        if ~isempty(poi),
            favconf(conf,lab)=-1;
        end;
    end;
end;
fclose(fid);

% rvec=rvec(1:ndat);
% wvec=wvec(1:ndat);
% rvec=rvec.*wvec/sum(wvec);
% fprintf(1,'Mean distance %4.2f Å\n',10*sum(rvec));
% fprintf(1,'Std. deviation %4.2f Å\n',10*std(rvec));
% fprintf(1,'%i conformer pairs\n',ndat);


inv_distr=ifft(dist).*ifft(vari);
distr=real(fft(inv_distr));

rax=rax(101:end-100);
%rax=rax(101:401);
dist=distr(201:end);