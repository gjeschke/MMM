function [rax,distr]=get_distribution_oligomer(NOpos,n,sig,threshold,normalize,re,options)

% NOpos - Nx4 array; columns 1 to 3 - XYZ coordinates for the NO bond midpoints;
%                     coloumn 4 - population of the corresponding rotamer
% n     multiplicity of the oligomer
% sig - width of Gaussian used to convolute resulting distribution,
%       defaults to 0.001 (no significant broadening)
% threshold - minimum population, below which rotamers are not considered
%            (not recommended in general case, default value 0)
% normalize - (optional) flag, true: distribution is normalized, false: it
%             is not, defaults to true
% re        - upper distance limit
% options   - options structure
%             .trajectory array of two flags, if both are true, NOpos1 and
%             NOpos2 are considered as coordinates of subsequent trajectory
%             frames, in that case only distances from the same frames are
%             combined
%
% Ye. Polyhach, 2009

distr = [];

if nargin<3 || isempty(sig)
    sig=0.001; % no significant broadening
    threshold=0;
end;
if nargin<4 || isempty(threshold)
    threshold=0;
end;
if nargin<5 || isempty(normalize)
    normalize=true;
end;


%warning off MATLAB:warning:FrequencyOutputObsolete;

ra=-5;
if nargin<6 || isempty(re),
    re=20;
    nr=501;
else
    re=re-ra;
    nr=round(20*(re-ra))+1;
end;

if nargin<7 || isempty(options),
    options.trajectory=[0,0];
end;

rax=linspace(ra,re,nr);

if nargin<2,
    rax=rax(101:end-100);
    return;
end;

vari=rax/sig;
vari=exp(-vari.^2); % broadening function

distr=zeros(size(rax));

[n_rot,~]=size(NOpos);   % returns number of rotamers at position 1

for k=1:n_rot
    if NOpos(k,4)>=threshold
        NOr=sqrt(sum(NOpos(k,1:2).^2)); % NO center radial distance for k-th rotamer, assuming z axis is symmetry axis 
        r = 2*NOr*sin(pi/n)/10;
        pop=NOpos(k,4); % weight for k-th rotamer in the pos. 1

        poi=1+round(nr*(r-ra)/(re-ra));
        %distr(poi)=distr(poi)+1;
        if poi<=nr && poi>=1,
        %keyboard
            distr(poi)=distr(poi)+pop;
        end;
    end
end

% Convolution with broadening function
inv_distr=ifft(distr).*ifft(vari);
distr=real(fft(inv_distr));

rax=rax(101:end-100);
%rax=rax(101:401);
distr=distr(201:end);
%distr=distr(201:501);

if normalize,
    if sum(distr)>1e-12*length(distr),
        distr=distr/sum(distr);
    else
        distr=1e-12*ones(size(distr));
    end;
end;
