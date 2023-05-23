function [rax,distr]=get_distribution(NOpos1,NOpos2,sig,threshold,normalize,re,options)

% NOpos1 - Nx4 array; columns 1 to 3 - XYZ coordinates for the NO bond midpoints;
%                     coloumn 4 - population of the corresponding rotamer
% NOpos2 - Mx4 array; the same. N and M can be different.
%           if NOpos1 or NOpos2 are not provided, only the default distance
%           axis is returned together with an empty distribution
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
end
if nargin<4 || isempty(threshold)
    threshold=0;
end
if nargin<5 || isempty(normalize)
    normalize=true;
end


%warning off MATLAB:warning:FrequencyOutputObsolete;

ra=-5;
if nargin<6 || isempty(re)
    re=20;
    nr=501;
else
    re=re-ra;
    nr=round(20*(re-ra))+1;
end

if nargin<7 || isempty(options)
    options.trajectory=[0,0];
end

rax=linspace(ra,re,nr);

if nargin<2
    rax=rax(101:end-100);
    return;
end

vari=rax/sig;
vari=exp(-vari.^2); % broadening function

distr=zeros(size(rax));

[n_rot_1,~]=size(NOpos1);   % returns number of rotamers at position 1
[n_rot_2,~]=size(NOpos2);   % -//- at position 2

missing = 0;
% keyboard
if prod(options.trajectory)
    for k=1:n_rot_1
        NO1=NOpos1(k,1:3); % NO center for k-th frame in the pos. 1
        pop1=NOpos1(k,4); % weight for k-th frame in the pos. 1    end;
        NO2=NOpos2(k,1:3); % NO center for k-th frame in the pos. 2
        pop2=NOpos2(k,4);    % weight for k-th frame in the pos. 2
        NO12=((sum((NO1-NO2).^2))^(1/2))/10; % NOpos1-NOpos2 distance (dipolar distance for k-kk pair) in nm!
        if NO12>=1
            poi=1+round(nr*(NO12-ra)/(re-ra));
            if poi<=nr && poi>=1
                distr(poi)=distr(poi)+pop1*pop2;
            else
                missing = missing + pop1*pop2;
            end
        else
            missing = missing + pop1*pop2;
        end
    end
else
    for k=1:n_rot_1
        if NOpos1(k,4)>=threshold
            NO1=NOpos1(k,1:3); % NO center for k-th rotamer in the pos. 1
            pop1=NOpos1(k,4); % weight for k-th rotamer in the pos. 1

    %         disp([pop1,k]);

            for kk=1:n_rot_2
                if NOpos2(kk,4)>=threshold
                    NO2=NOpos2(kk,1:3); % NO center for k-th rotamer in the pos. 1
                    pop2=NOpos2(kk,4);    % weight for k-th rotamer in the pos. 1

                    NO12=((sum((NO1-NO2).^2))^(1/2))/10; % NOpos1-NOpos2 distance (dipolar distance for k-kk pair) in nm!
                    if NO12>=0.1
                        poi=1+round(nr*(NO12-ra)/(re-ra));
                        if poi<=nr && poi>=1
                            distr(poi)=distr(poi)+pop1*pop2;
                        else
                            missing = missing + pop1*pop2;
                        end
                    else
                        missing = missing + pop1*pop2;
                    end
                end
            end
        end
    end
end
% keyboard
% Convolution with broadening function
inv_distr=ifft(distr).*ifft(vari);
distr=real(fft(inv_distr));

rax=rax(101:end-100);
%rax=rax(101:401);
distr=distr(201:end);
%distr=distr(201:501);

if normalize
    if sum(distr)>1e-12*length(distr)
        in_histogram = sum(distr)/(missing + sum(distr));
        distr = in_histogram*distr/sum(distr);
    else
        distr=1e-12*ones(size(distr));
    end
end
