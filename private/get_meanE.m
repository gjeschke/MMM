function meanE = get_meanE(R0, kD, Pr, r1, r2, diff, nR)
% nR = 300; %256; % 50

zero_ind = bitand(Pr(:,1) ~= 0, Pr(:,2) >= 1e-4);
Pr = Pr(zero_ind,:);
r1 = max(r1, min(Pr(:,1)));
r2 = min(r2, max(Pr(:,1)));
% disp(['Dist. boundaries set to [' num2str(r1) ' ' num2str(r2) ']']);

[Kmc, rlist] = discretize_Diffusion_inMean_Potential(diff, Pr, r1, r2, nR);
M = diag(kD.*(1 + (R0./rlist).^6));
Km = Kmc - M;

peq = null(Kmc);
peq = peq./sum(peq);

ntimes = 275; % ntimes = 70;
times(1:ntimes,1) = 0.00005.*(1.05.^(0:ntimes-1)); % times(1:ntimes,1) = 0.0001.*(1.2.^(0:ntimes-1));
for i = 1:length(times)
    dat(i) = sum(expm(Km.*times(i))*peq);
end
dat = [times, dat'];
delta_t = 0.00005; % delta_t = 0.0001;
times_new = (delta_t:delta_t:30)'; %times_new = (delta_t:delta_t:24)';
dat_new = interp1(dat(:,1),dat(:,2),times_new);
% tauDAint = delta_t .* cumsum(dat_new);
% tauDAint = tauDAint(end);
tauDAint = trapz(times_new,dat_new);

meanE = 1 - kD*tauDAint;

end