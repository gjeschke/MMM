% [Km, rlist] = discretizeDiffusionInMeanPotential(diff, peq, ri, rf, nR)
% discretizes diffusion problem
%
%
function [Km, rlist] = discretize_Diffusion_inMean_Potential(diff, peq, ri, rf, nR)

step = (rf - ri)./nR; % size of a step; set outer boundary to r=3a

x = (1:nR)';

% reduce peq to 1D, fewer points by interpolation
peq = interp1(peq(:,1),peq(:,2),linspace(ri,rf,nR)');

for i = 2:nR-1
    lp(i) = diff*(peq(i)+peq(i+1))./2./(step^2 * peq(i));
    lm(i) = diff*(peq(i)+peq(i-1))./2./(step^2 * peq(i));
end
lp(1) = diff*(peq(1)+peq(1+1))./2./(step^2 * peq(1));
lp(nR) = 0; % last element
lm(1) = 0; % first element
lm(nR) = diff*(peq(nR)+peq(nR-1))./2./(step^2 * peq(nR));

Km = zeros(nR, nR);
for i = 1:nR  
    % diagonal element
    Km(i,i) = -(lp(i) + lm(i));
    if i < nR
        Km(i,i+1) = lm(i+1);
    end
    if i > 1
        Km(i,i-1) = lp(i-1);
    end
end

rlist = 1:nR;
rlist = arrayfun(@rfunc, rlist);



    % local function
    function x = rfunc(a)

    x = ri + step.*(a-0.5);

    end

end