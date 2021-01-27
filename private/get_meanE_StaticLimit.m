function E_FRET = get_meanE_StaticLimit(R0,Pr)

% function E_FRET = getMeanEFastLimit(R0,Pr)
% return FRET efficiency in the static/rigid limit
%
% input:  R0       Foerster Radius in nm
%         Pr       Distance distribution P(r) in nm, must be size nx2
% output: E_FRET   expected FRET efficiency
% 
% Daniel Klose, 2019

[~, s] = size(Pr);
if s ~= 2
    error('P(r) has wrong size, must be nx2.');
end
zero_ind = bitand(Pr(:,1) ~= 0, Pr(:,2) >= 1e-4);
Pr = Pr(zero_ind,:);
Pr(:,2) = Pr(:,2)./trapz(Pr(:,1),Pr(:,2));

int = 1./(1+ (Pr(:,1)./R0).^6).*Pr(:,2);
E_FRET = trapz(Pr(:,1),int);

end