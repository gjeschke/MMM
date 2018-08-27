function profile_O2

global membrane_profiles

dax=linspace(-15,15,201);
% dax=linspace(1,16,16);
delta=zeros(1,201);

for k=1:length(dax),
    % delta(k)=delta_inv_T1(dax(k),membrane_profiles.O2.d0,membrane_profiles.O2.lambda,membrane_profiles.O2.mindel,membrane_profiles.O2.maxdel);
    delta(k)=delta_inv_T1(dax(k),membrane_profiles.O2.d0,membrane_profiles.O2.lambda,membrane_profiles.O2.mindel,membrane_profiles.O2.maxdel);
end;

figure(1); clf;
plot(dax,delta,'k');
