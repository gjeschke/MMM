function [pairs,popsum]=pairs_in_range(range,NOpos1,NOpos2)
% function pairs=pairs_in_range(range,NOpos1,NOpos2)
% analyzes which rotamer pairs contribute to distance distribution in a
% certain range
%
% range  - distance range (nm)
% NOpos1 - Nx5 array; columns 1 to 3 - XYZ coordinates for the NO bond midpoints;
%                     column 4 - population of the corresponding rotamer                     
%                     column 5 - rotamer number
%                     column 5 can be missing, in this case rotamer number
%                     is index (row number)
% NOpos2 - Mx5 array; the same. N and M can be different.
%
% pairs  - K*3 array; K is number of pairs contibuting in range
%                     column 1 - rotamer number for site 1
%                     column 2 - rotamer number for site 2
%                     column 3 - population of this pair
% popsum - double, sum of pair populations for the whole distribution 
%
% G. Jeschke, 2010
% adapted from get_distribution by Ye. Polyhach, 2009


%warning off MATLAB:warning:FrequencyOutputObsolete;

[n_rot_1,m]=size(NOpos1);   % returns number of rotamers at position 1
if m<5,
    indy=1:n_rot_1;
    NOpos1(:,5)=indy';
end;
[n_rot_2,m]=size(NOpos2);   % -//- at position 2
if m<5,
    indy=1:n_rot_2;
    NOpos2(:,5)=indy';
end;


poi=0;
popsum=0;
pairs=zeros(100000,3);
for k=1:n_rot_1
    NO1=NOpos1(k,1:3); % NO center for k-th rotamer in the pos. 1
    pop1=NOpos1(k,4); % weight for k-th rotamer in the pos. 1       
    for kk=1:n_rot_2
        NO2=NOpos2(kk,1:3); % NO center for k-th rotamer in the pos. 1
        pop2=NOpos2(kk,4);    % weight for k-th rotamer in the pos. 1
        popsum=popsum+pop1*pop2;
        NO12=((sum((NO1-NO2).^2))^(1/2))/10; % NOpos1-NOpos2 distance (dipolar distance for k-kk pair) in nm!
        if NO12>=range(1) && NO12<=range(2),
            poi=poi+1;
            pairs(poi,1)=NOpos1(k,5);
            pairs(poi,2)=NOpos2(kk,5);
            pairs(poi,3)=pop1*pop2;
        end
    end
end

pairs=pairs(1:poi,:);