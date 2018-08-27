function ang=angle(c1,c2,c3)
% function ang=angle(c1,c2,c3)
% function ang=angle(v1,v2)
%
% Compute angle for three points with coordinates c1,c2,c3 (three input 
% arguments) or two vectors v1, v2 (two arguments)
% 
% ang   angle in radians
%
% (c) 2009 Gunnar Jeschke
%

if nargin>2,
    v1=c1-c2;
    v2=c3-c2;
else
    v1=c1;
    v2=c2;
end;

v1=v1/norm(v1);
v2=v2/norm(v2);
ang=acos(sum(v1.*v2));
