function ang=dihedral(c1,c2,c3,c4)
% function ang=dihedral(c1,c2,c3,c4)
% function ang=dihedral(v1,v2,v3)
%
% Compute dihedral angle for four points with coordinates c1,c2,c3,c4 (four
% input arguments or three vectors v1, v2, v3 (three input arguments)
%
% ang   dihedral angle in units of radians
%
% (c) 2000 Gunnar Jeschke
%

if nargin>3,
    vec1=c1-c2;
    vec2=c4-c3;
    vec3=c2-c3;
else
    vec1=c1;
    vec2=c2;
    vec3=c3;
end;
vec3=vec3/norm(vec3);
vec1=vec1-sum(vec1.*vec3)*vec3;
vec2=vec2-sum(vec2.*vec3)*vec3;
vec1=vec1/norm(vec1);
vec2=vec2/norm(vec2);
vec4=cross(vec1,vec2);
if norm(vec4)<1e-6,
    sign0=1;
else
	vec4=vec4/norm(vec4);
    sign0=-sum(vec4.*vec3);
end;
ang=sign0*acos(sum(vec1.*vec2));
