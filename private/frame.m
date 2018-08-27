function dircos=frame(p1,p2,p3)
% function dircos=frame(p1,p2,p3)
%
% Direction cosines of a right-handed frame defined by three points p1 
% (origin), p2 (point on the x axis), and p3 (point in the xy plane, but 
% not on the x axis);
% if p3 is nearly on the x axis (collinear points), an empty matrix is
% returned, this exception occurs for angles between vectors (p2,p1) and
% (p3,p1) smaller than 5° (0.087266 rad)
%
% dircos    direction cosines (3x3 rotation matrix)
%           row vectors are transformed to the frame by rvec=rvec0*dircos'
%           column vectors are transformed by cvec=dircos*cvec
%
% G. Jeschke, 2009

collinear=0.087266;

dircos=[];

vec1=p2-p1;
vec1=vec1/norm(vec1);
vec2=p3-p1;
vec2=vec2/norm(vec2);
if acos(sum(vec1.*vec2))<collinear,
    return
end;
vec3=cross(vec1,vec2);
vec3=vec3/norm(vec3);
vec2=cross(vec3,vec1);
vec2=vec2/norm(vec2);
dircos=[vec1;vec2;vec3];
