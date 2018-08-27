function [point, sos, singular] = multilaterate(ref_points,dist)
% point = multilaterate(ref_points,distances)
%
% Determination of point coordinates by multilateration from three
% (trilateration) or more reference points and uncertain distances to these
% reference points
% 
% Input:
% ref_points    array [np,3] of Cartesian reference point coordinates
% dist          vector of the np distances to the reference points
% 
% Output:
% point         [1,3] row vector of Cartesian coordinates of the
%               multilaterated point, empty output for less than three
%               reference points, for three reference points, a [2,3] array
%               of two solutions or an empty array (no solution) may be
%               returned
% sos           sum of square distance errors in the multilateration case,
%               zero for successful trilateration, negative for unsuccesful
%               trilateration, zero for less than two reference points
% singular      flag that shows if a poorly conditioned matrix is
%               encountered in mutilateration (true) or not (false), this
%               happens for instance for multilateration from 4 points
%               within a plane, as far as I see this does not pose problems
%               for this particular task
%
% G. Jeschke, 25.4.2013

% the following uses the linearization of the problem by W. S. Murphy Jr.,
% Master thesis, Colorado School of Mines, 2007 for obtaining a starting
% coordinate for the standard non-linear least squares solution of the
% multilateration problem

testit = true;

sos = 0;
point = [];
singular = false;

[np,~] = size(ref_points);

if np<3,
    return;
end;

if np == 3 && ~testit, % trilateration, actually based on Wikipedia article
    % go to standard frame
    origin = ref_points(1,:);
    coor = ref_points - repmat(origin,3,1);
    xe = coor(2,:)/norm(coor(2,:)); % unit vector along new x direction
    xye = coor(3,:)/norm(coor(3,:)); % another unit vector in xy plane
    if abs(sum(xe.*xye))>acos(5*pi/180), % points are almost collinear
        sos = -1;
        return;
    end;
    ze = cross(xe,xye); % unit vector along new z direction
    ze = ze/norm(ze); % defensive programming
    ye = cross(ze,xe);
    ye = ye/norm(ye);
    Rp = [xe;ye;ze];
    coor1 = coor*Rp'; % coordinates in standard frame, point 1 origin, 
                      % point 2 on x axis, point 3 in xy plane
    d = coor1(2,1);
    i = coor1(3,1);
    j = coor1(3,2);
    x = (dist(1)^2 - dist(2)^2+d^2)/(2*d);
    y = (dist(1)^2 - dist(3)^2 + i^2 + j^2)/(2*j) - i*x/j;
    z = sqrt(dist(1)^2 - x^2 - y^2);
    if imag(z) > 10*eps, % no solution
        sos = -2;
        return;
    elseif z < 10*eps, % single solution
        point = [x,y,z]*Rp + origin;
    else % two solutions
        p1 = [x,y,real(z)]*Rp + origin;
        p2 = [x,y,-real(z)]*Rp + origin;
        point = [p1;p2];
    end;
elseif np == 3 && testit,
        % go to standard frame
    origin = ref_points(3,:);
    coor = ref_points - repmat(origin,3,1);
    xe = coor(1,:)/norm(coor(1,:)); % unit vector along new x direction
    xye = coor(2,:)/norm(coor(2,:)); % another unit vector in xy plane
    if abs(sum(xe.*xye))>acos(5*pi/180), % points are almost collinear
        sos = -1;
        return;
    end;
    ze = cross(xe,xye); % unit vector along new z direction
    ze = ze/norm(ze); % normal vector to the plane
    A = zeros(np-1,3);
    b = zeros(np-1,1);
    for k=2:np,
        A(k-1,:) = ref_points(k,:)-ref_points(1,:); % Eq. (2.9)
        dk1 = norm(ref_points(k,:)-ref_points(1,:)); % Eq. (2.4)
        b(k-1) = (dist(1)^2 - dist(k)^2 + dk1^2)/2; % Eq. (2.8)
    end;

    if cond(A) < 1e10
        x = A\b;
    else
        x = pinv(A)*b;
        singular = true;
    end;
    pt_0 = x' + ref_points(1,:);

    [point1, sos] = fminsearch(@sos_distances,pt_0,[],ref_points,dist);
    point2 = point1 + 2*dot(ref_points(1,:) - point1, ze)*ze;
    point = [point1;point2];
else % multilateration
    A = zeros(np-1,3);
    b = zeros(np-1,1);
    for k=2:np
        A(k-1,:) = ref_points(k,:)-ref_points(1,:); % Eq. (2.9)
        dk1 = norm(ref_points(k,:)-ref_points(1,:)); % Eq. (2.4)
        b(k-1) = (dist(1)^2 - dist(k)^2 + dk1^2)/2; % Eq. (2.8)
    end;

    if cond(A) < 1e10,
        x = A\b;
    else
        x = pinv(A)*b;
        singular = true;
    end;
    pt_0 = x' + ref_points(1,:);

    [point, sos] = fminsearch(@sos_distances,pt_0,[],ref_points,dist);
end;

function sos = sos_distances(pt,ref_points,distances)

[np,~] = size(ref_points);
diff = ref_points - repmat(pt,np,1);
dist_found = sqrt(sum(diff.^2,2));
sos = sum((dist_found-distances).^2);
