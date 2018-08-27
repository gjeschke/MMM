function [dmat,err]=network_metrize(lower_bounds,upper_bounds,dist,sigr,types,silent)
% [dmat,err]=metrize(lower_bounds,upper_bounds)
%
% Generates a pseudo-random metric space (distance matrix dmat), so that
% all distances lie within their respective bounds, mixed
% lower-bound/upper-bound and Gaussian distribution constraints are allowed
%
% error output to the command window is suppressed if input argument silent 
% is provided and is not zero
%
% algorithm as in G. M. Crippen, T. F. Havel, Distance Geometry and
% Molecular Conformation, Research Studies Press, Taunton, 1988, p. 253 f.
% with special handling for constraints given as Gaussian distributions
%
% lower_bounds  square matrix of lower distance bounds between network points
% upper_bounds  square matrix of upper distance bounds between network points
% dist          square matrix of mean distances for Gaussian constraints
% sigr          square matrix of standard deviations for Gaussian constraints
% types         square matrix of constraint types, 1 is Gaussian constraint
% silent        flag that suppresses error output to command window, if
%               present and true, defaults to true
%
% Error codes err:
% 0   no error
% 1   erroneous bounds 
% 2   at least one of the bound matrices is not square
% 3   bound matrices do not have the same size
%
% (c) G. Jeschke 2008-2013

if nargin<6,
    silent = true;
end;

err=0;
[m,n]=size(lower_bounds);
if m~=n,
    if ~silent, disp('### ERROR: Lower bound matrix not square. ###'); end;
    err=2;
    return;
end;
[m2,n2]=size(upper_bounds);
if m2~=n2,
    if ~silent, disp('### ERROR: Upper bound matrix not square. ###'); end;
    err=2;
    return;
end;
if n~=n2,
    if ~silent, disp('### ERROR: Bound matrices have different size. ###'); end;
    err=3;
    return;
end;

dmat=zeros(n,n);

for k=1:n-1,
    for kk=k+1:n,
%         disp(sprintf('%s%i%s%i%s','(',k,',',kk,')'));
%         if k==2 && kk==4
%             keyboard
%         end;
        [lower_bounds,upper_bounds,err]=triangle(lower_bounds,upper_bounds,silent);
        if err, 
            return; 
        end;
        if types(k,kk)==1,
            dmat(k,kk)=gauss_rand(dist(k,kk),sigr(k,kk));
        else
            dmat(k,kk)=random(lower_bounds(k,kk),upper_bounds(k,kk));
        end;
        dmat(kk,k)=dmat(k,kk);
        lower_bounds(k,kk)=dmat(k,kk);
        upper_bounds(k,kk)=dmat(k,kk);
        lower_bounds(kk,k)=dmat(kk,k);
        upper_bounds(kk,k)=dmat(kk,k);
    end;
end;

function x=random(xmin,xmax)
% provides a pseudo-random number conforming to a uniform distribution
% between xmin and xmax
x=xmin+rand*(xmax-xmin);

function x=gauss_rand(r,sr)
% provides a pseudo-random number conforming to a uniform distribution
% between xmin and xmax
x=r+sr*randn;