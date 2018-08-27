function [dmat,err]=metrize(lower_bounds,upper_bounds,silent,Gaussian)
% [dmat,err]=metrize(lower_bounds,upper_bounds,silent,Gaussian)
%
% Generates a pseudo-random metric space (distance matrix dmat), so that
% all distances lie within their respective bounds
%
% in Gaussian mode, lower and upper bounds are interpreted as being -/+ one
% standard deviation
%
% error output to the command window is suppressed if input argument silent 
% is provided and is not zero
%
% algorithm as in G. M. Crippen, T. F. Havel, Distance Geometry and
% Molecular Conformation, Research Studies Press, Taunton, 1988, p. 253 f.
%
% Error codes err:
% 0   no error
% 1   erroneous bounds 
% 2   at least one of the bound matrices is not square
% 3   bound matrices do not have the same size
%
% (c) G. Jeschke 2008, 2017 (Gaussian mode)

if nargin<3 || isempty(silent),
    silent=0;
end;

if nargin<4 || isempty(Gaussian),
    Gaussian=false;
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
        if err, return; end;
        if Gaussian
            dmat(k,kk)=random_n(lower_bounds(k,kk),upper_bounds(k,kk));
        else
            dmat(k,kk)=random(lower_bounds(k,kk),upper_bounds(k,kk));
        end
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

function x = random_n(xmin,xmax)
rmean = (xmin+xmax)/2;
stdr = xmax-rmean;
x = rmean + stdr*randn;
if x < xmin,
    x = xmin;
end;
if x > xmax
    x = xmax;
end;
        