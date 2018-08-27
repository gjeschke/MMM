function [ kappa ] = mode_collectivity( modes , masses )
%MODE_COLLECTIVITY Computes the collectivity of a normal mode
%   conforms to R. Bruschweiler, J. Chem. Phys. 1995, 102, 3396-3403
%
% modes     normal mode vector or matrix of normal mode vectors (modes are then
%           of the form u(:,k)), must have dimension [m,n] with m being a
%           multiple of 3
% masses    (optional) vector of nuclear masses, must have length m/3 if
%           present, otherwise output is empty
% kappa     collectivity of each input mode
%
% 
%
% G. Jeschke, 2011

kappa=[];

[m,n]=size(modes);

if m-3*floor(m/3)>eps,
    fprintf(2,'Warning: Normal mode vector length not a multiple of three in mode_collectivity.\n');
    fprintf(2,'Empty output returned.\n');
    return
end;

if nargin>1,
    N=length(masses);
    if N~=m/3,
        fprintf(2,'Warning: Mass vector and normal mode length do not macht in mode_collecitivity\n');
        fprintf(2,'Empty output returned.\n');
        return
    else
        [mm,nm]=size(masses);
        if nm<mm,
            masses=masses';
        end;
    end;
else
    masses=ones(1,m/3);
    N=m/3;
end;

kappa=zeros(1,n);
for kk=1:n
    mode=reshape(modes(:,kk),3,m/3);
    uinsq=sum(mode.^2,1)./masses;
    uinsq=uinsq/sum(uinsq);
    kappa(kk)=exp(-sum(uinsq.*log(uinsq)))/N;
end

