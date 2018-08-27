function u=reorient_ANM(u,Ca_coor0,Ca_coor)
% function u=reorient_ANM(u,Ca_coor0,Ca_coor)
%
% provides (approximate) reoriented eigenvectors for an anisotropic network
% model whose knot coordinates have changed from Ca_coor0 to Ca_coor
%
% u         eigenvectors before (input) and after (output) coordinate
%           change as columns of a matrix
% Ca_coor0  knot coordinates before change
% Ca_coor   knot coordinates after change
%
% for algorithm, see:
% W. Zheng, B. R. Brooks, Biophys. J. 2006, 90: 4327-4336.
%
% G. Jeschke, 2010

% return
%  
% Hessian=setup_ANM_poly(Ca_coor);
Hessian=setup_ANM_bonded(Ca_coor);
% opts.issym=1;
% opts.isreal=1;
% opts.tol=1e-12;
% opts.disp=0;
% [u0,D]=eigs(Hessian,17,'sm',opts);
% u(:,1:17)=u0;
[u,D]=eig(Hessian);
return

% return

[m,n]=size(Ca_coor0);
rotmats=zeros(m,3,3); % rotation matrices for all residues
SRN0=zeros(m,3); % set of rigid neighbors, old coordinates
SRN=zeros(m,3); % set of rigid neighbors, new coordinates
ndist=zeros(1,m);
for k=1:m,
    r1=Ca_coor0(k,:);
    nn=0; % number of neighbors within 10 Å
    for j=1:m,
        r2=Ca_coor0(j,:);
        r=norm(r1-r2);
        if r<=10,
            nn=nn+1;
            ndist(nn)=r;
            SRN0(nn,:)=r2;
            SRN(nn,:)=Ca_coor(j,:);
        end;
        if nn>=3,
            pp=nn;
            if pp>4, pp=4; end;
            [y,poi]=sort(ndist(1:nn));
            coor0=SRN0(poi(1:pp),:);
            coor=SRN(poi(1:pp),:);
            % [rms,coor2,transmat]=rmsd_superimpose(coor,coor0);
            rotmat=opt_rot(coor,coor0);
            rotmats(k,:,:)=rotmat;
        else
            rotmats(k,:,:)=eye(3);
        end;
    end;    
end;
for k=1:m,
    basnum=3*(k-1);
    rotmat=reshape(rotmats(k,:,:),3,3);
    coor=u(basnum+1:basnum+3,:);
    u(basnum+1:basnum+3,:)=rotmat*coor;
end;
