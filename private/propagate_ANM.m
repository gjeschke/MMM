function [fom,network,DEER] = propagate_ANM(coeff,modes,network0,DEER,coarse)
% function [fom,network,DEER] = propagate_ANM(coeff,modes,network0,DEER,coarse)
%
% propagates an anisotropic network model along a linear combination of a
% few selected modes and computes a figure of merit (fom) with respect to
% DEER restraints
%
% Input:
% coeff     vector of linear coefficients for the normal modes along which
%           the network is to be propagated
% modes     vector of mode numbers corresponding to coeff, the numbers
%           refer to the anisotropic network model of the current structure
%           in MMM: model.ANM(model.current_structure).u
% network0  matrix [n,3] of coordinates of network points (knots),
%           restraint variable DEER must refer labelled sites to network
%           knots (C_alpha coordinates of residues), these references are
%           given for the k-th restraint by DEER(k).res1 and DEER(k).res2
%           indexing of residues in network is supposed to conform to
%           model.coarse(model.current_structure).indices
% DEER      restraint variable, vector of n structures, propagate_ANM uses
%           the fields
%           DEER(k).res1    (read)
%           DEER(k).res2    (read)
%           DEER(k).xyz1    (read/write)
%           DEER(k).xyz2    (read/write)
%           DEER(k).r       (read)
%           DEER(k).sigr    (read)
% coarse    flag that requests fast, coarse restraint computation that is
%           a good approximation for small steps, optional, defaults to
%           false
%
% Output:
% fom       figure of merit, quantifies deviation of the propagated model
%           from restraints, maximum likelihood estimated based on mean
%           distance and standard deviation of mean distance for all
%           restraints
% network   the propagated network
% DEER      the DEER restraints variable with updated label coordinates
%           DEER(k).xyz1 and DEER(k).xyz2
%           
% G. Jeschke, 2010

global model

if nargin<5,
    coarse=false;
end;

fom = 1e6;
m=model.ANM(model.current_structure).residues;

% network propagation
network=network0; % initialize propagated network
for k=1:length(coeff),
    cmode=modes(k);
    evec=model.ANM(model.current_structure).u(:,cmode);
    mode=reshape(evec,3,m);
    network=network+coeff(k)*mode';
end;

[maxnum,m]=size(network);
scarce=0;
% update of label coordinates
for k=1:length(DEER),
    xyz=DEER(k).xyz1;
    l=DEER(k).res1;
    cindices=model.coarse(model.current_structure).indices; % actual indices of residues in network
    local_template=zeros(5,3);
    local_template_0=zeros(5,3);
    % make a local template to fit rotation and translation
    poi=0;
    if ~coarse,
        for kk=-2:2,
            if l+kk>0 && l+kk<=maxnum, % is addressed residue a network point?
                diff=cindices(l+kk,4)-cindices(l,4);
                if diff==kk, % is addressed residue part of a continuous segment?
                    poi=poi+1;
                    local_template(poi,:)=network(l+kk,:);
                    local_template_0(poi,:)=network0(l+kk,:);
                end;
            end;
        end;
    end;
    if poi>=3, % found sufficient number of points to determine local rotation and translation
        [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz=[xyz 1];
        xyz=transmat*xyz';
        xyz=xyz';
        DEER(k).xyz1=xyz(1:3);
    else
        DEER(k).xyz1=xyz+network(l,:)-network0(l,:);
        scarce=scarce+1;
    end;
    xyz=DEER(k).xyz2;
    l=DEER(k).res2;
    cindices=model.coarse(model.current_structure).indices; % actual indices of residues in network
    local_template=zeros(5,3);
    local_template_0=zeros(5,3);
    % make a local template to fit rotation and translation
    poi=0;
    if ~coarse,
        for kk=-2:2,
            if l+kk>0 && l+kk<=maxnum, % is addressed residue a network point?
                diff=cindices(l+kk,4)-cindices(l,4);
                if diff==kk, % is addressed residue part of a continuous segment?
                    poi=poi+1;
                    local_template(poi,:)=network(l+kk,:);
                    local_template_0(poi,:)=network0(l+kk,:);
                end;
            end;
        end;
    end;
    if poi>=3, % found sufficient number of points to determine local rotation and translation
        [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz=[xyz 1];
        xyz=transmat*xyz';
        xyz=xyz';
        DEER(k).xyz2=xyz(1:3);
    else
        DEER(k).xyz2=xyz+network(l,:)-network0(l,:);
        scarce=scarce+1;
    end;
end;

fom=0; % figure of merit
for k=1:length(DEER),
    r0=norm(DEER(k).xyz1-DEER(k).xyz2)/10;
    det=(r0-DEER(k).r)/DEER(k).sigr;
    fom=fom+det^2;
end;
fom=sqrt(fom/length(DEER));

model.ANM(model.current_structure).u=reorient_ANM(model.ANM(model.current_structure).u,network0,network);

% disp(sprintf('For %i labeled sites C_alpha-N-O vector rotation could not be specified.',scarce));