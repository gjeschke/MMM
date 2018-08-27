function [deer,bckg,bckg_k,mod_depth,dim,rmsd]=fit_primary(handles,texp,vexp,ff)
% [deer,bckg,bckg_k,mod_depth]=fit_primary(handles,texp,vexp,ff)
%
% Fits the primary DEER data defined by time axis texp and data trace vexp
% with a given simulated form factor ff by varying modulation depth and
% background decay constant
% for a stretched exponential background (model 'fractional'), background
% dimension is also varied
% the background model is specified in handles:
%
% deer  simulated primary data (fit to vexp) 

rmsd=1e6;
fit_depth=get(handles.checkbox_mod_depth,'Value');


contents = get(handles.popupmenu_bckg,'String');
bckg_mode=contents{get(handles.popupmenu_bckg,'Value')};

depth0=handles.mod_depth; % starting value for modulation depth
if ~isempty(handles.exp_depth),
    depth0=handles.exp_depth;
end;
        
[poly,s]=polyfit(texp,log(vexp),1); % linear fit of logarithm
v0=[depth0 -poly(1)];
switch bckg_mode,
    case 'fractal',
        if fit_depth
            v0=[depth0 -poly(1) 3];
            [v1,rmsd]=fminsearch(@rms_stretched_MMM,v0,[],texp,vexp,ff);
            mod_depth=v1(1);
            v1=v1(2:3);
        else
            v0=[-poly(1) 3];
            [v1,rmsd]=fminsearch(@rms_stretched_MMM,v0,[],texp,vexp,ff,depth0);
            mod_depth=depth0;
        end;
        bckg=(1-mod_depth)*decay_stretched_MMM(v1,texp);
        bckg_k=v1(1);
        dim=v1(2);
    case '2D',
        distr=handles.Pake_r; % 2D distribution has probability proportional to r
        distr=0.01*distr/sum(distr);
        logB=distr*handles.Pake_kernel;
        if fit_depth,
            [v1,rmsd]=fminsearch(@rmsnD_MMM,v0,[],texp,vexp,ff,handles.Pake_t,logB);
            mod_depth=v1(1);
            v1=v1(2);
        else
            v0=v0(2);
            [v1,rmsd]=fminsearch(@rmsnD_MMM,v0,[],texp,vexp,ff,handles.Pake_t,logB,depth0);
            mod_depth=depth0;
        end;
        bckg=(1-mod_depth)*decaynD_MMM(v1,texp,handles.Pake_t,logB);
        bckg_k=v1;
        dim=2;
    case '3D',
        distr=handles.Pake_r.^2; % 2D distribution has probability proportional to r
        distr=0.01*distr/sum(distr);
        logB=distr*handles.Pake_kernel;
        if fit_depth,
            [v1,rmsd]=fminsearch(@rmsnD_MMM,v0,[],texp,vexp,ff,handles.Pake_t,logB);
            mod_depth=v1(1);
            v1=v1(2);
        else
            v0=v0(2);
            [v1,rmsd]=fminsearch(@rmsnD_MMM,v0,[],texp,vexp,ff,handles.Pake_t,logB,depth0);
            mod_depth=depth0;
        end;
        bckg=(1-mod_depth)*decaynD_MMM(v1,texp,handles.Pake_t,logB);
        bckg_k=v1;
        dim=3;
    otherwise, % this includes background none
        if fit_depth,
            [mod_depth,rmsd]=fminsearch(@rmsd_depth,depth0,[],vexp,ff);
        else
            mod_depth=depth0;
            sim=ff*mod_depth+ones(size(ff))*(1-mod_depth);
            diff=sim-vexp;
            rmsd=sqrt(diff.^2/length(diff));
        end;
        bckg_k=0;
        bckg=(1-mod_depth)*ones(size(vexp));
        dim=3;
end;
bckg=real(bckg);

deer=ff*mod_depth+ones(size(ff))*(1-mod_depth);
deer=bckg.*deer/max(bckg);
