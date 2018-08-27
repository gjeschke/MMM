function [gobjects,x,y,z]=plot_ribbon(indices)
% function gobjects=plot_ribbon(indices)
%
% Plots ribbon model for a range of residues and returns
% array of graphics objects handles
%
% the handle graphics objects and graphics settings are stored in a
% structure variable residue_graphics for that residue with fields
% .objects  handles to Matlab graphics objects
% .mode     atom graphics mode, 0 none, 1 tube, 2 helical, 3 sheet-like
% .color    RGB color vector
% .opaque   opaqueness (1-transparency) 
%
% indices   indices to identify structure, coordinate set, chain, and residue
%
% gobjects  vector of handles to graphics objects of all alternate locations
%
% G. Jeschke, 2009

global model

gobjects=[];

[m,n]=size(indices);

plotted=zeros(100,5);
done=1;
for k=1:m,
    newflag=1;
    for kk=1:done,
        if sum(indices(1:3)==plotted(kk,1:3)) && indices(4)>= plotted(kk,4) && indices(4) <= plotted(kk,5) % check, whether this was already plotted
            newflag=0;
        end;
    end;
    if newflag
        nn=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);
        if indices(4)<=nn,
            % Determine range of residues with the same secondary structure
            sec=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).secondary; % secondary structure type
            range=indices(4);
            secn=sec;
            poi=indices(4);
            while poi>1 && secn==sec && ~model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(poi-1).hetflag,
                poi=poi-1;
                secn=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(poi).secondary;
                if secn==sec,
                    range=[poi range];
                end;
            end;
            secn=sec;
            poi=indices(4);
            while poi<length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info) && secn==sec && ~model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(poi+1).hetflag,,
                poi=poi+1;
                secn=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(poi).secondary;
                if secn==sec,
                    range=[range poi];
                end;
            end;
            if ~sum(plotted(1,:))
                plotted(1,1:3)=indices(1:3);
                plotted(1,4:5)=[range(1),range(end)];
            else
                done=done+1;
                plotted(done,1:3)=indices(1:3);
                plotted(done,4:5)=[range(1),range(end)];
            end;
            [objects,x,y,z]=plot_sec_element(indices(1:3),range,sec);
            gobjects=[gobjects objects];
        end;
    end;
end;

if isempty(gobjects),
    message.error=1;
    message.text='Nothing to plot.';
end;

