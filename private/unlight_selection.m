function unlight_selection
% resets highlight color to actual color for current selection

global model

indices=resolve_address('*');
if isempty(indices),
    return
end;

[m,n]=size(indices);

for ko=1:m, % loop over all objects
    idepth=length(find(indices(ko,:)>0)); % determine type of current object
    cindices=indices(ko,1:idepth);
    allindices=cindices;
    if length(cindices)<4,
        if length(cindices)<2,
            adr=mk_address(cindices(1));
            adr=sprintf('%s(:){:}:',adr);
            allindices=resolve_address(adr);
        elseif length(cindices)<3,
            adr=mk_address(cindices(1:2));
            adr=sprintf('%s{:}:',adr);
            allindices=resolve_address(adr);
        else
            adr=mk_address(cindices(1:3));
            adr=sprintf('%s:',adr);
            allindices=resolve_address(adr);
        end;
    end;
    if length(cindices)==4,
        [msg,allindices]=get_residue(cindices,'descendants');
    end;
    [ma,na]=size(allindices);
    for ka=0:ma,
        if ka>0, % ka=0 takes care of the (un)selected object, k>0 of its descendants
            cindices=allindices(ka,:);
        end;
        cindices=cindices(cindices>0);
        if length(cindices)==3, % selection of secondary structure elements
            if isfield(model.structures{cindices(1)}(cindices(2)).residues{cindices(3)},'secondary_graphics')
                allgraphics=model.structures{cindices(1)}(cindices(2)).residues{cindices(3)}.secondary_graphics;
                if ~isempty(allgraphics),
                    for k=1:length(allgraphics),
                        graphics=allgraphics(k);
                        if ~isempty(graphics),
                            if ~isempty(graphics.objects),
                                for kk=1:length(graphics.objects),
                                    if isprop(graphics.objects(kk),'FaceColor'),
                                            set(graphics.objects(kk),'FaceColor',graphics.color(1,:));
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
        if ~isempty(cindices)
            [msg,graphics]=get_object(cindices,'graphics');
            if ~isempty(graphics) && ~isempty(graphics.objects),
                if graphics.mode==1,
                    for k=1:length(graphics.objects),
                        if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                            if isprop(graphics.objects(k),'Color'),
                                set(graphics.objects(k),'Color',graphics.color(1,:));
                            end;
                        end;
                    end;
                elseif graphics.mode>1,
                    for k=1:length(graphics.objects),
                        if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                            if isprop(graphics.objects(k),'FaceColor'),
                                set(graphics.objects(k),'FaceColor',graphics.color(1,:));
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

% Secondary structure elements

if isfield(model,'secondary_slected') && ~isempty(model.secondary_selected),
    for k=1:length(model.secondary_selected),
        graphics=model.secondary_selected{k};
        if ~isempty(graphics.objects),
            for kk=1:length(graphics.objects),
                if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                    set(graphics.objects(kk),'FaceColor',graphics.color(1,:));
                end;
            end;
        end;
    end;
end;
