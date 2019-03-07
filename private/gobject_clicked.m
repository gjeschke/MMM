function gobject_clicked(cb,eventdata)
% gobject_clicked(cb,eventdata)
%
% Function executed when a user clicks on a cartoon graphics object

global model
global hMain

mouseclick=get(gcf,'SelectionType');

if ~strcmpi(mouseclick,'extend') && ~strcmpi(mouseclick,'alt')

    % Deselect old selection
    if isfield(model,'selected') && ~isempty(model.selected) % was anything selected
        hMain=cmd(hMain,'unselect *');
        model.secondary_selected={};
    end
end

if strcmpi(mouseclick,'open')
    add_msg_board('--- All deselected ---');
    return
end

poi=find(model.graphics_objects(1:model.graphics_lookup_pointer)==cb);
indices=model.graphics_lookup(poi,2:end);

if indices(1)<0 % special handling for indexed isosurfaces
    if isfield(model,'surfaces')
        dg = model.surfaces(-indices(1));
        invcol = [1,1,1] - dg.color;
        set(dg.gobjects,'FaceColor',invcol);
        pause(0.25);
        set(dg.gobjects,'FaceColor',dg.color);
        if strcmpi(mouseclick,'alt')
            camlookat(dg.gobjects);
            add_msg_board(sprintf('Zoomed in on isosurface %s.',dg.tag));
        else
            add_msg_board(sprintf('Click on isosurface %s with center of gravity (%4.1f, %4.1f, %4.1f) Å and highest probability (%4.1f, %4.1f, %4.1f) Å.',dg.tag,dg.cog,dg.mlp));
        end
    else
        add_msg_board('BUG: Clicked surface object was not registered.');
    end
    return
end


loctags0=':A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:';
loc_tag='';
if indices(6)>0 && indices(5)==0 % special handling for label graphics
    rotamer=indices(6);
    locid=mod(rotamer,26);
    if locid==0, locid=26; end
    locmid=floor((rotamer-1)/26);
    if locmid>0
        loc_tag=strcat(id2tag(locmid,loctags0),id2tag(locid,loctags0));
    else
        loc_tag=id2tag(locid,loctags0);
    end
    indices(6)=0;
end

if indices(4)<0 % index extension for secondary structure elements
    allgraphics=model.structures{indices(1)}(indices(2)).residues{indices(3)}.secondary_graphics;
    if ~isempty(allgraphics)
        graphics=allgraphics(-indices(4));
        if ~isempty(graphics.objects)
            if ~isempty(graphics.range)
                address=sprintf('%s%i',mk_address(indices(1:3)),graphics.range(1));
                if graphics.range(2)>graphics.range(1)
                    address=sprintf('%s-%i',address,graphics.range(2));
                end
                if strcmpi(mouseclick,'alt')
                    command=sprintf('unselect %s',address);
                    if ~isempty(model.secondary_selected)
                        poi=0;
                        for k=1:length(model.secondary_selected)
                            if graphics.objects~=model.secondary_selected{k}.objects,
                                poi=poi+1;
                                model.secondary_selected{poi}=model.secondary_selected{k};
                            end
                        end
                    end
                    model.secondary_selected=model.secondary_seleceted{1:poi};
                else
                    command=sprintf('select %s',address);
                    if isempty(model.secondary_selected)
                        model.secondary_selected{1}=graphics;
                        model.secondary_indices{1}=indices(1:4);
                    else
                        n=length(model.secondary_selected);
                        model.secondary_selected{n+1}=graphics;
                        model.secondary_indices{n+1}=indices(1:4);
                    end                
                end
                hMain=cmd(hMain,command);
            end
        end
    end
else
    indices=indices(indices>0);
    if ~isempty(indices) % store selection
        address=mk_address(indices);
        if length(indices)==4 && ~isempty(loc_tag)
            address=sprintf('%s.:%s',address,loc_tag);
        end
        if strcmpi(mouseclick,'alt')
            command=sprintf('unselect %s',address);
        else
            command=sprintf('select %s',address);
        end
        hMain=cmd(hMain,command);
    end
end

% % Display of current selection
% 
% full_info={'--- Selected ---'};
% poi=1;
% if ~strcmpi(mouseclick,'open') && model.selections>=1,
%     k0=1;
%     if model.selections>10,
%         full_info{2}='...';
%         poi=2;
%         k0=model.selections-9;
%     end;
%     for k=k0:model.selections,
%         [message,prev_info]=get_object(model.selected{k},'info');
%         poi=poi+1;
%         full_info{poi}=prev_info{1};
%     end;
% end;
% % if ~strcmpi(mouseclick,'alt')
% %     full_info{poi+1}='is:';
% %     for k=2:length(info_text),
% %         full_info{poi+k}=info_text{k};
% %     end;
% % end;
% set(hMain.text_message_board,'String',full_info);
