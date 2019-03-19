function highlight_selection
% sets highlight color to current selection

global graph_settings
global model
global hMain

maxx=-1e6;
minx=1e6;
maxy=-1e6;
miny=1e6;
maxz=-1e6;
minz=1e6;

handles=guidata(hMain.figure);
set(handles.menu_edit_annotation,'Enable','off');
set(handles.menu_analysis_context,'Enable','off');
set(handles.menu_analysis_distance,'Enable','off');
set(handles.menu_analysis_angle,'Enable','off');
set(handles.menu_analysis_dihedral,'Enable','off');
set(handles.menu_analysis_alignment,'Enable','off');
set(handles.uipushtool_annotation,'Enable','off');
set(handles.uipushtool_context,'Enable','off');
set(handles.uipushtool_distance,'Enable','off');
set(handles.uipushtool_angle,'Enable','off');
set(handles.uipushtool_dihedral,'Enable','off');

indices=resolve_address('*');

if isempty(indices),
    return
end;

[m,n]=size(indices);

poi=0;
for k=1:m, % test for selection of at least two chains
    cindices1=indices(k,:);
    cindices1=cindices1(cindices1>0);
    if length(cindices1)==2,
        poi=poi+1;
    end;
end;
if poi>=2,
    set(handles.menu_analysis_alignment,'Enable','on');
end;

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
    for ka=0:ma, % ka=0 takes care of the selected object, k>0 of its descendants
        if ka>0,
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
                                    if isprop(graphics.objects(k),'FaceColor'),
                                        set(graphics.objects(k),'FaceColor',graph_settings.selected_color);
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
                                set(graphics.objects(k),'Color',graph_settings.selected_color);
                            end;
                            xx = get(graphics.objects(k),'XData');
                            yy = get(graphics.objects(k),'YData');
                            zz = get(graphics.objects(k),'ZData');
                            maxx=max([maxx,max(xx)]);
                            minx=min([minx,min(xx)]);
                            maxy=max([maxy,max(yy)]);
                            miny=min([miny,min(yy)]);
                            maxz=max([maxz,max(zz)]);
                            minz=min([minz,min(zz)]);
                        end;
                    end;
                elseif graphics.mode>1,
                    for k=1:length(graphics.objects),
                        if ishandle(graphics.objects(k)) && graphics.objects(k)~=0,
                            if isprop(graphics.objects(k),'FaceColor'),
                                set(graphics.objects(k),'FaceColor',graph_settings.selected_color);
                                xx = get(graphics.objects(k),'XData');
                                yy = get(graphics.objects(k),'YData');
                                zz = get(graphics.objects(k),'ZData');
                                maxx=max([maxx,max(xx)]);
                                minx=min([minx,min(xx)]);
                                maxy=max([maxy,max(yy)]);
                                miny=min([miny,min(yy)]);
                                maxz=max([maxz,max(zz)]);
                                minz=min([minz,min(zz)]);
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;


% Secondary structure elements

if isfield(model,'secondary_selected') && ~isempty(model.secondary_selected),
    disp(length(model.secondary_selected));
    for k=1:length(model.secondary_selected),
        graphics=model.secondary_selected{k};
        if ~isempty(graphics.objects),
            for kk=1:length(graphics.objects),
                if ishandle(graphics.objects(kk)) && graphics.objects(kk)~=0,
                    if isprop(graphics.objects(kk),'FaceColor'),
                        set(graphics.objects(kk),'FaceColor',graph_settings.selected_color);
                    end;
                    xx = get(graphics.objects(k),'XData');
                    yy = get(graphics.objects(k),'YData');
                    zz = get(graphics.objects(k),'ZData');
                    maxx=max([maxx,max(xx)]);
                    minx=min([minx,min(xx)]);
                    maxy=max([maxy,max(yy)]);
                    miny=min([miny,min(yy)]);
                    maxz=max([maxz,max(zz)]);
                    minz=min([minz,min(zz)]);
                end;
            end;
        end;
    end;
end;

if minx<1e5 && maxx>1e-5,
    model.extent_selection=[minx,miny,minz,maxx,maxy,maxz];
end;

switch m
    case 1,
        set(handles.menu_analysis_context,'Enable','on');
        set(handles.uipushtool_context,'Enable','on');
        set(handles.menu_edit_annotation,'Enable','on');
        set(handles.uipushtool_annotation,'Enable','on');
    case 2,
        set(handles.menu_analysis_distance,'Enable','on');
        set(handles.uipushtool_distance,'Enable','on');
    case 3,
        set(handles.menu_analysis_angle,'Enable','on');
        set(handles.uipushtool_angle,'Enable','on');
    case 4,
        set(handles.menu_analysis_dihedral,'Enable','on');
        set(handles.uipushtool_dihedral,'Enable','on');
end;        
