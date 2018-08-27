function cindices=selection_context(noninteractive)
% function cindices=selection_context
%
% Indices of all objects nearby the selected object
% user is queried for the context radius, unless input parameter
% noninteractive is present and nonzero
%
% only one object may be selected otherwise output is empty
%
% noninteractive    optional flag for supperssing user queries, 0 user is
%                   queried, 1 user is not queried, default 0
% cindices          index matrix of nearby objects
%
% G. Jeschke, 2009

global hMain
global geometry_settings

handles=guidata(hMain.figure);

cindices=[];

if nargin<1,
    noninteractive=0;
end;

indices=resolve_address('*');
indices=indices(indices>0);

if isempty(indices),
    add_msg_board('Nothing selected => No context.');
    return
end;

[m,n]=size(indices);

if m>1,
    add_msg_board('### Warning ### Context is defined only for selection of a single object.');
    return
end;

if n<2,
    add_msg_board('### Warning ### Context is defined only within the same structure.');
    return
end;

switch n
    case 2, % chain context
        radius=geometry_settings.chain_context;
        title='Chain context radius';
        quest='Radius in Å for chain context search:';
    case 3, % chain model context        
        radius=geometry_settings.chain_context;
        title='Chain context radius';
        quest='Radius in Å for chain context search:';
    case 4, % residue context        
        radius=geometry_settings.residue_context;
        title='Residue context radius';
        quest='Radius in Å for residue context search:';
    case 5, % atom context        
        radius=geometry_settings.atom_context;
        title='Atom context radius';
        quest='Radius in Å for atom context search:';
    case 6, % location context        
        radius=geometry_settings.atom_context;
        title='Location context radius';
        quest='Radius in Å for location context search:';
end;

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(quest,title,1,{sprintf('%6.1f',radius)},options);

newrad=str2double(answer);

if isnan(newrad)
    add_msg_board('### ERROR ### Provide numerical value for context radius.');
    return
end;

if newrad<0,
    add_msg_board('### ERROR ### Context radius must be positive.');
    return
end;

radius=newrad;

[cindices,info]=context(indices,radius);

n=length(info);
if ischar(info),
    add_msg_board(info);
else
    for k=1:n,
        add_msg_board(info{k});
    end;
end;



