function [snum,message]=copy_structure(snum0,id,Ca_coor,modnum,snum1)
% function [snum,message]=copy_structure(snum0,id,Ca_coor,modnum,snum1)
%
% Makes a copy of an MMM structure, keeping only basic information
% the identifier id, if present, is defined as a synonym for the structure
% coordinates can be transformed so that they conform to
% updated C_alpha coordinates, secondary structure definitions are NOT
% automatically redetermined, no annotations are made for the new
% structure, information on sites and metal binding is not copied, but 
% information on missing atoms and residues is
%
% snum0         number of the structure to be copied
% id            identifier of the new structure, should not be a possible
%               PDB identifier, suggestion for first character: +
%               should be four characters long, may not be an existing
%               structure tag
% Ca_coor       new C_alpha coordinates, residue coordinates are updated by
%               fitting a change of the local frame of the residue between
%               the initial and new C_alpha trace
%               this is optional, if argument is given, it is an array
%               [N,3], where N must coincide with the number of C_alpha
%               atoms in the copied structure, otherwise a message is set
%               and no copy is made
% modnum        optional model (dataset) number, if present, only this
%               model is replaced by the new C_alpha coordinates in
%               structure snum0 or added to this structure, no new
%               structure is made; call without modnum for the first model
%               (dataset)
% snum1         number of the target structure, if only a single model is
%               added, this parameter must be provided, if modnum is
%               present
% snum          number of the new (copied) structure in the model, is
%               empty, if no copy was made, is generated automatically, if
%               modnum is missing, otherwise it is equal to snum1
% message       message.error
%                   0   no error
%                   1   structure to be copied does not exist, no copy made
%                   2   id does already exist, no copy made
%                   3   mismatch of C_alpha trace length, no copy made
%                   4   modnum is present, but snum1 is missing
%               message.text  clear text error message
%
% G. Jeschke, 2010

global model

message.error=0;
message.text='';
snum=[]; % empty structure number, if it fails

if isempty(model) || ~isfield(model,'structures') || isempty(model.structures),
    message.error=1;
    message.text='No model defined yet or no structure in model.';
    add_msg_board('Error while copying structure');
    add_msg_board(message.text);
    return
elseif snum0>length(model.structures)
    message.error=1;
    message.text='Structure to be copied does not exist.';
    add_msg_board('Error while copying structure');
    add_msg_board(message.text);
    return    
else
    snum=length(model.structures)+1;
end;
if nargin>4,
    snum=snum1;
end;

exists=tag2id(id,model.structure_tags);
if ~isempty(exists),
    message.error=2;
    message.text='Structure identifier already exists.';
    add_msg_board('Error while copying structure');
    add_msg_board(message.text);
    snum=[];
    return    
end;    

if nargin<4,
    modnum=1;
else
    if nargin<5,
        message.error=4;
        message.text='Must specify target structure if model number is specified'; 
        add_msg_board('Error while copying structure');
        add_msg_board(message.text);
        snum=[];
        return    
    end;
end;

if nargin>2 && ~isempty(Ca_coor),
    [mCa,n]=size(Ca_coor);
    if ~isfield(model,'coarse') || length(model.coarse)<snum || isempty(model.coarse(snum).Ca_coor),
        [Ca_coor0,masses,rindices]=coarse_residues(mk_address(snum0));
    else
        Ca_coor0=model.coarse(snum0).Ca_coor;
        rindices=model.coarse(snum0).indices;
    end;
    [mCa0,n]=size(Ca_coor0);
    if mCa0~=mCa,
        message.error=3;
        message.text='Mismatch of C_alpha trace lengths for coordinate transformation.';
        add_msg_board('Error while copying structure');
        add_msg_board(message.text);
        snum=[];
        return    
    end;
end;


if nargin<4,
    model.structures{snum}=model.structures{snum0};
    % redistribute information
    model.chain_tags{snum}=model.chain_tags{snum0};
    model.chain_ids{snum}=model.chain_ids{snum0};
    model.info{snum}.title=sprintf('Structure copied from %s.',model.info{snum0}.idCode);
    if nargin>2,
        model.info{snum}.title=sprintf('%s Coordinates updated (transformed structure).',model.info{snum}.title);
    end;
    model.info{snum}.remarks=model.info{snum0}.remarks;
    model.info{snum}.class=model.info{snum0}.class;
    model.info{snum}.depDate='';
    model.info{snum}.idCode=id;
    model.info{snum}.center=model.info{snum0}.center;
    model.info{snum}.atoms=model.info{snum0}.atoms;
    model.info{snum}.residues=model.info{snum0}.residues;
    model.info{snum}.B_range=model.info{snum0}.B_range;
    model.info{snum}.SSbonds=model.info{snum0}.SSbonds;
    model.structure_tags=sprintf('%s%s:',model.structure_tags,id);
    model.structure_ids=[model.structure_ids snum];

    if isfield(model.info{snum0},'missing'),
        model.info{snum}.missing=model.info{snum0}.missing;
    end;
    if isfield(model.info{snum0},'incomplete'),
        model.info{snum}.incomplete=model.info{snum0}.incomplete;
    end;
    if isfield(model.info{snum0},'cryst'),
        model.info{snum}.cryst=model.info{snum0}.cryst;
    end;
    model.info{snum}.authors=sprintf('Structure based on PDB ID %s and modified with MMM.',model.info{snum0}.idCode);
    if isfield(model.info{snum0},'keywords'),
        model.info{snum}.keywords=model.info{snum0}.keywords;
    end;
else
    chains=length(model.structures{snum0});
    for cnum=1:chains,
        model.structures{snum}(cnum).isotopes=model.structures{snum0}(cnum).isotopes;
        model.structures{snum}(cnum).conn=model.structures{snum0}(cnum).conn;
        model.structures{snum}(cnum).maxconn=model.structures{snum0}(cnum).maxconn;
        if isfield(model.structures{snum}(cnum),'dbref'),
            model.structures{snum}(cnum).dbref=model.structures{snum0}(cnum).dbref;
        end;
        model.structures{snum}(cnum).atoms{modnum}=model.structures{snum0}(cnum).atoms{1};
        model.structures{snum}(cnum).residues{modnum}=model.structures{snum0}(cnum).residues{1};
        model.structures{snum}(cnum).xyz{modnum}=model.structures{snum0}(cnum).xyz{1};
        model.structures{snum}(cnum).isotopes=model.structures{snum0}(cnum).isotopes;
        model.structures{snum}(cnum).Bfactor{modnum}=model.structures{snum0}(cnum).Bfactor{1};
        model.structures{snum}(cnum).Btensor{modnum}=model.structures{snum0}(cnum).Btensor{1};
    end;
end;

if nargin<3 || isempty(Ca_coor), % pure copy, no coordinate update
    return
end;

[m,n]=size(rindices);
transmat2=zeros(4);
for k=1:m,
    cindices=rindices(k,:);
    cindices(1)=snum;
    cindices(3)=modnum;
    local_template=zeros(5,3);
    local_template_0=zeros(5,3);
    trans=Ca_coor(k,:)-Ca_coor0(k,:);
    transmat1=affine('translation',trans); % shift to centroid of structure to be fitted
    % make a local template to fit rotation and translation
    poi=0;
    for kk=-2:2,
        if k+kk>0 && k+kk<=m, % is addressed residue a network point?
            diff=rindices(k+kk,4)-rindices(k,4);
            if diff==kk, % is addressed residue part of a continuous segment?
                poi=poi+1;
                local_template(poi,:)=Ca_coor(k+kk,:);
                local_template_0(poi,:)=Ca_coor0(k+kk,:);
            end;
        end;
    end;
    if poi>=3, % found sufficient number of points to determine local rotation and translation
        % rotmat=opt_rot(local_template(1:poi,:),local_template_0(1:poi,:));
        [rms,coor2b,transmat]=rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
%         transmat2(1:3,1:3)=rotmat'; % rotate to new model frame
%         transmat2(4,4)=1; 
        set_residue(cindices,'transform',transmat);
    else
        set_residue(cindices,'transform',transmat1);
    end;
end;
