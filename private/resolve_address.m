function [indices,message]=resolve_address(address)
% function [indices,message]=resolve_address(address)
%
% --- Object address resolution for MMM ---
% Tests whether 'address' is a syntactically correct address of the model
% contained in the global variable 'model' and whether it corresponds to an
% object or several objects within this model
% if so, the indices into the Matlab array of structures are returned for
% all matching objects, if not, an empty index array is returned
% any errors are reported in message
% the default structure and chain are defined in global variable 'defaults'
% if defaults is empty, the first chain of the first model is default
%
% resolve_address does not handle coarse-graining level identifiers
% these are ignored, if present
%
% global variables model and defaults must exist, otherwise error is thrown
%
% address   address of the form
%           [structure_tag](chain_tag){coordinate_set}residue#.atom_tag:location_tag
%           see Manual for list addressing and defaults
% indices   array [m,n] of index vectors of the m objects selected by
%           address, n corresponds to the hierarchy level:
%           n=1  whole structure
%           n=2  chain
%           n=3  coordinate set within chain
%           n=4  residue
%           n=5  atom
%           n=6  alternate location of atom
% message   error message structure with fields:
%           .error  error code: 0 no error, 1: address syntax error, 2: no
%                               object with this address exists
%           .text   text error message

global model
global chemistry

if ~isfield(model,'structures')
    message.error=13;
    message.text='No structure.';
    indices=[];
    return
end

% Definition of type synonyms
tsyn{1}.short='"Xaa"';
tsyn{1}.long='"Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Asx","Glx","Xaa"';
tsyn{2}.short='"Asx"';
tsyn{2}.long='"Asn","Asp"';
tsyn{3}.short='"Glx"';
tsyn{3}.long='"Gln","Glu"';
tsyn{4}.short='"R"';
tsyn{4}.long='"G","A"';
tsyn{5}.short='"Y"';
tsyn{5}.long='"T","C"';
tsyn{6}.short='"K"';
tsyn{6}.long='"G","T"';
tsyn{7}.short='"M"';
tsyn{7}.long='"A","C"';
tsyn{8}.short='"S"';
tsyn{8}.long='"G","C"';
tsyn{9}.short='"W"';
tsyn{9}.long='"A","T"';
tsyn{10}.short='"B"';
tsyn{10}.long='"G","T","C"';
tsyn{11}.short='"D"';
tsyn{11}.long='"G","A","T"';
tsyn{12}.short='"H"';
tsyn{12}.long='"A","T","C"';
tsyn{13}.short='"V"';
tsyn{13}.long='"G","A","C"';
tsyn{14}.short='"N"';
tsyn{14}.long='"G","A","T","C"';
tsyn{15}.short='"X"';
tsyn{15}.long='"G","A","T","C"';
tsyn{16}.short='"U"';
tsyn{16}.long='"T"';

type_synonyms=':';
for k=1:length(tsyn),
    type_synonyms=[type_synonyms tsyn{k}.short ':'];
end;

% Definition of element synonyms
esyn{1}.short='"X"';
esyn{1}.long='"F","Cl","Br","I"';
element_synonyms=':';
for k=1:length(esyn),
    element_synonyms=[element_synonyms esyn{k}.short ':'];
end;


max_objects=20000; % only for preallocation, selection of more objects 
                   % slows down address resolution, but it still works
m=0; % initialize object number
n=0; % initialize hierarchy level

indices=zeros(max_objects,6); % initialize index array (pre-allocation) 

message.error=0; % initialize message with no error code
message.text=''; % and no clear text message

if isempty(address) % empty address selects no object
    message=no_object;
    indices=[];
    return
else
    address=strtrim(address);
end;

% Handling for special addresses
if strcmp(address,'*') % returns adresses of all selected objects
    if ~isfield(model,'selections')
        indices=[];
        message.error=3;
        message.text='Nothing selected';
        return
    end;
    if model.selections<1,
        indices=[];
        message.error=3;
        message.text='Nothing selected';
        return
    end;
    for k=1:model.selections,
        indices(k,1:length(model.selected{k}))=model.selected{k};
    end;
    indices=indices(1:model.selections,:);
    return
end;

if strcmp(address,'!'),
    address='[*]';
end;

if strcmp(address,'#'),
    address=sprintf('[*](%s)',model.current_chain);
end;


% Remove non-address parts of input string
a_stop=strfind(address,';'); % check for end of address string
if a_stop==1 % return syntax error, if address is just a semicolon
    message=syntax_error;
    indices=[];
    return
end;
if ~isempty(a_stop),
    address=address(1:a_stop-1); % cut off semicolon and everything thereafter
end;

% Preprocess domain addressing

domain_request=false;
subdomain_request=false;
sub_tag='';
dom_tag='';
dom_limit=strfind(address,'|');
if ~isempty(dom_limit)
    domain_request=true;
    % keyboard
    if dom_limit(1)==1, % extension for domain addressing in current structure
        snum=model.current_structure;
        stag=id2tag(snum,model.structure_tags,model.structure_ids);
        address0=sprintf('[%s]%s',stag,address);
    else
        address0=address;
    end;
    if length(dom_limit)~=2,
        message=syntax_error;
        indices=[];
        return
    else
        if dom_limit(2)-dom_limit(1)<2,
            message=no_object;
            indices=[];
            return
        end;
        dom_tag=address(dom_limit(1)+1:dom_limit(2)-1);
        sec_start=strfind(dom_tag,'<');
        if ~isempty(sec_start),
            sec_stop=strfind(dom_tag,'>');
            if isempty(sec_stop) || sec_stop-sec_start<2,
                message=syntax_error;
                indices=[];
                return
            end;
            sub_tag=dom_tag(sec_start+1:sec_stop-1);
            subdomain_request=true;
            dom_tag=dom_tag(1:sec_start-1);
        end;
    end;
    address=address0;
end;

% Pre-process secondary structure string, if any

secondary_request=0;

t_start=strfind(address,'<');
t_mid=strfind(address,'>');
t_stop=strfind(address,'.');

if ~isempty(t_start)  % a secondary structure element request exists
    secondary_request=1;
    if isempty(t_mid) || t_mid<t_start,
        message=syntax_error;
        indices=[];
        return
    end;
    if t_mid-t_start<4, % no place for type and name
        message=no_object;
        indices=[];
        return
    end;
    if isempty(t_stop), % no dot separator between type and name
        message=syntax_error;
        indices=[];
        return
    end;
    if t_stop(1)-t_start~=2
        message=syntax_error;
        indices=[];
        return
    else
        sec_type=address(t_start+1);
    end;
    if t_mid-t_stop(1)<2,
        message=no_object;
        indices=[];
        return
    end;
    sec_tag=address(t_stop(1)+1:t_mid-1);
    if length(t_stop)<2,
        t_stop=length(address)+1;
    else
        t_stop=t_stop(2);
    end;
    if t_stop-t_mid<2
        residue_tag=':'; % no residue specification means all residues
    else
        residue_tag=address(t_mid+1:t_stop-1);
    end;
    if t_start>1,
        header=address(1:t_start-1);
    else
        header='';
    end;
    if t_stop<length(address),
        trailer=address(t_stop:end);
    else
        trailer='';
    end;
end;


% Look for and resolve structure tag
t_start=strfind(address,'['); 
t_stop=strfind(address,']');
if isempty(t_stop) || isempty(t_start) || t_stop-t_start==1 % no structure tag or empty tag []
    if isfield(model,'current_structure') % find out how default structure is defined and select it
        struct_ind=model.current_structure;
    else
        struct_ind=1;
    end;
    n=1;
else
    if length(t_start)>1 || length(t_stop)>1 || t_stop<t_start % all these are syntax errors
        message=syntax_error;
        indices=[];
        return
    end;
    structure_tag=address(t_start+1:t_stop-1); % cut out structure tag
    % Special handling for explicit reference to current structure
    if strcmpi(structure_tag,'*'),
        struct_ind=model.current_structure;
    else
        if isfield(model,'structure_tags') && isfield(model,'structure_ids') % check whether synonyms are defined
            tags=model.structure_tags;
            ids=model.structure_ids;
        else % generate tags string and code list
            tags=':';
            for k=1:length(model.structures),
                tags=[tags sprintf('%i:',k)];
            end;
            ids=1:length(model.structures);
        end;
        struct_ind=taglist2indices(upper(structure_tag),upper(tags),ids);
    end;
    if isempty(struct_ind) || sum(isnan(struct_ind)) % object does not exist
        message=no_object;
        indices=[];
        return;
    end;
    % Process domain request, if any
    if domain_request,
        indices=[];
        for sk=1:length(struct_ind),
            sind=struct_ind(sk);
            dobjects=[];
            if isfield(model,'domains')
                if sk<=length(model.domains)
                    for dk=1:length(model.domains{sind})
                        if strcmp(dom_tag,model.domains{sind}(dk).name),
                            dobjects=model.domains{sind}(dk).objects;
                            if subdomain_request
                                sub_ind=[];
                                if isfield(model.domains{sind}(dk),'subdomains'),
                                    for ssk=1:length(model.domains{sind}(dk).subdomains)
                                        if strcmp(sub_tag,model.domains{sind}(dk).subdomains{ssk}.name)
                                            sub_ind=model.domains{sind}(dk).subdomains{ssk}.indices;
                                        end;
                                    end;
                                end;
                                dobjects=dobjects(sub_ind',:);
                            end;
                        end;
                    end;
                end;
            end;
            [mm,nn]=size(dobjects);
            if mm>0 % there are domain objects for this structure
                if isempty(indices) % no other domain was yet stored
                    indices=dobjects;
                else % there are already objects stored and the hierarchy level must match
                    [mmm,nnn]=size(indices);
                    if nn~=nnn,
                        message.error=3;
                        message.text='Domain size mismatch';
                        indices=[];
                    else
                        indices=[indices;dobjects];
                    end;
                end;
            end;
        end;
        return % the domain request is completely processed here
    end;
    n=1; % hierarchy level of objects is at least structure
end;
indices=struct_ind';
m=length(struct_ind);
indices0=indices;
[m0,n0]=size(indices0);

% Look for and resolve chain tag
t_start=strfind(address,'('); 
t_stop=strfind(address,')');
if isempty(t_stop) || isempty(t_start) || t_stop-t_start==1 % no chain tag or empty tag ()
    if isfield(model,'chain_tags') && numel(struct_ind)==1 && struct_ind==model.current_structure, % find out how default chain is defined and select it
        chain_ind=tag2id(model.current_chain,model.chain_tags{struct_ind});
    else
        chain_ind=1;
    end;
    % Expanded index array for just the default chain, needed if address
    % refers to lower levels than chains
    indices=indices0;
    indices0=zeros(m0,2);
    for ipoi=1:m0, 
        indices0(ipoi,1)=indices(ipoi,1);
        indices0(ipoi,2)=chain_ind;
    end;
    n0=2;
    % Empty parentheses () are interpreted as reference to default chain
    if t_stop-t_start==1,
        indices=indices0;
        m=m0;
        n=2;
    end;
else
    if length(t_start)>1 || length(t_stop)>1 || t_stop<t_start % all these are syntax errors
        message=syntax_error;
        indices=[];
        return
    end;
    chain_tag=upper(address(t_start+1:t_stop-1)); % cut out chain tag and convert to upper case
    % Special handling for explicit reference to current chain
    if strcmpi(chain_tag,'*')
        chain_tag=model.current_chain;
    end;
    m=0;
    [m0,n0]=size(indices0);
    for ipoi=1:m0, % loop over all existing index vectors
        snum=indices0(ipoi,1);
        if isfield(model,'chain_tags')
            tags=model.chain_tags{snum};            
            if isfield(model,'chain_ids') % check whether synonyms are defined
                if ~isempty(model.chain_ids{snum})
                    ids=model.chain_ids{snum};
                end;
            else
                ids=1:length(model.structures{snum});
            end;
        else % generate tags string and code list
            tags=':';
            bas=double('A')-1;
            for k=1:length(model.structures{snum}),
                tags=[tags char(bas+k)];
            end;
            ids=1:length(model.structures{snum});
        end;
        chain_ind=taglist2indices(chain_tag,upper(tags),ids);
        if sum(isnan(chain_ind)) % something was addressed that does not exist
            message=syntax_error;
            indices=[];
            return
        end;
        if ~isempty(chain_ind) % check whether in the current structure chains are addressed
            for cpoi=1:length(chain_ind) % store all chains
                m=m+1;
                indices(m,1)=indices0(ipoi,1);
                indices(m,2)=chain_ind(cpoi);
            end;
        end;
    end;
    if m==0,
        message=no_object;
        indices=[];
        return
    end;
    indices=indices(1:m,1:2);
    n=2; % hierarchy level of objects is at least chain
    indices0=indices;
    [m0,n0]=size(indices0);
end;

% Look for and resolve coordinate set tag
t_start=strfind(address,'{'); 
t_stop=strfind(address,'}');
if isempty(t_stop) || isempty(t_start) || t_stop-t_start==1 % no coordinate set tag or empty tag {}
    coorset_ind=1;
    % Expanded index array for just the default set, needed if address
    % refers to lower levels than chains
    indices00=indices0;
    indices0=zeros(m0,3);
    for ipoi=1:m0, 
        indices0(ipoi,1:2)=indices00(ipoi,1:2);
        indices0(ipoi,3)=coorset_ind;
    end;
    n0=3;
    % Empty braces {} are interpreted as reference to first coordinate set
    if t_stop-t_start==1,
        indices=indices0;
        m=m0;
        n=3;
    end;
else
    indices=indices0;
    n=n0;
    m=m0;
    if length(t_start)>1 || length(t_stop)>1 || t_stop<t_start % all these are syntax errors
        message=syntax_error;
        indices=[];
        return
    end;
    coorset_tag=address(t_start+1:t_stop-1); % cut out coordinate set tag
    m=0;
    [m0,n0]=size(indices0);
    for ipoi=1:m0, % loop over all existing index vectors
        snum=indices0(ipoi,1);
        cnum=indices0(ipoi,2);
        sets=length(model.structures{snum}(cnum).residues);
        ids=1:sets;
        tags=':';
        for setpoi=1:sets,
            tags=sprintf('%s%i:',tags,setpoi);
        end;
        set_ind=taglist2indices(coorset_tag,tags,ids);
        if sum(isnan(set_ind)) % something was addressed that does not exist
            message=syntax_error;
            indices=[];
            return
        end;
        if ~isempty(set_ind) % check whether in the current structure chains are addressed
            for setpoi=1:length(set_ind) % store all sets
                m=m+1;
                indices(m,1:2)=indices0(ipoi,1:2);
                indices(m,3)=set_ind(setpoi);
            end;
        end;
    end;
    if m==0,
        message=no_object;
        indices=[];
        return
    end;
    indices=indices(1:m,1:3);
    n=3; % hierarchy level of objects is at least coordinate set
    indices0=indices;
    [m0,n0]=size(indices0);
end;

% Look for and resolve residue tag
if ~secondary_request, % if a secondary structure element is requested, the residue tag is already made
    head=strfind(address,']'); 
    head2=strfind(address,')');
    head3=strfind(address,'}');
    if isempty(head), head=0; end;
    if isempty(head2), head2=0; end;
    if isempty(head3), head3=0; end;
    if length(head)>1 || length(head2)>1 || length(head3)>1 % all these are syntax errors
        message=syntax_error;
        indices=[];
        return
    end;
    t_start=max([head,head2,head3]);
    if t_start==length(address), % address is fully resolved at structure or chain or coordinate set level
        return;
    end;
    t_stop=strfind(address,'.');
    if isempty(t_stop), 
        residue_tag=address(t_start+1:length(address)); % cut out residue tag
    elseif t_stop-t_start==1 % empty residue tag means all residues
        residue_tag=':';
    elseif t_stop<=t_start % dot comes too early
        message=syntax_error;
        indices=[];
        return    
    else
        residue_tag=address(t_start+1:t_stop-1);
    end;
end;

% if we get here, there is a request to address residues

indices=indices0;
n=n0;

[m0,n0]=size(indices0);
m=0;
residue_tag_0=residue_tag;
for ipoi=1:m0, % loop over all existing index vectors
    snum=indices0(ipoi,1); % structure index
    cnum=indices0(ipoi,2); % chain index
    mnum=indices0(ipoi,3); % coordinate set index
    residue_ind=[];
    % Handle secondary element request, if any
    if secondary_request % handle secondary structure request
        range=[];
        sec_tags=':';
        switch upper(sec_type)
            case 'H'
                if isfield(model.structures{snum}(cnum),'helix_defs'),
                    secnum=length(model.structures{snum}(cnum).helix_defs);
                    sec_ids=1:secnum;
                    for sk=1:secnum,
                        sec_tags=[sec_tags model.structures{snum}(cnum).helix_defs{sk}.name ':'];
                    end;
                    sec_ind=tag2id(sec_tag,sec_tags,sec_ids);
                    if ~isempty(sec_ind)
                        range=model.structures{snum}(cnum).helix_defs{sec_ind(length(sec_ind))}.range;
                    end;
                end;
            case 'E'
                if isfield(model.structures{snum}(cnum),'sheet_defs'),
                    secnum=length(model.structures{snum}(cnum).sheet_defs);
                    sec_ids=1:secnum;
                    for sk=1:secnum,
                        sec_tags=[sec_tags model.structures{snum}(cnum).sheet_defs{sk}.name ':'];
                    end;
                    sec_ind=tag2id(sec_tag,sec_tags,sec_ids);
                    if ~isempty(sec_ind)
                        range=model.structures{snum}(cnum).sheet_defs{sec_ind(length(sec_ind))}.range;
                    end;
                end;
            case 'L'
                if isfield(model.structures{snum}(cnum),'loop_defs'),
                    secnum=length(model.structures{snum}(cnum).loop_defs);
                    sec_ids=1:secnum;
                    for sk=1:secnum,
                        sec_tags=[sec_tags model.structures{snum}(cnum).loop_defs{sk}.name ':'];
                    end;
                    sec_ind=tag2id(sec_tag,sec_tags,sec_ids);
                    if ~isempty(sec_ind)
                        range=model.structures{snum}(cnum).loop_defs{sec_ind(length(sec_ind))}.range;
                    end;
                end;
        end;
        if isempty(range),
            residue_tag='0';
        else
            if strcmp(residue_tag_0,':')
                residue_tag=sprintf('%i-%i',range(1),range(2));
            elseif isempty(strfind(residue_tag,'"')) % apply offset to residue numbers, if number request
                sec_length=range(2)-range(1)+1;
                num_list=textscan(residue_tag,'%s','Delimiter',',-');
                items=length(num_list{1});
                for sk=1:items,
                    rel_num=str2double(num_list{1}{sk});
                    if rel_num>sec_length, rel_num=sec_length; end; % limit to length of secondary structure element
                    num_list{1}{sk}=sprintf('%i',rel_num+range(1)-1);
                end;
                res_tag_0=residue_tag;
                residue_tag=num_list{1}{1};
                sk=1;
                for kk=1:length(res_tag_0),
                    if strcmp(res_tag_0(kk),',')
                        sk=sk+1;
                        residue_tag=[residue_tag ',' num_list{1}{sk}];
                    end;
                    if strcmp(res_tag_0(kk),'-')
                        sk=sk+1;
                        residue_tag=[residue_tag '-' num_list{1}{sk}];
                    end;
                end;
            end;
        end;
    end;
    % check whether residues are selected by type
    quote=findstr(residue_tag,'"');
    if isempty(quote) % selection by residue number
        % generate tags string and code list
        tags=model.structures{snum}(cnum).residues{mnum}.residue_tags;
        ids=1:length(model.structures{snum}(cnum).residues{mnum}.info);
        residue_ind=taglist2indices(residue_tag,tags,ids);
        if sum(isnan(residue_ind)) % something was addressed that does not exist
            message=syntax_error;
            indices=[];
            return
        end;
    else % selection by residue type
        residue_ind=[];
        % treatment for type name synonyms
        nonsense=textscan(residue_tag,'%s','Delimiter',',');
        tagcells=nonsense{1};
        residue_tag='';
        for k=1:length(tagcells) 
            tag=strtrim(char(tagcells(k))); % reformat to char and remove leading/trailing spaces
            syn_id=tag2id(tag,type_synonyms);
            if ~isempty(syn_id),
                tag=tsyn{syn_id(1)}.long; % defensive programming, if synonym was accidentally defined twice
            end;
            residue_tag=[residue_tag tag ','];
        end;
        residue_tag=residue_tag(1:end-1);
        % resolve comma-separated list 
        nonsense=textscan(residue_tag,'%s','Delimiter',',');
        tagcells=nonsense{1};
        for k=1:length(tagcells) 
            tag=strtrim(char(tagcells(k))); % reformat to char and remove leading/trailing spaces
            quote=findstr(tag,'"');
            if length(quote)~=2 % non-matching quotes or more than two quotes between commas
                message=syntax_error;
                indices=[];
                return
            end;
            if quote(2)-quote(1)<2, % quotes follow each other without residue name
                message=syntax_error;
                indices=[];
                return
            end;
            tag=tag(quote(1)+1:quote(2)-1);
            for rnum=1:length(model.structures{snum}(cnum).residues{mnum}.info),
                within_range=1;
                if secondary_request,
                    if isempty(range),
                        within_range=0;
                    else
                        curr_num=model.structures{snum}(cnum).residues{mnum}.info(rnum).number;
                        if curr_num<range(1) || curr_num>range(2)
                            within_range=0;
                        end;
                    end;
                end;
                if within_range && strcmpi(tag,model.structures{snum}(cnum).residues{mnum}.info(rnum).name) % match?
                    if ~sum(find(residue_ind==rnum)), % protect against listing the same residue twice
                        residue_ind=[residue_ind rnum];
                    end;
                end;
            end;
        end;
    end;
    if ~isempty(residue_ind) % check whether in the current structure, chain, coordinate set residues are addressed
        for rpoi=1:length(residue_ind) % store all residues
            m=m+1;
            indices(m,1:3)=indices0(ipoi,1:3);
            indices(m,4)=residue_ind(rpoi);
        end;
    end;
end;
if m==0,
    message=no_object;
    indices=[];
    return;
end;
indices=indices(1:m,1:4);
n=4; % hierarchy level of objects is at least residue
indices0=indices;
[m0,n0]=size(indices0);
if secondary_request, % remove secondary request and substitute by placeholder residue information
     address=sprintf('%s1%s',header,trailer); % necessary to remove dot, or further address analysis will go astray
end;
% Look for and resolve atom tag

loc_tag=''; % define empty location tag, as atom tag resolution also provides the location tag
t_start=strfind(address,'.');
if isempty(t_start), % no atom addressing 
    return
elseif t_start==length(address) % dot only addresses all atoms
    all_atoms=1;
else
    all_atoms=0;
end;
atom_tag=address(t_start:end);
t_stop=strfind(atom_tag,':'); % check, whether alternate locations are addressed
if ~isempty(t_stop), % cut off alternate location tag
    loc_tag=atom_tag(t_stop:end);
    atom_tag=atom_tag(1:t_stop-1);
    if length(atom_tag)==1, % only dot
        all_atoms=1;
    end;
end; 
if all_atoms,
    atom_tag=':';
else
    atom_tag=strtrim(atom_tag(2:end)); % cut off leading dot and any spaces
end;

% if we get here, there is a request to address atoms

indices=indices0;
n=n0;

[m0,n0]=size(indices0);
m=0;
for ipoi=1:m0, % loop over all existing index vectors
    snum=indices0(ipoi,1); % structure index
    cnum=indices0(ipoi,2); % chain index
    mnum=indices0(ipoi,3); % coordinate set index
    rnum=indices0(ipoi,4); % residue index
    % check whether atoms are selected by element
    quote=findstr(atom_tag,'"');
    if isempty(quote) % selection by atom type
        % generate tags string and code list
        tags=model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_tags;
        ids=1:length(model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers);
        atom_ind=taglist2indices(atom_tag,tags,ids);
        if sum(isnan(atom_ind)) % something was addressed that does not exist
            message=syntax_error;
            indices=[];
            return
        end;
    else % selection by element type
        atom_ind=[];
        % treatment for type name synonyms
        nonsense=textscan(atom_tag,'%s','Delimiter',',');
        tagcells=nonsense{1};
        atom_tag='';
        for k=1:length(tagcells) 
            tag=strtrim(char(tagcells(k))); % reformat to char and remove leading/trailing spaces
            syn_id=tag2id(tag,element_synonyms);
            if ~isempty(syn_id),
                tag=esyn{syn_id(1)}.long; % defensive programming, if synonym was accidentally defined twice
            end;
            atom_tag=[atom_tag tag ','];
        end;
        atom_tag=atom_tag(1:end-1);
        % resolve comma-separated list 
        nonsense=textscan(atom_tag,'%s','Delimiter',',');
        tagcells=nonsense{1};
        for k=1:length(tagcells) 
            tag=strtrim(char(tagcells(k))); % reformat to char and remove leading/trailing spaces
            quote=strfind(tag,'"');
            if length(quote)~=2 % non-matching quotes or more than two quotes between commas
                message=syntax_error;
                indices=[];
                return
            end;
            if quote(2)-quote(1)<2, % quotes follow each other without residue name
                message=syntax_error;
                indices=[];
                return
            end;
            tag=tag(quote(1)+1:quote(2)-1);
            elm_id=tag2id(upper(tag),upper(chemistry.element_tags));
            for anum=1:length(model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers),
                if elm_id==model.structures{snum}(cnum).residues{mnum}.info(rnum).elements(anum) % match?
                    if ~sum(find(atom_ind==anum)), % protect against listing the same residue twice
                        atom_ind=[atom_ind anum];
                    end;
                end;
            end;
        end;
    end;
    if ~isempty(atom_ind) % check whether in the current object atoms are addressed
        for apoi=1:length(atom_ind) % store all atoms
            m=m+1;
            indices(m,1:4)=indices0(ipoi,1:4);
            indices(m,5)=atom_ind(apoi);
        end;
    end;
end;
if m==0,
    message=no_object;
    indices=[];
    return;
end;
indices=indices(1:m,1:5);
n=5; % hierarchy level of objects is at least atom
indices0=indices;
[m0,n0]=size(indices0);

if isempty(loc_tag)
    return
end;

if length(loc_tag)==1, % select all locations
    loc_tag=':';
else
    loc_tag=loc_tag(2:end);
end;

enforce_op=strfind(loc_tag,'!'); % check whether alternate location addressing is enforced
if isempty(enforce_op) 
    enforce=0;
else
    enforce=1;
    if length(loc_tag)==1 % only enforce operator '!'
        loc_tag=':'; % select all alternative locations
    else
        loc_tag=loc_tag(1:enforce_op-1); % cut off enforce operator (and everything thereafter)
    end;
end;
    
% if we get here, there is a request to address locations

indices=indices0;
n=n0;

[m0,n0]=size(indices0);
m=0;
for ipoi=1:m0, % loop over all existing index vectors
    snum=indices0(ipoi,1); % structure index
    cnum=indices0(ipoi,2); % chain index
    mnum=indices0(ipoi,3); % coordinate set index
    rnum=indices0(ipoi,4); % residue index
    anum=indices0(ipoi,5);
    % generate tags string and code list
    if ~isfield(model.structures{snum}(cnum).residues{mnum}.info(rnum),'location_tags') % no location identifier for this residue
        ma=1; na=1; % in this case there is only one location per atom
        if enforce,
            loc_ind=[]; % but enforced location selection => nothing selected
        else
            loc_ind=1; % selection not enforced, default index 1
        end;
    else        
        [ma,na]=size(model.structures{snum}(cnum).residues{mnum}.info(rnum).atom_numbers{anum});
        if enforce && ma==1, % only alternate location s requested, but there are none
            loc_ind=[];
        else
            tags=model.structures{snum}(cnum).residues{mnum}.info(rnum).location_tags;
            loc_num=length(strfind(tags,':'))-1;
            if loc_num<1, loc_num=1; end;
            ids=1:loc_num;
            loc_ind=taglist2indices(loc_tag,tags,ids);
            if sum(isnan(loc_ind)) % something was addressed that does not exist
                message=syntax_error;
                indices=[];
                return
            end;
            if isempty(loc_ind) && ~enforce, % no match, but alternate location is not to be enforced
                loc_ind=1;
            end;
        end;
    end;
    if ~isempty(loc_ind) % check whether in the current object atoms are addressed
        for lpoi=1:length(loc_ind) % store all locations
            if loc_ind(lpoi)<=ma, % address only existing locations
                m=m+1;
                indices(m,1:5)=indices0(ipoi,1:5);
                indices(m,6)=loc_ind(lpoi);
            end;
        end;
    end;
end;
if m==0,
    message=no_object;
    indices=[];
    return;
end;
indices=indices(1:m,1:6);
n=6; % hierarchy level of objects is location
indices0=indices;
[m0,n0]=size(indices0);


% The following is only defensive programming, should never be necessary

if n==0 || m==0, 
    message=no_object;
    indices=[];
    return
end;

indices=indices(1:m,1:n);

% --- Private functions ---

function message=syntax_error
% returns the error message structure for a syntax error

message.error=1;
message.text='Invalid address.';

function message=no_object
% returns the error message structure if no object is addressed

message.error=2;
message.text='Object does not exist.';

function indices=taglist2indices(taglist,tags,codelist)
% Given a tag string, possibly with list addressing, this function returns
% a list of indices, mapping of tags to indices is code by string tags and
% identifier list codelist, using function tag2id
% tag errors lead to an empty index list

indices=[]; % initialize empty return value

taglist=strtrim(taglist); % remove any leading or trailing spaces

% processing for address all operator ':'
if strfind(taglist,':') 
    if length(taglist)>1, % but for this the tag list must be exactly ':'
        return % no selection for ':' in combination with anything else
    else
        for k=1:length(codelist), % everything, but an object with a tag synonym should not be indexed twice
            if ~sum(indices==codelist(k))
                indices=[indices codelist(k)];
            end;
        end;
        return
    end;
end;

if isempty(taglist),
    message=no_object;
    indices=nan;
    return;
end;

% processing for comma-separated lists
nonsense=textscan(taglist,'%s','Delimiter',',');
tagcells=nonsense{1};
for k=1:length(tagcells) 
    tag=strtrim(char(tagcells(k))); % reformat to char and remove leading/trailing spaces
    if strfind(tag,'-') % check for range of residues (or structures)
        nonsense=textscan(tag,'%s','Delimiter','-');
        two_tags=nonsense{1};
        if length(two_tags)~=2,
            indices=[];
            return % wrong format of range (start or end missing or more than one '-')
        end;
        start_num=round(str2double(two_tags(1)));
        end_num=round(str2double(two_tags(2)));
        if isnan(start_num) || isnan(end_num) || start_num>end_num % wrong format of range, no numbers or end<start
            indices=[];
            return
        end;
        for kk=start_num:end_num,
            mtag=sprintf('%i',kk);
            id=tag2id(mtag,tags,codelist);
            if ~isempty(id) && ~sum(find(indices==id))
                indices=[indices id]; 
            end;
        end;
    else % processing for single tag (no range)
        id=tag2id(tag,tags,codelist);
        if ~isempty(id) && ~sum(find(indices==id)), 
            indices=[indices id]; 
        end;
    end;
end;


