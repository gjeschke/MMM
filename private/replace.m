function [repnum,msg]=replace(snum,replacements,selected,showit)
%function [repnum,msg]=replace(snum,replacements)
%   Replaces non-native amino acids by their native equivalents
%
% snum          number of the structure in which sidechains are to be replaced
%               if optional argument selected is true, snum is also
%               interpreted as an index array: residues matching an index
%               vector in that array are considered to be selected
% replacements  replacement list in lag list format, e.g. ':CSE:MSE:' to
%               replace selenocysteine by cysteine and selenomethionine by
%               methionine, ':' selects all non-native amino acids
% selected      optional flag, indicates that only selected residues are to
%               be replaced, 0 replace all matching residues, 1 replace
%               only selected residues, defaults to 0 (false)
%               both MMM's selection table and input indices snmum are
%               considered as selected objects
% showit        (optional) flag suppressing 3D model update, 0 no graphics
%               output, 1 graphics output, defaults to 1
%
% repnum        number of replaced amino acids
% msg           error message
%

global model
global residue_defs
global hMain


if nargin<3,
    selected=false;
end;

if nargin<4,
    showit=true;
end;

msg.error=0;
msg.text='No error.';

if strcmp(replacements,':::'),
    replacements=residue_defs.hettags;
end;

repnum=0;

if isempty(snum),
    msg.error=1;
    msg.text='No structure selected.';
else
    iindices=snum;
    snum=snum(1,1);
    [mi,ni]=size(iindices);
end;

stag=id2tag(snum,model.structure_tags);
chains=length(model.structures{snum});
add_msg_board(sprintf('Replacing non-native amino acids in %i chains of structure %s',chains,stag));

for cnum=1:chains,
    models=length(model.structures{snum}(cnum).residues);
    for mnum=1:models,
        info=model.structures{snum}(cnum).residues{mnum}.info;
        residues=length(info);
        for rnum=1:residues,
            rtag=info(rnum).name;
            id=tag2id(rtag,replacements);
            if ~isempty(id),
                indices=[snum cnum mnum rnum];
                if ~selected || indexmatch(indices,iindices) || isselected(indices),
                    adr=mk_address(indices,1);
                    add_msg_board(sprintf('Replacing sidechain of residue %s',adr));
                    if showit,
                        cmd(hMain,sprintf('hide %s',adr));
                    end;
                    info=model.structures{snum}(cnum).residues{mnum}.info;
                    repnum=repnum+1;
                    info=replace_residue(info,indices,upper(rtag));
                    model.structures{snum}(cnum).residues{mnum}.info=info;
                end;
            end;
        end;
    end;
    adr=mk_address([snum,cnum]);
    if showit,
        cmd(hMain,sprintf('hide %s',adr));
    end;
    message=consolidate_chain([snum,cnum]);
    if showit,
        cmd(hMain,sprintf('show %s ribbon',adr));
    end;
end;




    function info=replace_residue(info,indices,rtag)
        % replaces a single residue by the corresponding native residue
        rnum1=indices(4);
        atom_tags0=info(rnum1).atom_tags;
        elements0=info(rnum1).elements;
        atom_numbers0=info(rnum1).atom_numbers;
        switch rtag
            case {'MTN','R1A','IA1','R1B','R1F','R7A','V1A'}
                [atom_tags,elements,atom_numbers]=condense(atom_tags0,elements0,atom_numbers0,':C:N:O:CA:CB:SG:',':C:N:O:CA:CB:SG:');
                info(rnum).name='CYS';
                set_sequence(indices,'C','CYS');
                unset_mutation(indices,'CYS','cysteine');
            case 'CSE'
                [atom_tags,elements,atom_numbers]=condense(atom_tags0,elements0,atom_numbers0,':C:N:O:CA:CB:SE:',':C:N:O:CA:CB:SG:');
                info(rnum).name='CYS';
                set_sequence(indices,'C','CYS');
                unset_mutation(indices,'CYS','cysteine');
            case 'MSE'
                [atom_tags,elements,atom_numbers]=condense(atom_tags0,elements0,atom_numbers0,':C:N:O:CA:CB:CG:SE:CE:',':C:N:O:CA:CB:CG:SD:CE:');
                info(rnum).name='MET';
                set_sequence(indices,'M','MET');
                unset_mutation(indices,'MET','methionine');
        end;
        info(rnum).atom_tags=atom_tags;
        info(rnum).elements=elements;
        info(rnum).atom_numbers=atom_numbers;
        info(rnum).hetflag=false;
        adr0=mk_address(indices);
        if isfield(model,'labels'),
            labels0=model.labels;
            poi=0;
            for k=1:length(model.labels),
                if ~strcmp(model.labels(k).adr,adr0),
                    poi=poi+1;
                    labels0(poi)=labels0(k);
                end;
            end;
            labels0=labels0(1:poi);
            model.labels=labels0;
        end;
    end

    function [atom_tags,elements,atom_numbers]=condense(atom_tags0,elements0,atom_numbers0,maintain,newtags)
        % condenses the record of a residue to a abbreviated list of maintained atoms with new tags 
        
        atom_tags=newtags;
        nonsense=textscan(maintain,'%s','Delimiter',':');
        atag_list0=nonsense{1};
        atag_list=cell(1,length(atag_list0)-1);
        for k=2:length(atag_list0),
            atag_list{k-1}=atag_list0{k};
        end;
        elements=zeros(1,length(atag_list));
        atom_numbers=cell(1,length(atag_list));
        atlist=zeros(1,length(atag_list));
        for anum=1:length(atag_list), % condense lists
            atag=atag_list{anum};
            aind=tag2id(atag,atom_tags0);
            elements(anum)=elements0(aind);
            atnum=atom_numbers0{aind};
            atom_numbers{anum}=atnum(1,1);
            atlist(anum)=atnum(1,1);
        end;
        conn=model.structures{indices(1)}(indices(2)).conn;
        for anum=1:length(atlist), % correct connection (bond) list
            oconn=conn(atlist(anum),:); % old connections
            poi=0;
            nconn=zeros(size(oconn));
            for k=1:length(oconn),
                if sum(oconn(k)==atlist), % bonds only to atoms that still exist
                    poi=poi+1;
                    nconn(poi)=oconn(k);
                end;
            end;
            conn(atlist(anum),:)=nconn;
        end;
        model.structures{indices(1)}(indices(2)).conn=conn;
        coors=model.structures{indices(1)}(indices(2)).xyz{indices(3)};
        for k=1:length(elements),
            if elements(k)==34,
                elements(k)=16;
                atnum_Se=atom_numbers{k};
                xyz_Se=coors(atnum_Se(1,1),:);
                atnum_prev=atom_numbers{k-1};
                xyz_prev=coors(atnum_prev(1,1),:);
                vec=-0.1*(xyz_Se-xyz_prev); % S-C bond is 10% shorter than Se-C bond
                atnum_next=atom_numbers{k+1};
                xyz_next=coors(atnum_next(1,1),:);
                vec2=-0.1*(xyz_next-xyz_Se); % S-C bond is 10% shorter than Se-C bond
                for kk=k:length(elements),
                    atnum=atom_numbers{kk};
                    xyz=coors(atnum(1,1),:);
                    coors(atnum(1,1),:)=xyz+vec;
                end;
                for kk=k+1:length(elements),
                    atnum=atom_numbers{kk};
                    xyz=coors(atnum(1,1),:);
                    coors(atnum(1,1),:)=xyz+vec2;
                end;
            end;
        end;
        model.structures{indices(1)}(indices(2)).xyz{indices(3)}=coors;
    end

    function set_sequence(indices,slc,tlc)
        % updates single-letter code in the sequence record
    
        rnum2=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number;
        seq=model.structures{indices(1)}(indices(2)).sequence;
        if indices(4)>=1 && indices(4)<=length(seq),
            seq(rnum2)=slc;
            poi=2+4*(rnum2-1);
            model.structures{indices(1)}(indices(2)).restags(poi:poi+2)=tlc;
        end;
        model.structures{indices(1)}(indices(2)).sequence=seq;
    end

    function unset_mutation(indices,newtag,comment)
        model.structures{indices(1)}(indices(2)).modified=model.structures{indices(1)}(indices(2)).modified-1;
        if ~isfield(model.structures{indices(1)}(indices(2)),'mutations'),
            return;
        end;
        oldmutations=model.structures{indices(1)}(indices(2)).mutations;
        mutpoi=0;
        for k=1:length(oldmutations),
            if oldmutations(k).number~=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number,
                mutpoi=mutpoi+1;
                mutations(mutpoi)=oldmutations(k);
            elseif ~strcmpi(oldmutations(k).original,newtag)
                mutations(mutpoi)=oldmutations(k);
                mutations(mutpoi).modified=newtag;
                mutations(mutpoi).comment=comment;
            end;
        end;
        if exist('mutations','var'),
            model.structures{indices(1)}(indices(2)).mutations=mutations;
        else
            model.structures{indices(1)}(indices(2)).mutations=[];
        end;
        oldhet=model.structures{indices(1)}(indices(2)).het;
        mutpoi=0;
        for k=1:length(oldhet),
            if oldhet(k).number~=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(indices(4)).number,
                mutpoi=mutpoi+1;
                het(mutpoi)=oldhet(k);
            end;
        end;
        if exist('het','var'),
            model.structures{indices(1)}(indices(2)).het=het;
        else
            model.structures{indices(1)}(indices(2)).het=[];
        end;
    end


end

