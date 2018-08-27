function [transmat,target,moving,alpha]=opt_superposition_ensemble(move_adr,targ_adr,options)
% Optimum coordinate superposition of structure with address move_adr onto
% structure with address targ_adr
%
% move_adr  address of the moving structure
% targ_adr  address of the target structure, if the target structure has
%           several models, superposition is onto the mean coordinates
% options   structure with fields
%           .mode           atom types per residue to be superimposed
%                           'all', 'backbone' or 'Calpha'
%           .correspondence 1 superimpose only residues with addresses that
%                           differ only in the structure tag
%                           0 superimpose all residues
%           .selected       1 superimpose only selected residues, this
%                           requires that the same number of residues is
%                           selected in both structures
%                           0 superimpose all residues
%           .aligned        1 superimpose residues based on sequence
%                           alignment, only selected chains
%                           0 no sequence alignment attempt is made
%                           if chains are selected, only selected chains
%                           are superimposed, in that case, the same number
%                           of chains must be selected in both structures
%           .whole          1 superimpose whole structure by alignment, 
%                           requires same number of chains in both
%                           structures
% Warnings: 
% - out of the flags selected, aligned, and whole, only one should be set
% not all combinations are tested with ensemble mode, output vector alpha
% is defined only for suggested usage below
% - ensemble size must be the same for all chains in a multi-chain
% structure
% suggested usage:  .mode='Calpha'; correspondence=0; .selected=0;
%                   .aligned=1; .whole=1;
% options defaults to suggested usage
% 
% transmat  4x4 affine transformation matrix that transforms the moving
%           structure onto the (mean) target structure
% target    target structure mean Calpha coordinates, if
%           options.mode='Calpha', empty array otherwise
% moving    moving structure Calpha coordinates, if
%           options.mode='Calpha', empty array otherwise
% alpha     array of certainties of the motion vectors (target-moving) of
%           all Calpha atoms, if options.mode='Calpha', empty array otherwise
%           if the moving structure is not an ensemble, all entries are 1
%           for an ensemble, sum of the normalized scalar products between
%           the mean displacement vector and the displacement vectors of
%           the individual structures in the ensemble. The sum is
%           normalized to the number of structures in the ensemble.
%           ranges between 0 (completely uncertain direction) and 1 (fully
%           certain direction), alpha is empty array unless options
%           conforms to suggested usage, while transmat, target, and moving
%           are still defined
%           
% G. Jeschke, 2011

global model

if nargin<3,
    options.mode='Calpha';
    options.correspondence=0;
    options.selected=0;
    options.aligned=0;
    options.whole=1;
end;

alpha_defined=true;
if ~strcmpi(options.mode,'Calpha'),
    alpha_defined=false;
end;
if options.correspondence~=0,
    alpha_defined=false;
end;
if options.selected~=0,
    alpha_defined=false;
end;
if options.whole~=1,
    alpha_defined=false;
end;

transmat=eye(4,4);
target=[];
moving=[];
alpha=[];

% id_template=tag2id(targ_adr,model.structure_tags,model.structure_ids);
% id_trafo=tag2id(move_adr,model.structure_tags,model.structure_ids);

id_template=resolve_address(targ_adr);
id_trafo=resolve_address(move_adr);

if isempty(id_template),
    add_msg_board('ERROR: Starting structure not found (address wrong).');
    return
else 
    id_template=id_template(1);
end;
            
if isempty(id_trafo),
    add_msg_board('ERROR: End structure not found (address wrong).');
    return
else
    id_trafo=id_trafo(1);
end;

esize=1;
if options.selected || options.aligned || options.whole, % in all these cases arrays sel1 and sel2 of selected residues are generated
    if ~options.aligned && ~options.whole, % user selection
        [selections,msg]=resolve_address('*');
        [msel,nsel]=size(selections);
        poi1=0;
        poi2=0;
        sel1=zeros(msel,4);
        sel2=zeros(msel,4);
        if msel>0,
            for k=1:msel,
                cindices=selections(k,:);
                cindices=cindices(cindices>0);
                if length(cindices)==4,
                    if selections(k,1)==id_template,
                            poi1=poi1+1;
                            sel1(poi1,:)=cindices;
                    end;
                    if selections(k,1)==id_trafo,
                        poi2=poi2+1;
                        sel2(poi2,:)=cindices;
                    end;
                end;
            end;
        end;
    end;
    cid1=model.chain_ids(id_template);
    cid1=cid1{1};
    cid2=model.chain_ids(id_trafo);
    cid2=cid2{1};
    if options.whole, % try to align whole structures
        if length(cid1)~=length(cid2),
            add_msg_board('ERROR: Different number of chains in target and moving structure.');
            add_msg_board('Deactivate "whole structure".');
            return
        end;
        sel1=zeros(2000,4);
        sel2=sel1;
        psel=0;
        for k=1:length(cid1),
            seqs{1}=model.structures{id_template}(cid1(k)).sequence;
            seqs{2}=model.structures{id_trafo}(cid2(k)).sequence;
            sindices=[id_template,cid1(k);id_trafo,cid2(k)];
            [message,inname]=align_sequences(seqs,sindices,true);
            if message.error,
                add_msg_board('ERROR: MUSCLE sequence alignment failed.');
                add_msg_board(message.text);
                add_msg_board('Deactivate "whole structure".');
                return
            end;
            alignment=get_multiple_clustal(inname);
            [asel1,asel2,esize]=select_aligned(alignment,sindices);
            [msel,nsel]=size(asel1);
            sel1(psel+1:psel+msel,:)=asel1;
            sel2(psel+1:psel+msel,:)=asel2;
            psel=psel+msel;
        end;
        sel1=sel1(1:psel,:);
        sel2=sel2(1:psel,:);
        poi1=psel;
        poi2=psel;
    elseif options.aligned, % we have to establish if there is only one chain per structure
                            % or if the same number of chains per structure is selected
        if length(cid1)==1 && length(cid2)==1,
            seqs{1}=model.structures{id_template}(cid1).sequence;
            seqs{2}=model.structures{id_trafo}(cid2).sequence;
            sindices=[id_template,cid1;id_trafo,cid2];
            [message,inname]=align_sequences(seqs,sindices,true);
            if message.error,
                add_msg_board('ERROR: MUSCLE sequence alignment failed.');
                add_msg_board(message.text);
                add_msg_board('Deactivate "whole structure".');
                return
            end;
            alignment=get_multiple_clustal(inname);
            [sel1,sel2]=select_aligned(alignment,sindices);
            [poi1,nsel]=size(sel1);
            poi2=poi1;
        else
            [selections,msg]=resolve_address('*');
            [msel,nsel]=size(selections);
            csel1=zeros(msel,2);
            csel2=zeros(msel,2);
            poi1=0;
            poi2=0;
            for k=1:msel,
                cindices=selections(k,:);
                cindices=cindices(cindices>0);
                if length(cindices)==2,
                    if selections(k,1)==id_template,
                        poi1=poi1+1;
                        csel1(poi1,:)=cindices;
                    end;
                    if selections(k,1)==id_trafo,
                        poi2=poi2+1;
                        csel2(poi2,:)=cindices;
                    end;
                end;
            end;
            csel1=csel1(1:poi1,:);
            csel2=csel2(1:poi2,:);
            if poi1~=poi2,
                add_msg_board('ERROR: Different number of chains selected in target and moving structure.');
                add_msg_board('Deactivate "aligned" to fit whole structure or change selection.');
                return
            end;
            sel1=zeros(2000,4);
            sel2=sel1;
            psel=0;
            for k=1:poi1,
                seqs{1}=model.structures{csel1(k,1)}(csel1(k,2)).sequence;
                seqs{2}=model.structures{csel2(k,1)}(csel2(k,2)).sequence;
                sindices=[csel1(k,:);csel2(k,:)];
                [message,inname]=align_sequences(seqs,sindices,true);
                if message.error,
                    add_msg_board('ERROR: MUSCLE sequence alignment failed.');
                    add_msg_board(message.text);
                    add_msg_board('Deactivate "whole structure".');
                    return
                end;
                alignment=get_multiple_clustal(inname);
                [asel1,asel2]=select_aligned(alignment,sindices);
                [msel,nsel]=size(asel1);
                sel1(psel+1:psel+msel,:)=asel1;
                sel2(psel+1:psel+msel,:)=asel2;
                psel=psel+msel;
            end;
            sel1=sel1(1:psel,:);
            sel2=sel2(1:psel,:);
            poi1=psel;
            poi2=psel;
        end;
    end;
    if poi1==0,
        add_msg_board('ERROR: No residues selected or aligned in target structure.');
        add_msg_board('Deactivate "only selected" or "aligned" to fit whole structure.');
        return
    end;
    if poi2==0,
        add_msg_board('ERROR: No residues selected or aligned in structure to be fitted.');
        add_msg_board('Deactivate "only selected" or "aligned" to fit whole structure.');
        return
    end;
    sel1=sel1(1:poi1,:);
    sel2=sel2(1:poi2,:);
    coor0=zeros(20000,3);
    coor1=coor0;
    if esize>1,
        coore=zeros(esize,20000,3);
    end;
    poi=0;
    test_corr=options.correspondence;
    if ~test_corr && poi1~=poi2,
        add_msg_board('ERROR: Different number of selected residues and no correspondence check.');
        return
    end;
    for k1=1:poi1,
        cindices1=sel1(k1,:);
        [stag1,ctag1,modelnum1,resnum1,icode1]=mk_address_parts(cindices1);
        cindices2=[];
        if test_corr,
            for k2=1:poi2,
                tindices=sel2(k2,:);
                [stag2,ctag2,modelnum2,resnum2,icode2]=mk_address_parts(tindices);
                if strcmpi(ctag1,ctag2) && resnum1==resnum2 && strcmpi(icode1,icode2),
                    cindices2=sel2(k2,:);
                end;
            end;
        else
            cindices2=sel2(k1,:);
        end;
        if isempty(cindices2),
        else
            adr1=mk_address(cindices1);
            switch options.mode
                case 'all'
                    adr1=sprintf('%s.:',adr1);
                case 'backbone'
                    adr1=sprintf('%s.N,CA,C,O',adr1);
                case 'Calpha'
                    adr1=sprintf('%s.CA',adr1);
            end;
            [indices,message]=resolve_address(adr1);

            [m,n]=size(indices);
            for k=1:m,
                cindices0=indices(k,:);
                cindices0=cindices0(1:5);
                address=mk_address(cindices0);
                postpoi=findstr(address,'.');
                post=address(postpoi:end);
                [stag2,ctag2,modelnum2,resnum2,icode2]=mk_address_parts(cindices2);
                adr2=sprintf('[%s](%s){%i}%i%s%s',stag2,ctag2,modelnum2,resnum2,icode2,post);
                [cindices,msg]=resolve_address(adr2);
                if msg.error==0 && ~isempty(cindices),
                    poi=poi+1;
                    [m2,n2]=size(cindices);
                    [msg,cc0]=get_atom(cindices0,'coor');
                    [msg,cc1]=get_atom(cindices,'coor');
                    if ~isempty(cc0) && ~isempty(cc1),
                        coor0(poi,:)=cc0;
                        coor1(poi,:)=cc1;
                    else
                        poi=poi-1;
                    end;
                    if esize>1,
                        for em=1:esize,
                            cindicese=cindices0;
                            cindicese(3)=em;
                            [msg,cce]=get_atom(cindicese,'coor');
                            if ~isempty(cce)
                                coore(em,poi,:)=cce;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    coor0=coor0(1:poi,:);
    coor1=coor1(1:poi,:);
    if esize>1,
        coore=coore(:,1:poi,:);
    end;
    add_msg_board(sprintf('%i selected atoms are fitted',poi));
else
    adr1=sprintf('[%s](:).',targ_adr);
    switch options.mode
        case 'all'
            adr1=sprintf('%s:',adr1);
        case 'backbone'
            adr1=sprintf('%sN,CA,C,O',adr1);
        case 'Calpha'
            adr1=sprintf('%sCA',adr1);
    end;

    [indices,message]=resolve_address(adr1);

    [m,n]=size(indices);
    
    coor0=zeros(m,3);
    coor1=coor0;

    poi=0;
    poi2=0;
    for k=1:m,
        cindices0=indices(k,:);
        cindices0=cindices0(1:5);
        address=mk_address(cindices0);
        prepoi=findstr(address,'[');
        pre=address(1:prepoi);
        postpoi=findstr(address,']');
        post=address(postpoi:end);
        adr2=[pre move_adr post];
        [cindices,msg]=resolve_address(adr2);
        if msg.error==0 && ~isempty(cindices),
            poi=poi+1;
            [m2,n2]=size(cindices);
            if m2>1,
                poi2=poi2+1;
            end;
            [msg,cc0]=get_atom(cindices0,'coor');
            [msg,cc1]=get_atom(cindices,'coor');
            if ~isempty(cc0) && ~isempty(cc1),
                coor0(poi,:)=cc0;
                coor1(poi,:)=cc1;
            else
                poi=poi-1;
            end;
        end;
    end;
    coor0=coor0(1:poi,:);
    coor1=coor1(1:poi,:);
    add_msg_board(sprintf('%i atoms were selected',m));
    add_msg_board(sprintf('of which %i atoms were found in 2nd structure',poi));
    if poi2>0,
        add_msg_board(sprintf('with %i of them being ambiguous',poi2));
    end;
end;

if esize>1,
    % make average ensemble coordinates
    coor0=zeros(size(coor0));
    [m,n]=size(coor0);
    for em=1:esize,       
        coor=reshape(coore(em,:,:),m,n);
        coor0=coor0+coor;
    end;
    coor0=coor0/esize;
end;

[m,n]=size(coor0);
alpha=ones(m,1);
coor2b=[];
if m>0 && n>0,
    [rmsd,coor0b,transmat]=rmsd_superimpose(coor1,coor0);
    add_msg_board(sprintf('RMSD of both structures is: %4.2f Å',rmsd));
    if esize>1,
        alpha=zeros(size(alpha));
        motion=coor0b-coor1;
        for em=1:esize,
            coor0c=coore(em,:,:);
            coor0c=reshape(coor0c,m,n);
            [rmsd,coor0d]=rmsd_superimpose(coor1,coor0c);
            motionc=coor0d-coor1;
            for k=1:m,
                move=motion(k,:)/norm(motion(k,:));
                movec=motionc(k,:)/norm(motionc(k,:));
                alpha(k)=alpha(k)+dot(move,movec);
            end;
        end;
        alpha=alpha/esize;
        alpha=(pi/2-acos(alpha))/(pi/2);
        for k=1:length(alpha),
            if alpha(k)<0, alpha(k)=0; end;
        end;
        alpha=alpha.^2;
    end;
else
    add_msg_board('Warning! No transformation performed as selections do not match.');
    transmat=[];
end;

% store Calpha coordinates for move display
if strcmpi(options.mode,'Calpha'),
    moving=coor1;
    target=coor0b;
    if ~alpha_defined,
        alpha=[];
    end;
else
    moving=[];
    target=[];
    alpha=[];
end;


function [sel1,sel2,esize]=select_aligned(alignment,sindices)
% select only residues that are aligned

global model

% determine ensemble size for moving structure
esize=length(model.structures{sindices(1,1)}(sindices(1,2)).residues);

seq1=alignment(1).sequence;
seq2=alignment(2).sequence;
sel1=zeros(length(seq1),4);
sel2=sel1;
poi1=0;
poi2=0;
pois=0;
for k=1:length(seq1),
    match1=false;
    match2=false;
    if char(seq1(k))~='-',
        poi1=poi1+1;
        match1=true;
    end;
    if char(seq2(k))~='-',
        poi2=poi2+1;
        match2=true;
    end;
    if match1 && match2
        tag1=sprintf('%i',poi1);
        id1=tag2id(tag1,model.structures{sindices(1,1)}(sindices(1,2)).residues{1}.residue_tags);
        tag2=sprintf('%i',poi2);
        id2=tag2id(tag2,model.structures{sindices(2,1)}(sindices(2,2)).residues{1}.residue_tags);
        if ~isempty(id1) && ~isempty(id2),
            pois=pois+1;
            sel1(pois,:)=[sindices(1,:) 1 id1];
            sel2(pois,:)=[sindices(2,:) 1 id2];
        end;
    end;
end;
sel1=sel1(1:pois,:);
sel2=sel2(1:pois,:);

