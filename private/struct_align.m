function struct_align(fname,ind1,ind2,r_match)
% function struct_align(snum1,snum2)
%
% attempts to create and write out a sequence alignment based on a
% structural alignment of the two chains with indices ind1 and ind2
% the alignemnt is written to files <fname>.ali (PIR format) and
% <fname>.clu (CLUSTAL-like alignment format)
%
% fname     basis file name (no extension) for output files
% ind1      index vector of length 2 (optionally 3) for structure 1
% ind2      index vector of length 2 (optionally 3) for structure 2
% r_match   C_alpha-C_alpha distance in Å that constitutes a structural match
%
% G. Jeschke, 2011

global residue_defs
global model

if nargin<4,
    r_match=10; % C_alpha-C_alpha distance in Å that constitutes a structural match
end;


if length(ind1)<3, 
    ind1=[ind1 1]; 
end;

if length(ind2)<3, 
   ind2=[ind2 1]; 
end;

adr1=[mk_address(ind1) '.CA'];
indices1=resolve_address(adr1);

adr2=[mk_address(ind2) '.CA'];
indices2=resolve_address(adr2);

if isempty(indices1) || isempty(indices2), 
    return; 
end; % may happen for DNA, RNA, generally non-proteins
[m1,n]=size(indices1); % m1 is the number of residues
[m2,n]=size(indices2);

% extract coordinates and single-letter codes
Ca_coor1=zeros(m1,3);
Ca_coor2=zeros(m2,3);
seq1=zeros(1,m1);
seq2=zeros(1,m2);
num1=zeros(1,m1);
num2=zeros(1,m2);

poi=0;
for k=1:m1,
    info=model.structures{indices1(k,1)}(indices1(k,2)).residues{indices1(k,3)}.info(indices1(k,4));
    if info.type==1, % only amino acids are considered
        poi=poi+1;
        [msg,coor]=get_atom(indices1(k,:),'coor');
        Ca_coor1(poi,:)=coor;
        tlc=model.structures{indices1(k,1)}(indices1(k,2)).residues{indices1(k,3)}.info(indices1(k,4)).name;
        seq1(poi)=tag2id(upper(tlc),upper(residue_defs.restags),residue_defs.single_letter_code);
        num1(poi)=model.structures{indices1(k,1)}(indices1(k,2)).residues{indices1(k,3)}.info(indices1(k,4)).number;
    end;
end;
seq1=char(seq1(1:poi));
Ca_coor1=Ca_coor1(1:poi,:);
num1=num1(1:poi);
m1=poi;
[snum1,rind1]=sort(num1);
seq1=seq1(rind1);
Ca_coor1=Ca_coor1(rind1,:);

poi=0;
for k=1:m2,
    info=model.structures{indices2(k,1)}(indices2(k,2)).residues{indices2(k,3)}.info(indices2(k,4));
    if info.type==1, % only amino acids are considered
        poi=poi+1;
        [msg,coor]=get_atom(indices2(k,:),'coor');
        Ca_coor2(poi,:)=coor;
        tlc=model.structures{indices2(k,1)}(indices2(k,2)).residues{indices2(k,3)}.info(indices2(k,4)).name;
        seq2(poi)=tag2id(upper(tlc),upper(residue_defs.restags),residue_defs.single_letter_code);
        num2(poi)=model.structures{indices2(k,1)}(indices2(k,2)).residues{indices2(k,3)}.info(indices2(k,4)).number;
    end;
end;
seq2=char(seq2(1:poi));
Ca_coor2=Ca_coor2(1:poi,:);
num2=num2(1:poi);
m2=poi;
[snum2,rind2]=sort(num2);
seq2=seq2(rind2);
Ca_coor2=Ca_coor2(rind2,:);

eseq1='';
eseq2='';
poi1=0;
poi2=0;
cpoi1=0;
cpoi2=0;
epoi=0;
shift=0;
while poi1<m1,
    poi1=poi1+1;
    coor1=Ca_coor1(poi1,:);
    cdiff2=Ca_coor2-repmat(coor1,m2,1);
    dist=sqrt(sum(cdiff2.^2,2));
    [match,poi2]=min(dist);
    if match<r_match,
        coor2=Ca_coor2(poi2,:);
        cdiff1=Ca_coor1-repmat(coor2,m1,1);
        dist=sqrt(sum(cdiff1.^2,2));
        [match2,mpoi2]=min(dist);
        if mpoi2==poi1, % the two residues are within structural matching range
                        % and no other residue in the other structure is
                        % closer than the paired residue          
            epoi=epoi+1;
            cshift=poi1-poi2;
            while cshift<shift,
                eseq1=[eseq1 '-'];
                cpoi2=cpoi2+1;
                eseq2=[eseq2 seq2(cpoi2)];
                shift=shift-1;
            end;
            while cshift>shift,
                eseq2=[eseq2 '-'];
                cpoi1=cpoi1+1;
                eseq1=[eseq1 seq1(cpoi1)];
                shift=shift+1;
            end;
            eseq1=[eseq1 seq1(poi1)];
            cpoi1=cpoi1+1;
            eseq2=[eseq2 seq2(poi2)];
            cpoi2=cpoi2+1;           
        else
            eseq1=[eseq1 seq1(poi1)];
            eseq2=[eseq2 '-'];
            shift=shift+1;
        end;
    else
        eseq1=[eseq1 seq1(poi1)];
        eseq2=[eseq2 '-'];
        shift=shift+1;
    end;
end;

