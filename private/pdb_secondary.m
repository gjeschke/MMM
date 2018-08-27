function message=pdb_secondary(snum,cnum)
% function message=pdb_secondary(sid,cid)
%
% Creates secondary structure definitions in an MMM model
% from the information that was present in the original PDB file,
% any user-defined definitions are overwritten!
% the transmembrane flag (TM) for helices and strands is set to 0
%
% snum      structure index
% cnum      chain index
%
% message   error message structure with fields
%           .error
%           .text
%
% G. Jeschke, 2009

global model

message.error=0;
message.text='';

success=0;

loops=0;
helices=0;
strands=0;

model.structures{snum}(cnum).loop_defs={};
model.structures{snum}(cnum).helix_defs={};
model.structures{snum}(cnum).sheet_defs={};

sequence=model.structures{snum}(cnum).sequence;
resnum=0;
% it is assumed
% that PDB helix definitions are usually nonsense
    indices=[snum,cnum,1];
    n=length(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info);
    CA=zeros(n,3);
    O=zeros(n,3);
    poi=0;
    cpoi=0;
    if ~isempty(model.structures{indices(1)}(indices(2)).residues{indices(3)}.info),
        starter=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(1).number;
        sec=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(1).secondary;
        for k=1:n,
            info=model.structures{indices(1)}(indices(2)).residues{indices(3)}.info(k);
            oldnum=resnum;
            resnum=info.number;
            if resnum>length(sequence),
                resnum=oldnum;
                break
            end;
            CA_id=tag2id('CA',info.atom_tags);
            O_id=tag2id('O',info.atom_tags);
            secend=1;
            if ~isempty(CA_id) && ~isempty(O_id),
                poi=poi+1;
                secend=0;
            end;
            if info.secondary~=sec && poi~=0,
                cpoi=cpoi+1;
                switch sec
                    case 0
                        loops=loops+1;
                        model.structures{snum}(cnum).loop_defs{loops}.name=sprintf('%i',loops);
                        model.structures{snum}(cnum).loop_defs{loops}.range=[starter,resnum-1];
                    case 1
                        helices=helices+1;
                        model.structures{snum}(cnum).helix_defs{helices}.name=sprintf('%i',helices);
                        model.structures{snum}(cnum).helix_defs{helices}.range=[starter,resnum-1];
                        model.structures{snum}(cnum).helix_defs{helices}.TM=0;
                    case 2
                        strands=strands+1;
                        model.structures{snum}(cnum).sheet_defs{strands}.name=sprintf('%i',strands);
                        model.structures{snum}(cnum).sheet_defs{strands}.range=[starter,resnum-1];
                        model.structures{snum}(cnum).sheet_defs{strands}.TM=0;
                end;
                starter=resnum;
                if secend,
                    poi=0;
                else
                    poi=1;
                end;
                sec=info.secondary;
                % disp(sprintf('Residue %i, secondary ID: %i',k,sec));
            end;
        end;
        if poi>0,
            switch sec
                case 0
                    loops=loops+1;
                    model.structures{snum}(cnum).loop_defs{loops}.name=sprintf('%i',loops);
                    model.structures{snum}(cnum).loop_defs{loops}.range=[starter,resnum];
                case 1
                    helices=helices+1;
                    model.structures{snum}(cnum).helix_defs{helices}.name=sprintf('%i',helices);
                    model.structures{snum}(cnum).helix_defs{helices}.range=[starter,resnum];
                    model.structures{snum}(cnum).helix_defs{helices}.TM=0;
                case 2
                    strands=strands+1;
                    model.structures{snum}(cnum).sheet_defs{strands}.name=sprintf('%i',strands);
                    model.structures{snum}(cnum).sheet_defs{strands}.range=[starter,resnum];
                    model.structures{snum}(cnum).sheet_defs{strands}.TM=0;
            end;
        end;
        model.structures{snum}(cnum).loop_defs=model.structures{snum}(cnum).loop_defs(1:loops);
        model.structures{snum}(cnum).helix_defs=model.structures{snum}(cnum).helix_defs(1:helices);
        model.structures{snum}(cnum).sheet_defs=model.structures{snum}(cnum).sheet_defs(1:strands);
    end;
% end;

if ~success,
    message.error=1;
    message.text='No secondary structure defined.';
end;
