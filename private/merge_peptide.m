function indices = merge_peptide(chain_model_1,res_range_1,chain_model_2,res_range_2,tag,indices)
% indices = merge_peptide(chain_model_1,res_range_1,chain_model_2,res_range_2,indices)
%
% Merges two peptide sections into a single chain
%
% chain_model_1     indices [1,3] of chain model of first section
% res_range_1       residue numbers [1,2] of the first section
% chain_model_2     indices [1,3] of chain model of first section
% res_range_2       residue numbers [1,2] of the first section
% tag               chain tag, defaults to A
% indices           chain model indices for the merged model (optional,
%                   otherwise a new structure is generated)
%
% G. Jeschke, 5.2.2019

global model

if ~exist('tag','var')
    tag = 'A';
end

if ~exist('indices','var')
    indices = ones(1,3);
    indices(1) = length(model.structures) + 1;
    model.structure(indices(1)).
end

model.structure(chain_id).name=tag;
chain_tags=[chain_tags tag ':'];
terminated_chains=[terminated_chains tag ':'];
chain_ids=[chain_ids chain_id];
structure(chain_id).seqtype=0;
structure(chain_id).modified=0;
structure(chain_id).nonstandard=0;
structure(chain_id).helices=0;
structure(chain_id).strands=0;
structure(chain_id).xyz{1}=zeros(max_atoms,3);
structure(chain_id).isotopes=zeros(max_atoms,2,'single');
structure(chain_id).Bfactor{1}=zeros(1,max_atoms);
structure(chain_id).Btensor{1}=zeros(max_atoms,6,'int32');
structure(chain_id).atoms{1}=0;
structure(chain_id).conn=zeros(max_atoms,30,'int32');
structure(chain_id).maxconn=1;
structure(chain_id).residues{1}.residue_tags=':';
structure(chain_id).residues{1}.info=[];

