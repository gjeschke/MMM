function [resaxis,s,npN,npC,meanct,uncert,msg] = analyze_local_order(indices,options)
% [resaxis,s,npN,npC] = analyze_local_order(indices)
%
% Compute local order parameters for a peptide chain
%
% all residues that have an atom of type CA are considered amino acids, all
% amino acid residues must also have atoms of type N and C
%
% indices   array (1,2) of structure and chain index
% options   further processing options, defualts to no further processing,
%           structure with fields
%           .superimpose    if true, conformers are superimposed at the
%                           residue with maximal local order parameter
% resaxis   residue number axis
% s         vector of local order parameters with range [0,1] for all 
%           residues, NaN, if undefined for this residue
% npN       vector of persistence numbers towards N terminus for all 
%           residues, NaN, if undefined for this residue
% npC       vector of persistence numbers towards C terminus for all 
%           residues, NaN, if undefined for this residue
% meanct    ensemble mean of cos(theta_n) for all residues
% uncert    uncertainty estimate (error bar) for s
% msg       error message with fields
%           msg.error   error code, 0 for no error
%           msg.text    error message text
%
% G. Jeschke, 6.1.2019

global model

msg.error = 0;
msg.text = 'OK';

if ~exist('options','var') || ~isfield(options,'superimpose')
    options.superimpose = true;
end

snum = indices(1);
cnum = indices(2);

models = length(model.structures{snum}(cnum).residues);
residues = length(model.structures{snum}(cnum).residues{1}.info);
    
s = nan(1,residues);
npN = nan(1,residues);
npC = nan(1,residues);
uncert = nan(2,residues);

if models < 2 % if there is only one conformation, local order parameters are undefined
    msg.error = 1;
    msg.text = 'No ensemble';
    return
end

aa = zeros(1,residues);
resaxis = zeros(1,residues);
CA_traces = cell(1,models);
CA_traces_LF = cell(1,models);

rng(13);
selection = rand(1,models);
% extract Calpha traces
for mnum=1:models
    aapoi = 0;
    curr_CA_trace = zeros(residues,3);
    for rnum=1:residues
        adr = mk_address([snum,cnum,mnum,rnum]);
        indCA = resolve_address([adr '.CA']);
        [~,coorCA] = get_object(indCA,'coor');
        if isempty(coorCA)
            continue;
        end
        aapoi = aapoi + 1; % if the Calpha atom is present, the residue is considered an amino acid
        if mnum == 1
            aa(aapoi) = rnum;
            resaxis(aapoi) = model.structures{snum}(cnum).residues{mnum}.info(rnum).number;
        else
            if aa(aapoi) ~= rnum
                msg.error = 2;
                msg.text = 'Residue mismatch between conformations';
            end
        end
        curr_CA_trace(aapoi,:) = coorCA;
    end
    curr_CA_trace = curr_CA_trace(1:aapoi,:);
    CA_traces{mnum} = curr_CA_trace;
end
aa = aa(1:aapoi);
resaxis = resaxis(1:aapoi);
s = s(1:aapoi);
uncert = uncert(:,1:aapoi);
npN = npN(1:aapoi);
npC = npC(1:aapoi);


% make order parameters
meanct = zeros(aapoi);
for ii = 1:aapoi
    part1 = 0;
    part2 = 0;
    mean_CA_t1 = zeros(size(curr_CA_trace));
    mean_CA_t2 = zeros(size(curr_CA_trace));
    rnum = aa(ii);
    % transform to CA traces to their local frames and compute mean CA
    % trace
    for mnum = 1:models
        curr_CA_trace = CA_traces{mnum};
        coorCA = curr_CA_trace(ii,:);
        adr = mk_address([snum,cnum,mnum,rnum]);
        indN = resolve_address([adr '.N']);
        [~,coorN] = get_object(indN,'coor');
        if isempty(coorN)
            msg.error = 3;
            msg.text = sprintf('No N atom for residue %s',adr);
            return
        end
        indC = resolve_address([adr '.C']);
        [~,coorC] = get_object(indC,'coor');
        if isempty(coorC)
            msg.error = 4;
            msg.text = sprintf('No C atom for residue %s',adr);
            return
        end
        coor0 = [coorN;coorCA;coorC];
        transmat = get_trafo(coor0);
        curr_CA_trace = affine_trafo_coor(curr_CA_trace,transmat);
        CA_traces_LF{mnum} = curr_CA_trace;
        if mnum == 1
            mean_CA_trace = curr_CA_trace;
        else
            mean_CA_trace = mean_CA_trace + curr_CA_trace;
        end
        if selection(mnum) < 0.5
            mean_CA_t1 = mean_CA_t1 + curr_CA_trace;
            part1 = part1 + 1;
        else
            mean_CA_t2 = mean_CA_t2 + curr_CA_trace;
            part2 = part2 + 1;            
        end
    end
    mean_CA_trace = mean_CA_trace/models;
    mean_CA_t1 = mean_CA_t1/part1;
    mean_CA_t2 = mean_CA_t2/part2;
    % compute the ensemble average of the cos theta_n
    mean_ctheta = zeros(1,aapoi);
    mean_ctheta_p1 = zeros(1,aapoi);
    mean_ctheta_p2 = zeros(1,aapoi);    
    naxis = zeros(1,aapoi);
    for k = 1:aapoi
        naxis(k) = k - ii;
        r_mean = mean_CA_trace(k,:);
        r_mean_p1 = mean_CA_t1(k,:);
        r_mean_p2 = mean_CA_t2(k,:);
        for mnum = 1:models
            curr_CA_trace = CA_traces_LF{mnum};
            r_curr = curr_CA_trace(k,:);
            if k == ii
                ctheta = 1;
                cth_p1 = 1;
                cth_p2 = 1;
            else
                ctheta = sum(r_mean.*r_curr)/(norm(r_mean)*norm(r_curr));
                cth_p1 = sum(r_mean_p1.*r_curr)/(norm(r_mean_p1)*norm(r_curr));
                cth_p2 = sum(r_mean_p2.*r_curr)/(norm(r_mean_p2)*norm(r_curr));
            end
            mean_ctheta(k) = mean_ctheta(k) + ctheta;
            if selection(mnum) < 0.5
                mean_ctheta_p1(k) = mean_ctheta_p1(k) + cth_p1;
            else
                mean_ctheta_p2(k) = mean_ctheta_p2(k) + cth_p2;
            end
        end
    end
    mean_ctheta = mean_ctheta/models;
    mean_ctheta_p1 = mean_ctheta_p1/part1;
    mean_ctheta_p2 = mean_ctheta_p2/part2;
    meanct(ii,:) = mean_ctheta;
    s(ii) = mean(mean_ctheta);
    uncert(1,ii) = s(ii) - mean(mean_ctheta_p1);
    if uncert(1,ii) <= 0 
        uncert(1,ii) = 1e-6;
    end
    uncert(2,ii) = mean(mean_ctheta_p2) - s(ii);
    if uncert(2,ii) <= 0 
        uncert(2,ii) = 1e-6;
    end
    n = 0;
    failed = true;
    while failed && ii-n > 1
        n = n + 1;
        if mean_ctheta(ii-n) <= 1/exp(1)
            failed = false;
        end
    end
    if ~failed
        npN(ii) = n;
    end
    n = 0;
    failed = true;
    while failed && ii+n < length(mean_ctheta)
        n = n + 1;
        if mean_ctheta(ii+n) <= 1/exp(1)
            failed = false;
        end
    end
    if ~failed
        npC(ii) = n;
    end
end

for ii = 1:aapoi
    uncert(:,ii) = sort(uncert(:,ii));
end

if options.superimpose
    [maxs,ii] = max(s);
    add_msg_board(sprintf('Maximum order parameter is %5.3f at residue %i',maxs,resaxis(ii)));
end

for mnum = 1:models
    curr_CA_trace = CA_traces{mnum};
    coorCA = curr_CA_trace(ii,:);
    adr = mk_address([snum,cnum,mnum,rnum]);
    indN = resolve_address([adr '.N']);
    [~,coorN] = get_object(indN,'coor');
    if isempty(coorN)
        msg.error = 3;
        msg.text = sprintf('No N atom for residue %s',adr);
        return
    end
    indC = resolve_address([adr '.C']);
    [~,coorC] = get_object(indC,'coor');
    if isempty(coorC)
        msg.error = 4;
        msg.text = sprintf('No C atom for residue %s',adr);
        return
    end
    coor0 = [coorN;coorCA;coorC];
    transmat = get_trafo(coor0);
    coor = model.structures{snum}(cnum).xyz{mnum};
    model.structures{snum}(cnum).xyz{mnum} = affine_trafo_coor(coor,transmat);
end

function transmat = get_trafo(coor)

orig = coor(2,:);
coor = coor - repmat(orig,3,1);
x = coor(1,:)-coor(2,:); 
x = x/norm(x);    % unit vector along x
yp = coor(3,:)-coor(2,:); 
yp = yp/norm(yp);
z = cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z = z/norm(z);
y = cross_rowvec(z,x); % real (corrected) y axis
Rp = [x;y;z];
transmat1 = zeros(4);
transmat1(1:3,1:3) = Rp;
transmat1(4,4) = 1;
transmat2 = eye(4);
transmat2(1:3,4) = -orig';
transmat = transmat1*transmat2;