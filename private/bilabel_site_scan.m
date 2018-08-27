function [msg,diagnostics,sites,Cu_coor] = bilabel_site_scan(indices,label,fname,native,torsionpot)
% bilabel_site_scan(indices,label,fname,native,torsionpot)
%
% Performs a site scan of a chain model for attaching two-pronged labels,
% i.e. labels that are linked to two residues
%
% indices   chain (model) indices, [1x3] or [1x2], if only structure and
%           chain are specified, the scan extends to all models
% label     label identifier
% fname     name of the output file (.txt is appended)
% native    optional flag, if present and true, only native His pairs are
%           considered, independent of their separation in the sequence
% torsionpot  optional torsion potential (J/mol) for ligand attachment
%
% msg       error message msg.error (number), msg.text
% diagnostics diagnostic information on all sites
%
% G. Jeschke, 28.3.2018

global model
global hMain

msg.error = 0;
msg.text = 'OK';

CA_CA_threshold = 10;
max_seq_dist = 5;

if ~exist('fname','var') || isempty(fname)
    fname = 'verbose_bilabel_site_scan';
end

if ~exist('native','var') || isempty(native)
    native = false;
end

if ~exist('torsionpot','var')
    torsionpot = [];
end


fname = [fname '.txt'];

fid = fopen(fname,'wt');

if length(indices) < 2 || length(indices) > 3
    msg.error = 1;
    msg.text = 'Only chain oder chain model indices are allowed';
    return
end

if length(indices) == 2
    modind = 1:length(model.structures{indices(1)}(indices(2)).xyz);
else
    modind = indices(3);
end

diagnostics = zeros(10000,6);
if nargout > 2
    sites = cell(1,10000);
end
if nargout > 3
    Cu_coor = cell(1,10000);
end
% poi_vec = zeros(1,max_seq_dist);
% for k = 1:max_seq_dist
%     Cu_mean{k} = zeros(10000,3);
% end
poi = 0;
tic,
for kmod = modind % loop over all chain models
    info = model.structures{indices(1)}(indices(2)).residues{kmod}.info;
    [~,coor1] = get_object([indices(1:2) kmod],'xyz_heavy');
%     cofm = mean(coor1,1);
    [~,elements] = get_object([indices(1:2) kmod],'elements_heavy');
    ecoor = [elements,coor1];
    fprintf(fid,'Site1             Site2          CA-CA (Å) rotamers  part.fct. rmsd (Å) d(bb)(Å)\n');
    for kr1 = 1:length(info)-1
        if native && ~strcmpi(info(kr1).name,'HIS')
            continue
        end
        adr1 = mk_address([indices(1:2) kmod kr1]);
        [msg1,CA1] = get_object(sprintf('%s.CA',adr1),'coor');
        if msg1.error == 0
            if native
                scan_end = length(info);
            else
                scan_end = kr1 + max_seq_dist;
                if scan_end > length(info)
                    scan_end = length(info);
                end
            end
            for kr2 = kr1+1:scan_end
                if native && ~strcmpi(info(kr2).name,'HIS')
                    continue
                end
                adr2 = mk_address([indices(1:2) kmod kr2]);
                [msg1,CA2] = get_object(sprintf('%s.CA',adr2),'coor');
                if msg1.error == 0
                    dist = norm(CA1 - CA2);   
                    if dist < CA_CA_threshold
                        midpoint = (CA1+CA2)/2;
                        [popcoor,~,part_fct,mean_Cu_coor,rmsd] = two_pronged_label(adr1,adr2,label,ecoor,torsionpot);
                        [mrot,~] = size(popcoor);
                        if mrot > 0
                            dist_bb = norm(mean_Cu_coor-midpoint);
                            poi = poi + 1;
                            diagnostics(poi,1) = dist;
                            diagnostics(poi,2) = mrot;
                            diagnostics(poi,3) = part_fct;
                            diagnostics(poi,4) = rmsd;
                            diagnostics(poi,5) = dist_bb;
                            diagnostics(poi,6) = kr2 - kr1;
                            diagnostics(poi,7) = kr1;
                            diagnostics(poi,8) = kr2;
                            if nargout > 2
                                sites{poi} = [adr1 '|' adr2];
                            end
                            if nargout > 3
                                Cu_coor{poi} = popcoor(:,1:4);
                            end
                            fprintf(fid,'%16s%16s%8.1f%10i%10.4f%8.2f%10.2f\n',adr1,adr2,dist,mrot,part_fct,rmsd,dist_bb);
%                             orig = (CA1+CA2)/2;
%                             x = CA2-CA1;
%                             x = x/norm(x);
%                             yp = orig - cofm;
%                             yp = yp/norm(yp);
%                             z = cross_rowvec(x,yp);
%                             z = z/norm(z);
%                             y = cross_rowvec(z,x);
%                             Rp = [x;y;z];
%                             rel_Cu_coor = (mean_Cu_coor-orig)*Rp';
%                             seq_dist = kr2 - kr1;
%                             curr_Cu_mean = Cu_mean{seq_dist};
%                             curr_poi = poi_vec(seq_dist);
%                             curr_Cu_mean(curr_poi+1,:) = rel_Cu_coor;
%                             poi_vec(seq_dist) = poi_vec(seq_dist) + 1;
%                             Cu_mean{seq_dist} = curr_Cu_mean;
                        end
                    end
                end
            end
        end
    end
end
runtime = toc;
add_msg_board(sprintf('Bilabel site scan for %s completed in %8.1f s',label,runtime));
diagnostics = diagnostics(1:poi,:);
if nargout > 2
    sites = sites(1:poi);
end
if nargout > 3
    Cu_coor = Cu_coor(1:poi);
end
% for k = 1:max_seq_dist
%     curr_Cu_mean = Cu_mean{k};
%     curr_Cu_mean = curr_Cu_mean(1:poi_vec(k),:);
%     Cu_mean{k} = curr_Cu_mean;
% end
% save site_scan_diagnostics diagnostics Cu_mean


fclose(fid);

hMain.report_file = fname;
report_editor;
