function labeling_analysis(sites,html_name,bin_name)
% Analyzes lists of labels by spatial restrictions and by predicted pair
% distance distributions
% the types of analyzed pairs and homooligomer analysis by symmetry are
% determined by settings in hMain
% distance analysis is only performed up to 10000 pairs
%
% G. Jeschke, 2009
% St. Stoll, 2011 (speed-up of analyze_distribution.m)

global hMain
global model

max_pairs_ascii=1024; % maximum number of residue pairs in ASCII distance matrix output

% check if there is binding site and metal center information
if isfield(model,'keywords'),
    mc=tag2id('metal centers',model.keywords);
    bs=tag2id('binding sites',model.keywords);
else
    mc=[];
    bs=[];
end;

bin_output=3; % determines extent of binary output, 0 no output, 1 single-site output, 2 distance output,
                % 3 distance and single-site output

if nargin<3,
    bin_output=0;
    % bin_name=sprintf('site_scan_%s',datestr(now,'yyyy-mm-dd_HH-MM-SS'));
end;

html_output=2; % determines extent of HTML output, 0 no output, 1 only single-site, 2 also distances

partition_function_threshold=0.05; % threshold, where warning about possible labeling difficulties is displayed
mirror_threshold=1e-3; % threshold (difference of cosines) for mirror symmetry of turning angles

if html_output,
    wfile=fopen(html_name,'w');
    if wfile<0 && bin_output==0,
        add_msg_board('### ERROR ### HTML file could not be saved.');
        return
    end;
    [pathstr, name] = fileparts(html_name);
    stat_name=fullfile(pathstr,[name '_statistics.dat']);
    stat_file=fopen(stat_name,'w');
    if stat_file<0,
        add_msg_board('--- Warning --- Labeling statistics file could not be saved.');
        save_stat=0;
    else
        save_stat=1;
        fprintf(stat_file,'%%  #  \t PDB-ID \t chain \t model \t residue \t type \t rotamers  \t partition fct.  \t std. dev. (nm)\n');
    end;
    indiv_stat_name=fullfile(pathstr,[name '_individual_statistics.dat']);
    indiv_stat_file=fopen(indiv_stat_name,'w');
    if indiv_stat_file<0,
        add_msg_board('--- Warning --- Individual labeling statistics file could not be saved.');
        save_indiv_stat=0;
    else
        save_indiv_stat=1;
        fprintf(indiv_stat_file,'%%  #  \t PDB-ID \t chain \t model \t residue \t type \t rotamer  \t character    \t population\n');
    end;
    dist_name=fullfile(pathstr,[name '_distances.dat']);
end;

% Create header

if html_output,
    title='MMM spin labeling site scan';
    fprintf(wfile,'%s\n','<html>');
    fprintf(wfile,'%s\n','<!-- $Created automatically by MMM $ -->');
    fprintf(wfile,'%s\n','<head>');
    fprintf(wfile,'%s','<title>');
    fprintf(wfile,'%s',title);
    fprintf(wfile,'%s\n','</title>');
    fprintf(wfile,'%s%s%s\n','<h1>',title,'</h1>');
    fprintf(wfile,'%s\n','<body bgcolor=#fffff0>');
    fprintf(wfile,'%s\n','<P>');

    fprintf(wfile,'For missing residues no rotamer had sufficiently low energy.\n<P>\n<P>\n');
end;

chains=length(sites);

% Jump section

for kc=1:chains,
    if ~isempty(sites(kc).residue),
        mark{kc}=mk_address(sites(kc).residue(1).indices(1:3));
    end;
end;
    
if html_output,
    fprintf(wfile,'<A NAME="jump_station">\n');
    fprintf(wfile,'<h2>Overview</h2>\n');
    for kc=1:chains,
        if isempty(sites(kc).residue),
            continue;
        end;
        fprintf(wfile,'Spatial restrictions in chain model <A HREF="#spatial%i">%s </A><P>\n',kc,mark{kc});
        fprintf(wfile,'Labeling statistics for chain model <A HREF="#statistics%i">%s </A><P>\n',kc,mark{kc});
        if hMain.site_scan_intra>0 && html_output>=2,
            fprintf(wfile,'Distance analysis within chain model <A HREF="#intra%i">%s </A>\n',kc,mark{kc});
            fprintf(wfile,'; <A HREF="#intra_DEER%i">DEER </A>\n',kc);
            fprintf(wfile,'; <A HREF="#intra_CW%i">CW </A>\n',kc);
            fprintf(wfile,'; <A HREF="#intra_DEERlong%i">long (DEER)</A>\n',kc);
            fprintf(wfile,'; <A HREF="#intra_short%i">too short</A>\n',kc);
            fprintf(wfile,'; <A HREF="#intra_long%i">too long</A>\n',kc);
            fprintf(wfile,'; <A HREF="#intra_tight%i">too tight</A><P>\n',kc);
        end;
    end;
    if hMain.site_scan_inter>0 && chains>1 && html_output>=2,
        for kc1=1:chains-1,
            for kc2=kc1+1:chains,
                fprintf(wfile,'Distance analysis <A HREF="#inter%ito%i">between residues in chain models %s and %s</A>\n',kc1,kc2,mark{kc1},mark{kc2});
                fprintf(wfile,'; <A HREF="#inter_DEER%ito%i">DEER </A>\n',kc1,kc2);
                fprintf(wfile,'; <A HREF="#inter_CW%ito%i">CW </A>\n',kc1,kc2);
                fprintf(wfile,'; <A HREF="#inter_DEERlong%ito%i">DEER long </A>\n',kc1,kc2);
                fprintf(wfile,'; <A HREF="#inter_short%ito%i">too short </A>\n',kc1,kc2);
                fprintf(wfile,'; <A HREF="#inter_long%ito%i">too long</A>\n',kc1,kc2);
                fprintf(wfile,'; <A HREF="#inter_tight%ito%i">too tight</A><P>\n',kc1,kc2);
            end;
        end;
    end;
    if hMain.site_scan_homooligomer && html_output>=2,
        for kc=1:chains,
            fprintf(wfile,'Distance analysis for %i-fold homooligomer of chains <A HREF="#oligo%i">%s </A>\n',hMain.site_scan_multiplicity,kc,mark{kc});
            fprintf(wfile,'; <A HREF="#oligo_DEER%i">DEER </A>\n',kc);
            fprintf(wfile,'; <A HREF="#oligo_CW%i">CW </A>\n',kc);
            fprintf(wfile,'; <A HREF="#oligo_DEERlong%i">DEER long </A>\n',kc);
            fprintf(wfile,'; <A HREF="#oligo_short%i">too short </A>\n',kc);
            fprintf(wfile,'; <A HREF="#oligo_long%i">too long</A>\n',kc);
            fprintf(wfile,'; <A HREF="#oligo_tight%i">too tight</A><P>\n',kc);
        end;
    end;
end;

all_res=0;
for kc=1:chains,
    all_res=all_res+length(sites(kc).residue);
end;
if all_res<=max_pairs_ascii,
    dist_mat=zeros(all_res,all_res);
    dist_file=fopen(dist_name,'w');
    if dist_file<0,
        add_msg_board('### ERROR ### Distance matrix file could not be saved.');
        return
    end;
    mk_dist_matrix=1;
else
    mk_dist_matrix=0;
end;
% analysis of spatial restrictions
all_poi=0;
for kc=1:chains,
    % initialize counters for statistics
    min_rmsd_loop=1e6;
    max_rmsd_loop=0;
    mean_rmsd_loop=0;
    n_loop=0;
    min_rmsd_helix=1e6;
    max_rmsd_helix=0;
    mean_rmsd_helix=0;
    n_helix=0;
    min_rmsd_strand=1e6;
    max_rmsd_strand=0;
    mean_rmsd_strand=0;
    n_strand=0;
    if html_output,
        fprintf(wfile,'<A NAME="spatial%i">\n',kc);
        fprintf(wfile,'<h2> Spatial restrictions for sites in chain model %s</h2>\n',mark{kc});
        fprintf(wfile,'(sorted by ascending spatial distribution of NO group)<P>\n');
        fprintf(wfile,'(very tight positions, where labeling may fail, are displayed at the end)<P>\n');
    end;
    residues=length(sites(kc).residue);
    rmsd=zeros(1,residues);
    rmsdz=rmsd;
    sec=rmsd;
    partition_function=zeros(1,residues);
    num=zeros(1,residues);
    indices=zeros(residues,4);
    for kr=1:residues,
        num(kr)=length(sites(kc).residue(kr).rotamers);
        ind=sites(kc).residue(kr).indices;
        indices(kr,:)=ind;
        sec(kr)=model.structures{ind(1)}(ind(2)).residues{ind(3)}.info(ind(4)).secondary;
        rmsd(kr)=NOpos_rmsd(sites(kc).residue(kr).NOpos);
        rmsdz(kr)=NOpos_rmsdz(sites(kc).residue(kr).NOpos);
        partition_function(kr)=sites(kc).residue(kr).partition_function;
        all_poi=all_poi+1;
        if mk_dist_matrix,
            dist_mat(all_poi,all_poi)=rmsd(kr);
        end;
        if save_stat,
            [stag,ctag,modelnum,resnum,icode,atag,ltag,resname]=mk_address_parts(ind);
            fprintf(stat_file,'%4i   \t %s    \t %s    \t %2i  \t%6i %s \t %s \t %6i   \t %12.5f   \t %10.3f\n',all_poi,stag,ctag,modelnum,resnum,icode,resname,num(kr),partition_function(kr),rmsd(kr));
        end;
        if save_indiv_stat,
            [stag,ctag,modelnum,resnum,icode,atag,ltag,resname]=mk_address_parts(ind);
            character_table=rotamer_char_table(sites(kc).library);
            all_pops=sites(kc).residue(kr).NOpos(:,4);
            [sort_pops,pop_poi] = sort(all_pops,1,'descend');
            for krot=1:length(sort_pops),
                rotnum=sites(kc).residue(kr).NOpos(pop_poi(krot),5);
                character=character_table{rotnum};
                fprintf(indiv_stat_file,'%4i   \t %s    \t %s    \t %2i  \t%6i %s \t %s \t %6i   \t %12s \t %5.3f\n',all_poi,stag,ctag,modelnum,resnum,icode,resname,rotnum,character,sort_pops(krot));
            end;
        end;
    end;
    [rmsd_sort,poi]=sort(rmsd,2,'ascend');
    poi2=zeros(size(poi));
    tight=0;
    for kpoi=1:residues,
        if partition_function(poi(kpoi))<partition_function_threshold,
            tight=tight+1;
            poi2(tight)=poi(kpoi);
        else
            radr=mk_address(indices(poi(kpoi),:));
            sectype='?';
            kr=poi(kpoi);
            switch sec(poi(kpoi))
                case 0
                    sectype='loop';
                    n_loop=n_loop+1;
                    mean_rmsd_loop=mean_rmsd_loop+rmsd(kr)^2;
                    if rmsd(kr)<min_rmsd_loop, min_rmsd_loop=rmsd(kr); end;
                    if rmsd(kr)>max_rmsd_loop, max_rmsd_loop=rmsd(kr); end;
                case 1
                    sectype='helix';
                    n_helix=n_helix+1;
                    mean_rmsd_helix=mean_rmsd_helix+rmsd(kr)^2;
                    if rmsd(kr)<min_rmsd_helix, min_rmsd_helix=rmsd(kr); end;
                    if rmsd(kr)>max_rmsd_helix, max_rmsd_helix=rmsd(kr); end;
                case 2
                    sectype='strand';
                    n_strand=n_strand+1;
                    mean_rmsd_strand=mean_rmsd_strand+rmsd(kr)^2;
                    if rmsd(kr)<min_rmsd_strand, min_rmsd_strand=rmsd(kr); end;
                    if rmsd(kr)>max_rmsd_strand, max_rmsd_strand=rmsd(kr); end;
            end;
            if html_output,
                ind=sites(kc).residue(kr).indices;
                name=model.structures{ind(1)}(ind(2)).residues{ind(3)}.info(ind(4)).name;
                name(2:end)=lower(name(2:end));
                fprintf(wfile,'<P><h3>Residue %s (%s) labeled by %s in %s</h3>\n',radr,name,sites(kc).residue(kr).label,sectype);
                fprintf(wfile,'<P>NO position r.m.s.d. %4.2f nm\n',rmsd(poi(kpoi)));
                if hMain.z_analysis,
                    fprintf(wfile,'<P>NO z position r.m.s.d. %4.2f nm\n',rmsdz(poi(kpoi)));
                end;
                fprintf(wfile,'<P>Number of rotamers: %i with partition function %7.5f<P>\n',num(poi(kpoi)),partition_function(poi(kpoi)));
                if ~isempty(bs) || ~isempty(mc), % test for binding site or metal center
                    rindices=resolve_address(radr);
                    [msg,anno]=get_annotations(rindices);
                    if ~isempty(bs) && ~isempty(anno),
                        if ~isempty(find(anno.keywords==bs, 1)),
                            fprintf(wfile,'<P><b>Warning:</b> This residue is involved in a binding site<P>\n');
                        end;
                    end;
                    if ~isempty(mc),
                        [msg,aindices]=get_object(radr,'children');
                        [maa,naa]=size(aindices);
                        for kaa=1:maa,
                            [msg,anno]=get_annotations(aindices(kaa,:));
                            if ~isempty(anno),
                                if ~isempty(find(anno.keywords==mc, 1)),
                                    aadr=mk_address(aindices(kaa,:));
                                    paa=findstr(aadr,'.');
                                    if ~isempty(paa) && paa<length(aadr),
                                        atom_tag=aadr(paa+1:end);
                                        fprintf(wfile,'<P><b>Warning:</b> Atom %s of this residue is involved in metal coordination<P>\n',atom_tag);
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    if n_loop>0,
        mean_rmsd_loop=sqrt(mean_rmsd_loop/n_loop);
    end;
    if n_helix>0,
        mean_rmsd_helix=sqrt(mean_rmsd_helix/n_helix);
    end;
    if n_strand>0,
        mean_rmsd_strand=sqrt(mean_rmsd_strand/n_strand);
    end;
    for kpoi=1:tight,
        radr=mk_address(indices(poi2(kpoi),:));
        sectype='?';
        switch sec(poi2(kpoi))
            case 0
                sectype='loop';
            case 1
                sectype='helix';
            case 2
                sectype='strand';
        end;
        if html_output,
            ind=sites(kc).residue(poi2(kpoi)).indices;
            name=model.structures{ind(1)}(ind(2)).residues{ind(3)}.info(ind(4)).name;
            name(2:end)=lower(name(2:end));
            fprintf(wfile,'<P><h3>Residue %s (%s) labeled by %s in %s</h3>\n',radr,name,sites(kc).residue(kr).label,sectype);
            fprintf(wfile,'<P>NO position r.m.s.d. %4.2f nm\n',rmsd(poi2(kpoi)));
            fprintf(wfile,'<P>Number of rotamers: %i with partition function %7.5f<P>\n',num(poi2(kpoi)),partition_function(poi2(kpoi)));
            fprintf(wfile,'Be careful. Labeling might fail here or distort structure.<P>\n');
            if ~isempty(bs) || ~isempty(mc), % test for binding site or metal center
                rindices=resolve_address(radr);
                [msg,anno]=get_annotations(rindices);
                if ~isempty(bs) && ~isempty(anno),
                    if ~isempty(find(anno.keywords==bs, 1)),
                        fprintf(wfile,'<P><b>Warning:</b> This residue is involved in a binding site<P>\n');
                    end;
                end;
                if ~isempty(mc),
                    [msg,aindices]=get_object(radr,'children');
                    [maa,naa]=size(aindices);
                    for kaa=1:maa,
                        [msg,anno]=get_annotations(aindices(kaa,:));
                        if ~isempty(anno),
                            if ~isempty(find(anno.keywords==mc, 1)),
                                aadr=mk_address(aindices(kaa,:));
                                paa=findstr(aadr,'.');
                                if ~isempty(paa) && paa<length(aadr),
                                    atom_tag=aadr(paa+1:end);
                                    fprintf(wfile,'<P><b>Warning:</b> Atom %s of this residue is involved in metal coordination<P>\n',atom_tag);
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    if html_output,
        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
        fprintf(wfile,'<A NAME="statistics%i">\n',kc);
        fprintf(wfile,'<h2> Labeling statistics for sites in chain model %s</h2>\n',mark{kc});
        fprintf(wfile,'<P>In %i residues there are %i tight sites where labeling might fail.\n',residues,tight);
        if n_loop>0,
            fprintf(wfile,'<h3>Loop sites</h3>\n');
            fprintf(wfile,'<p>Mean    r.m.s.d: %4.2f nm\n',mean_rmsd_loop);
            fprintf(wfile,'<p>Minimum r.m.s.d: %4.2f nm\n',min_rmsd_loop);
            fprintf(wfile,'<p>Maximum r.m.s.d: %4.2f nm\n',max_rmsd_loop);
        end;
        if n_helix>0,
            fprintf(wfile,'<h3>Helix sites</h3>\n');
            fprintf(wfile,'<p>Mean    r.m.s.d: %4.2f nm\n',mean_rmsd_helix);
            fprintf(wfile,'<p>Minimum r.m.s.d: %4.2f nm\n',min_rmsd_helix);
            fprintf(wfile,'<p>Maximum r.m.s.d: %4.2f nm\n',max_rmsd_helix);
        end;
        if n_strand>0,
            fprintf(wfile,'<h3>Strand sites</h3>\n');
            fprintf(wfile,'<p>Mean    r.m.s.d: %4.2f nm\n',mean_rmsd_strand);
            fprintf(wfile,'<p>Minimum r.m.s.d: %4.2f nm\n',min_rmsd_strand);
            fprintf(wfile,'<p>Maximum r.m.s.d: %4.2f nm\n<p>\n',max_rmsd_strand);
        end;
        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
    end;
    if mod(bin_output,2)==1,
        outname=sprintf('%s_%s.mat',bin_name,mark{kc});
        save(outname,'num','indices','sec','rmsd','rmsdz','partition_function');
    end;
end;

if save_stat,
    fclose(stat_file);
end;

if save_indiv_stat,
    fclose(indiv_stat_file);
end;

% analysis of intrachain distances (if requested)
if hMain.site_scan_intra>0,
    for kc=1:chains,
        residues=length(sites(kc).residue);
        if residues>=2,
            if html_output>=2,
                fprintf(wfile,'<A NAME="intra%i">\n',kc);
                fprintf(wfile,'<h2>Intrachain distances in chain model %s</h2>\n',mark{kc});
                fprintf(wfile,'(sorted in categories by ascending relative distribution width)<P>\n');
            end;
            pairs=residues*(residues-1)/2; % number of spin pairs
            if pairs<=15000,
                r_mean=zeros(1,pairs);
                r_stddev=zeros(1,pairs);
                ind1=zeros(pairs,4);
                ind2=zeros(pairs,4);
                pfmin=zeros(pairs,1);
                poi=0;
                for kr1=1:residues-1,
                    NOpos1=sites(kc).residue(kr1).NOpos;
                    for kr2=kr1+1:residues,
                        poi=poi+1;
                        ind1(poi,:)=sites(kc).residue(kr1).indices;
                        ind2(poi,:)=sites(kc).residue(kr2).indices;
                        pf1=sites(kc).residue(kr1).partition_function;
                        pf2=sites(kc).residue(kr2).partition_function;
                        pfmin(poi)=min([pf1,pf2]);
                        NOpos2=sites(kc).residue(kr2).NOpos;
                        [rm,sr]=analyze_distribution(NOpos1,NOpos2);
                        r_mean(poi)=rm;
                        r_stddev(poi)=sr;
                    end;
                end;
                if bin_output>=2,
                    outname=sprintf('%s_intra_%s.mat',bin_name,mark{kc});
                    save(outname,'pairs','ind1','ind2','pfmin','r_mean','r_stddev');
                end;
                [rel_width,poi]=sort(r_stddev./r_mean,2,'ascend');
                category0=zeros(1,pairs);
                cpoi0=0;
                type0='At least one site probably too tight for labeling';
                category1=zeros(1,pairs);
                cpoi1=0;
                type1='Distance too short for proper measurement (<r> < 8 Å)';
                category2=zeros(1,pairs);
                cpoi2=0;
                type2='Most favorable for CW EPR measurements (0.8 nm <= <r> < 1.8 nm)';
                category3=zeros(1,pairs);
                cpoi3=0;
                type3='Most favorable for DEER measurements (1.8 nm <= <r> < 4.5 nm)';
                category4=zeros(1,pairs);
                cpoi4=0;
                type4='DEER measurement possible, but rather imprecise (4.5 nm <= <r> < 7.0 nm)';
                category5=zeros(1,pairs);
                cpoi5=0;
                type5='Distance probably too long for useful measurement (<r> = 7.0 nm)';
                for kpoi=1:pairs,
                    if pfmin(poi(kpoi))<partition_function_threshold,
                        cpoi0=cpoi0+1;
                        category0(cpoi0)=poi(kpoi);
                    elseif r_mean(poi(kpoi))<0.8,
                        cpoi1=cpoi1+1;
                        category1(cpoi1)=poi(kpoi);
                    elseif r_mean(poi(kpoi))>=0.8 && r_mean(poi(kpoi))<1.8
                        cpoi2=cpoi2+1;
                        category2(cpoi2)=poi(kpoi);
                    elseif r_mean(poi(kpoi))>=1.8 && r_mean(poi(kpoi))<4.5
                        cpoi3=cpoi3+1;
                        category3(cpoi3)=poi(kpoi);
                    elseif r_mean(poi(kpoi))>=4.5 && r_mean(poi(kpoi))<7.0
                        cpoi4=cpoi4+1;
                        category4(cpoi4)=poi(kpoi);
                    elseif r_mean(poi(kpoi))>=7.0
                        cpoi5=cpoi5+1;
                        category5(cpoi5)=poi(kpoi);
                    end;    
                end;
                category0=category0(1:cpoi0);
                category1=category1(1:cpoi1);
                category2=category2(1:cpoi2);
                category3=category3(1:cpoi3);
                category4=category4(1:cpoi4);
                category5=category5(1:cpoi5);
                if html_output>=2,
                    if ~isempty(category3),
                        fprintf(wfile,'<A NAME="intra_DEER%i">\n',kc);
                        fprintf(wfile,'<h3>%s</h3>\n',type3);
                        list_distances(wfile,ind1,ind2,r_mean,r_stddev,category3);
                        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                    end;
                    if ~isempty(category2),
                        fprintf(wfile,'<A NAME="intra_CW%i">\n',kc);
                        fprintf(wfile,'<h3>%s</h3>\n',type2);
                        list_distances(wfile,ind1,ind2,r_mean,r_stddev,category2);
                        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                    end;
                    if ~isempty(category4),
                        fprintf(wfile,'<A NAME="intra_DEERlong%i">\n',kc);
                        fprintf(wfile,'<h3>%s</h3>\n',type4);
                        list_distances(wfile,ind1,ind2,r_mean,r_stddev,category4);
                        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                    end;
                    if ~isempty(category1),
                        fprintf(wfile,'<A NAME="intra_short%i">\n',kc);
                        fprintf(wfile,'<h3>%s</h3>\n',type1);
                        list_distances(wfile,ind1,ind2,r_mean,r_stddev,category1);
                        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                    end;
                    if ~isempty(category5),
                        fprintf(wfile,'<A NAME="intra_long%i">\n',kc);
                        fprintf(wfile,'<h3>%s</h3>\n',type5);
                        list_distances(wfile,ind1,ind2,r_mean,r_stddev,category5);
                        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                    end;
                    if ~isempty(category0),
                        fprintf(wfile,'<A NAME="intra_tight%i">\n',kc);
                        fprintf(wfile,'<h3>%s</h3>\n',type0);
                        list_distances(wfile,ind1,ind2,r_mean,r_stddev,category0);
                        fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                    end;
                end;
            else
                if html_output>=2,
                    fprintf(wfile,'<p>More than 10000 pairs of residues in chain %s. Distance analysis skipped.<P>\n',mark{kc});
                end;
            end;
        else
            if html_output>=2,
                fprintf(wfile,'<p>Less than two residues in chain model %s. No distance analysis.<P>\n',mark{kc});
            end;
        end;
    end;
end;

if hMain.site_scan_inter>0 && chains>1,
    for kc1=1:chains-1,
        residues1=length(sites(kc1).residue);
        for kc2=kc1+1:chains,
            residues2=length(sites(kc2).residue);
            pairs=residues1*residues2; % number of spin pairs
            if pairs<=15000,
                if pairs>0,
                    if html_output>=2,
                        fprintf(wfile,'<A NAME="inter%ito%i">\n',kc1,kc2);
                        mode='sites';
                        if hMain.site_scan_inter==1, mode='equivalent sites'; end;
                        fprintf(wfile,'<h2>Distances <A NAME="inter%ito%i">between %s in chain models %s and %s</A></h2>\n',kc1,kc2,mode,mark{kc1},mark{kc2});
                        fprintf(wfile,'(sorted in categories by ascending relative distribution width)<P>\n');
                    end;
                    r_mean=zeros(1,pairs);
                    r_stddev=zeros(1,pairs);
                    ind1=zeros(pairs,4);
                    ind2=zeros(pairs,4);
                    pfmin=zeros(1,pairs);
                    poi=0;
                    for kr1=1:residues1,
                        NOpos1=sites(kc1).residue(kr1).NOpos;
                        for kr2=1:residues2,
                            if sites(kc1).residue(kr1).indices(4)==sites(kc2).residue(kr2).indices(4) || hMain.site_scan_inter>1,
                                poi=poi+1;
                                ind1(poi,:)=sites(kc1).residue(kr1).indices;
                                ind2(poi,:)=sites(kc2).residue(kr2).indices;
                                pf1=sites(kc1).residue(kr1).partition_function;
                                pf2=sites(kc2).residue(kr2).partition_function;
                                pfmin(poi)=min([pf1,pf2]);
                                NOpos2=sites(kc2).residue(kr2).NOpos;
                                [rm,sr]=analyze_distribution(NOpos1,NOpos2);
                                r_mean(poi)=rm;
                                r_stddev(poi)=sr;
                            end;
                        end;
                    end;
                    pairs=poi;
                    ind1=ind1(1:pairs,:);
                    ind2=ind2(1:pairs,:);
                    r_mean=r_mean(1:pairs);
                    r_stddev=r_stddev(1:pairs);
                    if bin_output>=2,
                        outname=sprintf('%s_inter_%s_%s.mat',bin_name,mark{kc1},mark{kc2});
                        save(outname,'pairs','ind1','ind2','pfmin','r_mean','r_stddev');
                    end;
                    [rel_width,poi]=sort(r_stddev./r_mean,2,'ascend');
                    category0=zeros(1,pairs);
                    cpoi0=0;
                    type0='At least one site probably too tight for labeling';
                    category1=zeros(1,pairs);
                    cpoi1=0;
                    type1='Distance too short for proper measurement (<r> < 8 Å)';
                    category2=zeros(1,pairs);
                    cpoi2=0;
                    type2='Most favorable for CW EPR measurements (0.8 nm <= <r> < 1.8 nm)';
                    category3=zeros(1,pairs);
                    cpoi3=0;
                    type3='Most favorable for DEER measurements (1.8 nm <= <r> < 4.5 nm)';
                    category4=zeros(1,pairs);
                    cpoi4=0;
                    type4='DEER measurement possible, but rather imprecise (4.5 nm <= <r> < 7.0 nm)';
                    category5=zeros(1,pairs);
                    cpoi5=0;
                    type5='Distance probably too long for useful measurement (<r> = 7.0 nm)';
                    for kpoi=1:pairs,
                        if pfmin(poi(kpoi))<partition_function_threshold,
                            cpoi0=cpoi0+1;
                            category0(cpoi0)=poi(kpoi);
                        elseif r_mean(poi(kpoi))<0.8,
                            cpoi1=cpoi1+1;
                            category1(cpoi1)=poi(kpoi);
                        elseif r_mean(poi(kpoi))>=0.8 && r_mean(poi(kpoi))<1.8
                            cpoi2=cpoi2+1;
                            category2(cpoi2)=poi(kpoi);
                        elseif r_mean(poi(kpoi))>=1.8 && r_mean(poi(kpoi))<4.5
                            cpoi3=cpoi3+1;
                            category3(cpoi3)=poi(kpoi);
                        elseif r_mean(poi(kpoi))>=4.5 && r_mean(poi(kpoi))<7.0
                            cpoi4=cpoi4+1;
                            category4(cpoi4)=poi(kpoi);
                        elseif r_mean(poi(kpoi))>=7.0
                            cpoi5=cpoi5+1;
                            category5(cpoi5)=poi(kpoi);
                        end;    
                    end;
                    category0=category0(1:cpoi0);
                    category1=category1(1:cpoi1);
                    category2=category2(1:cpoi2);
                    category3=category3(1:cpoi3);
                    category4=category4(1:cpoi4);
                    category5=category5(1:cpoi5);
                    if html_output>=2,
                        if ~isempty(category3),
                            fprintf(wfile,'<A NAME="inter_DEER%ito%i">\n',kc1,kc2);
                            fprintf(wfile,'<h3>%s</h3>\n',type3);
                            list_distances(wfile,ind1,ind2,r_mean,r_stddev,category3);
                            fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                        end;
                        if ~isempty(category2),
                            fprintf(wfile,'<A NAME="inter_CW%ito%i">\n',kc1,kc2);
                            fprintf(wfile,'<h3>%s</h3>\n',type2);
                            list_distances(wfile,ind1,ind2,r_mean,r_stddev,category2);
                            fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                        end;
                        if ~isempty(category4),
                            fprintf(wfile,'<A NAME="inter_DEERlong%ito%i">\n',kc1,kc2);
                            fprintf(wfile,'<h3>%s</h3>\n',type4);
                            list_distances(wfile,ind1,ind2,r_mean,r_stddev,category4);
                            fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                        end;
                        if ~isempty(category1),
                            fprintf(wfile,'<A NAME="inter_short%ito%i">\n',kc1,kc2);
                            fprintf(wfile,'<h3>%s</h3>\n',type1);
                            list_distances(wfile,ind1,ind2,r_mean,r_stddev,category1);
                            fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                        end;
                        if ~isempty(category5),
                            fprintf(wfile,'<A NAME="inter_long%ito%i">\n',kc1,kc2);
                            fprintf(wfile,'<h3>%s</h3>\n',type5);
                            list_distances(wfile,ind1,ind2,r_mean,r_stddev,category5);
                            fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                        end;
                        if ~isempty(category0),
                            fprintf(wfile,'<A NAME="inter_tight%ito%i">\n',kc1,kc2);
                            fprintf(wfile,'<h3>%s</h3>\n',type0);
                            list_distances(wfile,ind1,ind2,r_mean,r_stddev,category0);
                            fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                        end;
                    end;
                else
                    if html_output>=2,
                        fprintf(wfile,'No site pair combines chain models %s and %s. No distance analysis.<P>\n',mark{kc1},mark{kc2});
                    end;
                end;
            else
                if html_output>=2,
                    fprintf(wfile,'More than 10000 pairs of residues  between chain %s and %s. Distance analysis skipped.<P>\n',mark{kc1},mark{kc2});
                end;
            end;
        end;
    end;
end;

if hMain.site_scan_homooligomer && hMain.site_scan_multiplicity>1,
    mult=hMain.site_scan_multiplicity;
    % make set of turning angles with corresponding populations
    % consider mirror symmetry of distances
    ang=zeros(1,mult);
    cang=0;
    apop=zeros(1,mult);
    apop(1)=1;
    poi=1;
    dang=2*pi/mult;
    for k=2:mult,
        cang=cang+dang;
        known=0;
        for kk=1:poi,
            if abs(cos(cang)-cos(ang(kk)))< mirror_threshold, known=1; apop(kk)=apop(kk)+1; break; end;
        end;
        if ~known, poi=poi+1; ang(poi)=cang; apop(poi)=1; end;
    end;
    ang=ang(1:poi);
    apop=apop(1:poi);
    for kc=1:chains,
        residues=length(sites(kc).residue);
        if residues>=1,
            if html_output>=2,
                fprintf(wfile,'<A NAME="oligo%i">\n',kc);
                fprintf(wfile,'<h2>Homooligomer distances for chain model %s with multiplicity %i</h2>\n',mark{kc},hMain.site_scan_multiplicity);
                fprintf(wfile,'(sorted in categories by ascending relative distribution width)<P>\n');
            end;
            pairs=residues*(length(ang)-1); % number of spin pairs
            r_mean=zeros(1,pairs);
            r_stddev=zeros(1,pairs);
            angles=zeros(1,pairs);
            populations=zeros(1,pairs);
            ind1=zeros(pairs,4);
            pfmin=zeros(1,pairs);
            poi=0;
            for kr1=1:residues,
                NOpos1=sites(kc).residue(kr1).NOpos;
                [mm,nn]=size(NOpos1);
                for ka=2:length(ang),
                    NOpos2=NOpos1;
                    transmat=affine('rotz',ang(ka));
                    for kk=1:mm,
                        NOpos2(kk,1:3)=affine_trafo_point(NOpos1(kk,1:3),transmat);
                    end;
                    poi=poi+1;
                    angles(poi)=180*ang(ka)/pi;
                    populations(poi)=apop(ka);
                    ind1(poi,:)=sites(kc).residue(kr1).indices;
                    pfmin(poi)=sites(kc).residue(kr1).partition_function;
                    [rm,sr]=analyze_distribution(NOpos1,NOpos2);
                    r_mean(poi)=rm;
                    r_stddev(poi)=sr;
                end;
            end;
            if bin_output>=2,
                outname=sprintf('%s_oligo_%s.mat',bin_name,mark{kc});
                save(outname,'pairs','ind1','ind2','pfmin','r_mean','r_stddev');
            end;
            [rel_width,poi]=sort(r_stddev./r_mean,2,'ascend');
            category0=zeros(1,pairs);
            cpoi0=0;
            type0='Site probably too tight for labeling';
            category1=zeros(1,pairs);
            cpoi1=0;
            type1='Distance too short for proper measurement (<r> < 8 Å)';
            category2=zeros(1,pairs);
            cpoi2=0;
            type2='Most favorable for CW EPR measurements (0.8 nm <= <r> < 1.8 nm)';
            category3=zeros(1,pairs);
            cpoi3=0;
            type3='Most favorable for DEER measurements (1.8 nm <= <r> < 4.5 nm)';
            category4=zeros(1,pairs);
            cpoi4=0;
            type4='DEER measurement possible, but rather imprecise (4.5 nm <= <r> < 7.0 nm)';
            category5=zeros(1,pairs);
            cpoi5=0;
            type5='Distance probably too long for useful measurement (<r> = 7.0 nm)';
            for kpoi=1:pairs,
                if pfmin(poi(kpoi))<partition_function_threshold,
                    cpoi0=cpoi0+1;
                    category0(cpoi0)=poi(kpoi);
                elseif r_mean(poi(kpoi))<0.8,
                    cpoi1=cpoi1+1;
                    category1(cpoi1)=poi(kpoi);
                elseif r_mean(poi(kpoi))>=0.8 && r_mean(poi(kpoi))<1.8
                    cpoi2=cpoi2+1;
                    category2(cpoi2)=poi(kpoi);
                elseif r_mean(poi(kpoi))>=1.8 && r_mean(poi(kpoi))<4.5
                    cpoi3=cpoi3+1;
                    category3(cpoi3)=poi(kpoi);
                elseif r_mean(poi(kpoi))>=4.5 && r_mean(poi(kpoi))<7.0
                    cpoi4=cpoi4+1;
                    category4(cpoi4)=poi(kpoi);
                elseif r_mean(poi(kpoi))>=7.0
                    cpoi5=cpoi5+1;
                    category5(cpoi5)=poi(kpoi);
                end;    
            end;
            category0=category0(1:cpoi0);
            category1=category1(1:cpoi1);
            category2=category2(1:cpoi2);
            category3=category3(1:cpoi3);
            category4=category4(1:cpoi4);
            category5=category5(1:cpoi5);
            if html_output>=2,
                if ~isempty(category3),
                    fprintf(wfile,'<A NAME="oligo_DEER%i">\n',kc);
                    fprintf(wfile,'<h3>%s</h3>\n',type3);
                    list_distances(wfile,ind1,ind1,r_mean,r_stddev,category3,angles,populations);
                    fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                end;
                if ~isempty(category2),
                    fprintf(wfile,'<A NAME="oligo_CW%i">\n',kc);
                    fprintf(wfile,'<h3>%s</h3>\n',type2);
                    list_distances(wfile,ind1,ind1,r_mean,r_stddev,category2,angles,populations);
                    fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                end;
                if ~isempty(category4),
                    fprintf(wfile,'<A NAME="oligo_DEERlong%i">\n',kc);
                    fprintf(wfile,'<h3>%s</h3>\n',type4);
                    list_distances(wfile,ind1,ind1,r_mean,r_stddev,category4,angles,populations);
                    fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                end;
                if ~isempty(category1),
                    fprintf(wfile,'<A NAME="oligo_short%i">\n',kc);
                    fprintf(wfile,'<h3>%s</h3>\n',type1);
                    list_distances(wfile,ind1,ind1,r_mean,r_stddev,category1,angles,populations);
                    fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                end;
                if ~isempty(category5),
                    fprintf(wfile,'<A NAME="oligo_long%i">\n',kc);
                    fprintf(wfile,'<h3>%s</h3>\n',type5);
                    list_distances(wfile,ind1,ind1,r_mean,r_stddev,category5,angles,populations);
                    fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                end;            
                if ~isempty(category0),
                    fprintf(wfile,'<A NAME="oligo_tight%i">\n',kc);
                    fprintf(wfile,'<h3>%s</h3>\n',type0);
                    list_distances(wfile,ind1,ind1,r_mean,r_stddev,category0,angles,populations);
                    fprintf(wfile,'<A HREF="#jump_station">To overview </A>\n');
                end;
            end;
        end;
    end;
end;

fclose(wfile);

if mk_dist_matrix,
    poi1=0;
    for kc1=1:chains,
        residues1=length(sites(kc1).residue);
        for kr1=1:residues1,
            NOpos1=sites(kc1).residue(kr1).NOpos;
            poi1=poi1+1;
            poi2=0;
            for kc2=1:chains,
                residues2=length(sites(kc2).residue);
                for kr2=1:residues2,
                    poi2=poi2+1;
                    if poi2>poi1,
                        NOpos2=sites(kc2).residue(kr2).NOpos;
                        [rm,sr]=analyze_distribution(NOpos1,NOpos2);
                        dist_mat(poi1,poi2)=rm;
                        dist_mat(poi2,poi1)=sr;
                    end;
                end;
            end;
        end;
    end;
    for k1=1:all_res,
        for k2=1:all_res,
            fprintf(dist_file,'%6.3f\t',dist_mat(k1,k2));
        end;
        fprintf(dist_file,'\n');
    end;
    fclose(dist_file);
end;
                
function rmsd=NOpos_rmsd(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
xmean=sum(NOall(:,1).*pop);
ymean=sum(NOall(:,2).*pop);
zmean=sum(NOall(:,3).*pop);
dx=(NOall(:,1)-xmean);
dy=(NOall(:,2)-ymean);
dz=(NOall(:,3)-zmean);
nNO=length(dx);
rmsd=sqrt(0.005+nNO*sum(dx.^2.*pop+dy.^2.*pop+dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function rmsd=NOpos_rmsdz(NOall)
% in nm(!)

pop=NOall(:,4);
pop=pop/sum(pop);
zmean=sum(NOall(:,3).*pop);
dz=(NOall(:,3)-zmean);
nNO=length(dz);
rmsd=sqrt(0.005/3+nNO*sum(dz.^2.*pop)/(nNO-1))/10; % divided by 10 for Å -> nm

function list_distances(wfile,ind1,ind2,r_mean,r_stddev,category,angles,apop)

global model

rel_width=r_stddev./r_mean;

for k=1:length(category),
    adr1=mk_address(ind1(category(k),:));
    adr2=mk_address(ind2(category(k),:));
    ind1a=ind1(category(k),:);
    ind2a=ind2(category(k),:);
    name1=model.structures{ind1a(1)}(ind1a(2)).residues{ind1a(3)}.info(ind1a(4)).name;
    name1(2:end)=lower(name1(2:end));
    name2=model.structures{ind2a(1)}(ind2a(2)).residues{ind2a(3)}.info(ind2a(4)).name;
    name2(2:end)=lower(name2(2:end));
    if nargin<7,
        fprintf(wfile,'Distance between residues %s (%s) and %s (%s) is %4.2f nm<P>\n',adr1,name1,adr2,name2,r_mean(category(k)));
    else
        fprintf(wfile,'Distance between residues %s (%s) at rotation angle %4.1f° occurs %i times and is %4.2f nm<P>\n',adr1,name1,angles(category(k)),apop(category(k)),r_mean(category(k)));
    end;
    fprintf(wfile,'standard deviation is %4.2f nm and relative width %3.1f %%<P>\n',r_stddev(category(k)),100*rel_width(category(k)));
end;

function [rm,sr]=analyze_distribution(NOpos1,NOpos2)

pop1=NOpos1(:,4);
n1=length(pop1);
xyz1=NOpos1(:,1:3)/10; % divided by 10 for Å -> nm
pop2=NOpos2(:,4);
n2=length(pop2);
xyz2=NOpos2(:,1:3)/10; % divided by 10 for Å -> nm

xyz1sq = repmat(sum(xyz1.^2,2),1,n2);
xyz2sq = repmat(sum(xyz2.^2,2),1,n1).';
rij = sqrt(abs(xyz1sq + xyz2sq - 2*xyz1*xyz2.'));
rij = rij(:);
popij = pop1*pop2.';
popij = popij(:);

rm=sum(rij.*popij)/sum(popij);
dr=rij-rm;
sr=sqrt(0.01+n1*n2/(sum(popij)*(n1*n2-1))*sum(popij.*dr.^2));
