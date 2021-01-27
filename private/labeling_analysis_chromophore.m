function labeling_analysis_chromophore(sites,html_name,bin_name)
% Analyses lists of labels by spatial restrictions and by predicted pair
% distance distributions.
% The types of analyzed pairs and homooligomer analysis by symmetry are
% determined by settings in hMain
% Distance analysis is only performed up to 15000 pairs.

global hMain
global model

% check if there is binding site and metal center information
if isfield(model,'keywords')
  MetalCenter = tag2id('metal centers',model.keywords);
  BindingSite = tag2id('binding sites',model.keywords);
else
  MetalCenter = [];
  BindingSite = [];
end

Options.matOutput=3; % determines extent of binary output, 0 no output, 1 single-site output, 2 distance output,
% 3 distance and single-site output

if nargin<3,
  Options.matOutput=0;
  % bin_name=sprintf('site_scan_%s',datestr(now,'yyyy-mm-dd_HH-MM-SS'));
end;

Options.HtmlOutputLevel = 2; % determines extent of HTML output, 0 no output, 1 only single-site, 2 also distances
Options.max_pairs_ascii = 1024; % maximum number of residue pairs in ASCII distance matrix output
Options.mirror_threshold=1e-3; % threshold (difference of cosines) for mirror symmetry of turning angles
Options.maxNumberOfPairsHtml = 15000;

cssStyle = sprintf([...
  '  body {background-color:#ffffff;}\n'...
  '  td {padding-left:3px;padding-right:3px;padding-top:0px;padding-bottom:0px;margin:0;}\n'...
  '  h4 {margin-top:6px;margin-bottom:2px;}\n'...
  '  .pdb {color:#999;}\n'...
  '  .chain {color:#070;}\n'...
  '  .model {font-size:80%%;}\n'...
  '  .distcategory {font-weight:bold;}\n'...
    ]);


Options.DistanceOutput = (Options.HtmlOutputLevel>=2);

if (Options.HtmlOutputLevel>0)
  wfile = fopen(html_name,'w');
  if (wfile<0) && (Options.matOutput==0)
    add_msg_board('### ERROR ### HTML file could not be saved.');
    return
  end
  fw = @(varargin)fprintf(wfile,varargin{:});
  [pathstr, name] = fileparts(html_name);
  stat_name = fullfile(pathstr,[name '_statistics.dat']);
  stat_file = fopen(stat_name,'w');
  if (stat_file<0)
    add_msg_board('--- Warning --- Labeling statistics file could not be saved.');
    saveStatistics=0;
  else
    saveStatistics=1;
    fprintf(stat_file,'%%  #  \t PDB-ID \t chain \t model \t residue \t type \t rotamers  \t partition fct.  \t std. dev. (nm)\n');
  end
  indiv_stat_name=fullfile(pathstr,[name '_individual_statistics.dat']);
  indiv_stat_file=fopen(indiv_stat_name,'w');
  if indiv_stat_file<0,
    add_msg_board('--- Warning --- Individual labeling statistics file could not be saved.');
    saveIndivStatistics=0;
  else
    saveIndivStatistics=1;
    fprintf(indiv_stat_file,'%%  #  \t PDB-ID \t chain \t model \t residue \t type \t rotamer  \t character    \t population\n');
  end
  dist_name=fullfile(pathstr,[name '_distances.dat']);
end

  % Header, title, style sheet
if (Options.HtmlOutputLevel>0)
  fw('<html>\n');
  fw('<!-- $Created automatically by MMM $ -->\n');
  fw('<head>\n');
  fw('<title>MMM fluorescence labeling analysis</title>\n');
  fw('<style type="text/css">\n%s</style>',cssStyle);
  fw('</head>\n');
  fw('<body>\n');
  fw('<h1>MMM fluorescence labeling site scan</h1>\n');
  fw('<p>Sites without rotamers of sufficiently low energy are omitted.\n<P>\n<P>\n');
end

nChains = length(sites);
for kc = nChains:-1:1
  nResidues(kc) = length(sites(kc).residue);
  ChainName{kc} = mk_address(sites(kc).residue(1).indices(1:3));
  ChainNameHtml{kc} = makeHtmlChainName(ChainName{kc});
end

% Prepare distance matrix file
%-----------------------------------------------------------------
nTotalResidues = sum(nResidues);
if (nTotalResidues<=Options.max_pairs_ascii)
  dist_file = fopen(dist_name,'w');
  if (dist_file<0)
    add_msg_board('### ERROR ### Distance matrix file could not be saved.');
    saveDistanceMatrix = 0;
  else
    dist_mat = zeros(nTotalResidues,nTotalResidues);
    saveDistanceMatrix = 1;
  end
else
  saveDistanceMatrix = 0;
end

% Overview section
%---------------------------------------------------------------------------
if (Options.HtmlOutputLevel>=1)
  fw('<a name="overview">\n');
  fw('<h2>Overview</h2>\n');  
  fw('<p>Site analysis:<br>');
  for kc = 1:nChains
    if nResidues(kc)==0, continue; end
    fw('&nbsp;&nbsp;Chain %s: <a href="#site%i">%i sites</a> ',ChainNameHtml{kc},kc,nResidues(kc));
  end
end

if Options.DistanceOutput
  fw('<p>Distance analysis (intra-chain):<br>');
  if (hMain.site_scan_intra>0)
    for kc = 1:nChains
      fw('&nbsp;&nbsp;Chain %s: ',ChainNameHtml{kc});
      fw('<A HREF="#intra%i_FRET">FRET</A> | ',kc);
%       fw('<A HREF="#intra%i_CW">CW</A> | ',kc);
      fw('<A HREF="#intra%i_FRETlong">long (FRET)</A> | ',kc);
      fw('<A HREF="#intra%i_short">too short</A> | ',kc);
      fw('<A HREF="#intra%i_long">too long</A> | ',kc);
      fw('<A HREF="#intra%i_tight">too tight</A><br>\n',kc);
    end
  end
  if (hMain.site_scan_inter>0) && (nChains>1)
    fw('<p>Distance analysis (inter-chain):<br>');
    for kc1 = 1:nChains-1
      for kc2 = kc1+1:nChains
        fw('&nbsp;&nbsp;<A HREF="#inter%ito%i"></a>Chains %s and %s: ',kc1,kc2,ChainNameHtml{kc1},ChainNameHtml{kc2});
        fw('<A HREF="#inter%ito%i_FRET">FRET</A> | ',kc1,kc2);
%         fw('<A HREF="#inter%ito%i_CW">CW</A> | ',kc1,kc2);
        fw('<A HREF="#inter%ito%i_FRETlong">FRET long</A> | ',kc1,kc2);
        fw('<A HREF="#inter%ito%i_short">too short</A> | ',kc1,kc2);
        fw('<A HREF="#inter%ito%i_long">too long</A> | ',kc1,kc2);
        fw('<A HREF="#inter%ito%i_tight">too tight</A><br>\n',kc1,kc2);
      end
    end
  end
  if (hMain.site_scan_homooligomer)
    for kc=1:nChains
      fw('Distance analysis for %i-fold homooligomer of chains <A HREF="#oligo%i">%s </A>\n',hMain.site_scan_multiplicity,kc,ChainName{kc});
      fw('; <A HREF="#oligo%i_FRET">FRET</A> | ',kc);
%       fw('; <A HREF="#oligo%i_CW">CW</A> | ',kc);
      fw('; <A HREF="#oligo%i_FRETlong">FRET long</A> | ',kc);
      fw('; <A HREF="#oligo%i_short">too short</A> | ',kc);
      fw('; <A HREF="#oligo%i_long">too long</A> | ',kc);
      fw('; <A HREF="#oligo%i_tight">too tight</A>\n',kc);
    end
  end
end

% Analysis of spatial restrictions
%--------------------------------------------------------------------------
idx = 0;
for kc = 1:nChains
  clear nRotamers indices rmsd rmsdz partition_function secondaryType
  for kr = 1:nResidues(kc)
    nRotamers(kr) = length(sites(kc).residue(kr).rotamers);
    ind = sites(kc).residue(kr).indices;
    indices(kr,:) = ind;
    secondaryType(kr) = model.structures{ind(1)}(ind(2)).residues{ind(3)}.info(ind(4)).secondary;
    rmsd(kr) = NOpos_rmsd(sites(kc).residue(kr).NOpos);
    rmsdz(kr) = NOpos_rmsdz(sites(kc).residue(kr).NOpos);
    partition_function(kr) = sites(kc).residue(kr).partition_function;
    
    idx = idx+1;
    if saveDistanceMatrix
      dist_mat(idx,idx) = rmsd(kr);
    end
    
    if saveStatistics
      [stag,ctag,modelnum,resnum,icode,atag,ltag,resname]=...
        mk_address_parts(ind);
      fprintf(stat_file,'%4i   \t %s    \t %s    \t %2i  \t%6i %s \t %s \t %6i   \t %12.5f   \t %10.3f\n',...
        idx,stag,ctag,modelnum,resnum,icode,resname,nRotamers(kr),partition_function(kr),rmsd(kr));
    end
    
    if saveIndivStatistics
      [stag,ctag,modelnum,resnum,icode,atag,ltag,resname]=mk_address_parts(ind);
      character_table=rotamer_char_table(sites(kc).library);
      all_pops=sites(kc).residue(kr).NOpos(:,4);
      [sort_pops,pop_poi] = sort(all_pops,1,'descend');
      for krot=1:length(sort_pops),
        rotnum=sites(kc).residue(kr).NOpos(pop_poi(krot),5);
        character=character_table{rotnum};
        fprintf(indiv_stat_file,'%4i   \t %s    \t %s    \t %2i  \t%6i %s \t %s \t %6i   \t %12s \t %5.3f\n',...
          idx,stag,ctag,modelnum,resnum,icode,resname,rotnum,character,sort_pops(krot));
      end
    end
    
  end
  
  % Sort with increasing rmsd
  [rmsd,idx1] = sort(rmsd,'ascend');
  nRotamers = nRotamers(idx1);
  indices = indices(idx1,:);
  secondaryType = secondaryType(idx1);
  rmsdz = rmsdz(idx1);
  partition_function = partition_function(idx1);
  tightSite = partition_function<0.05;
  
  % Loop/Helix/Strand labeling statistics
  if (Options.HtmlOutputLevel>=1)    
    fw('<a name="site%i">\n',kc);
    fw('<h2>Site analysis %s</h2>\n<p>',ChainNameHtml{kc});
    secStr = {'loop','helix','strand'};
    for secTyp = 0:2
      n = sum(secondaryType==secTyp);
      if n==0, continue; end
      rmsd_ = rmsd(secondaryType==secTyp);
      fw('%i %s sites, rmsd min/mean/max %4.2f/%4.2f/%4.2f nm<br>\n',...
        n,secStr{secTyp+1},min(rmsd_),sqrt(sum(rmsd_.^2)/n),max(rmsd_));
    end
  end
  
  % List of labeled sites with statistics for each site
  if (Options.HtmlOutputLevel>0)
    fw('<a name="spatial%i">\n',kc);
    fw('<table style:"border:0px;">\n');
    fw('<tr><td colspan="2">Residue<td>label</td><td>location</td><td>center rmsd</td><td>rotamers</td><td>partition function</td></tr>\n');
    for kr = [find(~tightSite) find(tightSite)]    
      ind = indices(kr,:);
      radr = mk_address(ind);
      residueIdx = str2num(radr(strfind(radr,'}')+1:end));
      residueName = model.structures{ind(1)}(ind(2)).residues{ind(3)}.info(ind(4)).name;
      residueName(2:end) = lower(residueName(2:end));
      switch secondaryType(kr)
        case 0, sectype='loop';
        case 1, sectype='helix';
        case 2, sectype='strand';
        otherwise, sectype = '?';
      end
      
      commStr = '';
      if ~isempty(BindingSite) || ~isempty(MetalCenter), % test for binding site or metal center
        rindices=resolve_address(radr);
        [msg,anno]=get_annotations(rindices);
        if ~isempty(BindingSite) && ~isempty(anno)
          if ~isempty(find(anno.keywords==BindingSite, 1)),
            commStr = ' (involved in a binding site)';
          end
        end
        if ~isempty(MetalCenter)
          [msg,aindices]=get_object(radr,'children');
          [maa,naa]=size(aindices);
          for kaa=1:maa
            [msg,anno]=get_annotations(aindices(kaa,:));
            if isempty(anno), continue; end
            if isempty(find(anno.keywords==MetalCenter,1)), continue; end
            aadr = mk_address(aindices(kaa,:));
            paa = strfind(aadr,'.');
            if ~isempty(paa) && paa<length(aadr),
              atom_tag=aadr(paa+1:end);
              commStr = sprintf(' (Atom %s is involved in metal coordination)',atom_tag);
            end
          end
        end
      end
      
      if tightSite(kr), ti=' (tight, labeling might fail)'; else ti=''; end
      
      fw('<tr><td align="right">%3d</td><td>%s</td><td>%s</td><td>%s</td><td>%4.2f nm</td><td>%i</td><td>%7.5f %s</td></tr>\n',...
        residueIdx,residueName,sites(kc).residue(kr).label,sectype,rmsd(kr),nRotamers(kr),partition_function(kr),[ti commStr]);
      if hMain.z_analysis
        %fw(', NO z position r.m.s.d. %4.2f nm\n',rmsdz(poi(kpoi)));
      end
      
    end
    fw('</table>\n');
    fw('<a href="#overview">To overview</A>\n');
  end
    
  % Save binary data
  if mod(Options.matOutput,2)==1,
    outname=sprintf('%s_%s.mat',bin_name,ChainName{kc});
    save(outname,'num','indices','sec','rmsd','rmsdz','partition_function');
  end
  
end

if saveStatistics, fclose(stat_file); end
if saveIndivStatistics, fclose(indiv_stat_file); end

if saveDistanceMatrix
  idx1 = 0;
  for kc1 = 1:nChains
    for kr1 = 1:nResidues(kc1)
      NOpos1 = sites(kc1).residue(kr1).NOpos;
      idx1 = idx1+1;
      idx2 = 0;
      for kc2 = 1:nChains
        for kr2 = 1:nResidues(kc2)
          idx2 = idx2+1;
          if idx2>idx1
            NOpos2 = sites(kc2).residue(kr2).NOpos;
            [rm,sr] = analyze_distribution(NOpos1,NOpos2);
            dist_mat(idx1,idx2) = rm;
            dist_mat(idx2,idx1) = sr;
          end
        end
      end
    end
  end
  for k1 = 1:nTotalResidues
    fprintf(dist_file,'%6.3f\t',dist_mat(k1,:));
    fprintf(dist_file,'\n');
  end
  fclose(dist_file);
end


% Intrachain distances
%--------------------------------------------------------------------------
if (hMain.site_scan_intra>0)
  for kc = 1:nChains
    if nResidues(kc)<2
      if Options.DistanceOutput
        fw('<p>Less than two residues in chain model %s. No distance analysis.\n',ChainName{kc});
      end
      continue
    end
    nPairs = nResidues(kc)*(nResidues(kc)-1)/2;
    if nPairs>Options.maxNumberOfPairsHtml
      if Options.DistanceOutput
        fw('<p>Too many pairs of residues in chain %s. Distance analysis skipped.\n',ChainName{kc});
      end
      continue
    end
    
    r_mean=zeros(1,nPairs);
    r_stddev=zeros(1,nPairs);
    ind1=zeros(nPairs,4);
    ind2=zeros(nPairs,4);
    pfmin=zeros(nPairs,1);
    idx=0;
    for kr1 = 1:nResidues(kc)-1
      pf1 = sites(kc).residue(kr1).partition_function;
      NOpos1 = sites(kc).residue(kr1).NOpos;
      for kr2 = kr1+1:nResidues(kc)
        idx = idx+1;
        ind1(idx,:) = sites(kc).residue(kr1).indices;
        ind2(idx,:) = sites(kc).residue(kr2).indices;
        pf2 = sites(kc).residue(kr2).partition_function;
        pfmin(idx) = min([pf1,pf2]);
        NOpos2 = sites(kc).residue(kr2).NOpos;
        [r_mean(idx),r_stddev(idx)] = analyze_distribution(NOpos1,NOpos2);
      end
    end
    
    if (Options.matOutput>=2)
      matFileName = sprintf('%s_intra_%s.mat',bin_name,ChainName{kc});
      save(matFileName,'pairs','ind1','ind2','pfmin','r_mean','r_stddev');
    end
    
    if (Options.DistanceOutput)
      anchorName = sprintf('intra%i',kc);
      fw('<a name="%s">\n',anchorName);
      fw('<h2>Intrachain distances %s</h2>\n',ChainNameHtml{kc});
      list_distances(anchorName,wfile,ind1,ind2,r_mean,r_stddev,pfmin);
      fprintf(wfile,'<a href="#overview">To overview </a>\n');
    end
    
  end
end

% Analysis of interchain distances
%--------------------------------------------------------------------------
if (hMain.site_scan_inter>0)
  for kc1 = 1:nChains-1
    for kc2 = kc1+1:nChains
      nPairs = nResidues(kc1)*nResidues(kc2);
      if (nPairs==0)
        if (Options.DistanceOutput)
          fw('%s --- %s: No inter-chain site pair. No distance analysis.\n',ChainNameHtml{kc1},ChainNameHtml{kc2});
        end
        continue
      end
      if (nPairs>Options.maxNumberOfPairsHtml)
        if (Options.DistanceOutput)
          fw('%s --- %s: Too many inter-chain site pairs. Distance analysis skipped.<P>\n',ChainName{kc1},ChainName{kc2});
        end
        continue
      end
            
      r_mean=zeros(1,nPairs);
      r_stddev=zeros(1,nPairs);
      ind1=zeros(nPairs,4);
      ind2=zeros(nPairs,4);
      pfmin=zeros(1,nPairs);
      iPair=0;
      for kr1 = 1:nResidues(kc1)
        pf1 = sites(kc1).residue(kr1).partition_function;
        NOpos1=sites(kc1).residue(kr1).NOpos;
        for kr2 = 1:nResidues(kc2)
          if sites(kc1).residue(kr1).indices(4)==sites(kc2).residue(kr2).indices(4) || hMain.site_scan_inter>1,
            iPair = iPair+1;
            ind1(iPair,:) = sites(kc1).residue(kr1).indices;
            ind2(iPair,:) = sites(kc2).residue(kr2).indices;
            pf2 = sites(kc2).residue(kr2).partition_function;
            pfmin(iPair) = min([pf1,pf2]);
            NOpos2 = sites(kc2).residue(kr2).NOpos;
            [r_mean(iPair),r_stddev(iPair)] = analyze_distribution(NOpos1,NOpos2);
          end
        end
      end
      nPairs=iPair;
      ind1=ind1(1:nPairs,:);
      ind2=ind2(1:nPairs,:);
      r_mean=r_mean(1:nPairs);
      r_stddev=r_stddev(1:nPairs);
      pfmin=pfmin(1:nPairs);
      
      if (Options.matOutput>=2)
        matFileName=sprintf('%s_inter_%s_%s.mat',bin_name,ChainName{kc1},ChainName{kc2});
        save(matFileName,'pairs','ind1','ind2','pfmin','r_mean','r_stddev');
      end
      
      if (Options.DistanceOutput)
        fw('<a name="inter%ito%i">\n',kc1,kc2);
        mode='sites';
        if hMain.site_scan_inter==1, mode='equivalent sites'; end
        fw('<h2>Distances <a name="inter%ito%i">between %s in chain models %s and %s</A></h2>\n',kc1,kc2,mode,ChainNameHtml{kc1},ChainNameHtml{kc2});
        %fw('(sorted in categories by ascending relative distribution width)<P>\n');
        list_distances(sprintf('inter%ito%i',kc1,kc2),wfile,ind1,ind2,r_mean,r_stddev,pfmin);
        fprintf(wfile,'<a href="#overview">To overview </a>\n');
      end
            
    end
  end
end

% Analysis of homooligomer distances
%--------------------------------------------------------------------------
if hMain.site_scan_homooligomer && (hMain.site_scan_multiplicity>1)
  multiplicity = hMain.site_scan_multiplicity;
  % make set of turning angles with corresponding populations
  % consider mirror symmetry of distances
  ang=zeros(1,multiplicity);
  apop=zeros(1,multiplicity);
  apop(1)=1;
  cang=0;
  dang=2*pi/multiplicity;
  idx=1;
  for k=2:multiplicity
    cang=cang+dang;
    known=0;
    for kk=1:idx
      if abs(cos(cang)-cos(ang(kk)))< Options.mirror_threshold
        known=1;
        apop(kk)=apop(kk)+1;
        break;
      end
    end
    if ~known, idx=idx+1; ang(idx)=cang; apop(idx)=1; end;
  end;
  ang=ang(1:idx);
  apop=apop(1:idx);
  
  for kc=1:nChains
    if (nResidues(kc)>=1)
      nPairs = nResidues(kc)*(length(ang)-1); % number of spin pairs
      
      r_mean=zeros(1,nPairs);
      r_stddev=zeros(1,nPairs);
      angles=zeros(1,nPairs);
      populations=zeros(1,nPairs);
      ind1=zeros(nPairs,4);
      pfmin=zeros(1,nPairs);
      idx=0;
      for kr1=1:nResidues(kc),
        NOpos1=sites(kc).residue(kr1).NOpos;
        [mm,nn]=size(NOpos1);
        for ka=2:length(ang),
          NOpos2=NOpos1;
          transmat=affine('rotz',ang(ka));
          for kk=1:mm,
            NOpos2(kk,1:3)=affine_trafo_point(NOpos1(kk,1:3),transmat);
          end;
          idx=idx+1;
          angles(idx)=180*ang(ka)/pi;
          populations(idx)=apop(ka);
          ind1(idx,:) = sites(kc).residue(kr1).indices;
          pfmin(idx) = sites(kc).residue(kr1).partition_function;
          [r_mean(idx),r_stddev(idx)]=analyze_distribution(NOpos1,NOpos2);
        end
      end
      
      if (Options.matOutput>=2)
        matFileName=sprintf('%s_oligo_%s.mat',bin_name,ChainName{kc});
        save(matFileName,'pairs','ind1','ind2','pfmin','r_mean','r_stddev');
      end
            
      if (Options.DistanceOutput)
        anchorName = sprintf('oligo%i',kc);
        fw('<a name="%s">\n',anchorName);
        fw('<h2>Homooligomer distances for chain model %s with multiplicity %i</h2>\n',ChainName{kc},hMain.site_scan_multiplicity);
        fw('(sorted in categories by ascending relative distribution width)<P>\n');
        list_distances(anchorName,wfile,ind1,ind1,r_mean,r_stddev,pfmin,angles,populations);
        fprintf(wfile,'<a href="#overview">To overview</a>\n');
      end
      
    end
  end
end
fclose(wfile);

%==========================================================================

function rmsd = NOpos_rmsd(NOall)
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

function list_distances(distanceStr,wfile,ind1,ind2,r_mean,r_stddev,pfmin,angles,apop)

global model
anchor{6} = '';
anchor{1} = '';
anchor{2} = '';
anchor{3} = '';
anchor{4} = '';
anchor{5} = '';
categoryDescr{6}='Distance probably too long for useful measurement (over 10.0 nm)';
categoryDescr{1}='At least one site probably too tight for labeling';
categoryDescr{2}='Distance too short for proper measurement (less than 1.5 nm)';
categoryDescr{3}='intermediate (should not occur)'; % 'Most favorable for CW EPR measurements (0.8 nm ... 1.8 nm)';
categoryDescr{4}='Most favorable for FRET measurements (1.5 nm ... 8.0 nm)';
categoryDescr{5}='FRET measurement possible, but rather imprecise (8.0 nm ... 10.0 nm)';
categoryStr = {'tight','short','intmediate','FRET','FRETlong','long'};
category(r_mean<1.5) = 2;
% category(r_mean>=1.5) = 3;
category(r_mean>=1.5) = 4;
category(r_mean>=8.0) = 5;
category(r_mean>=10.0) = 6;
category(pfmin<0.05) = 1;
fprintf(wfile,'<table>\n');
fprintf(wfile,'<tr><td>Site 1</td><td>Site 2</td><td>distance (nm)</td><td>rmsd (nm)</td><td>rel.width</td></tr>\n');
for iCat = [4 3 5 2 6 1]
  idx = find(category==iCat);
  if isempty(idx), continue; end
  
  [dummy,iidx] = sort(r_mean(idx));
  idx = idx(iidx);
  
  linkStr = sprintf('<a name="%s_%s"></a>\n',distanceStr,categoryStr{iCat});
  fprintf(wfile,'<tr><td colspan=5>%s<span class="distcategory">%s</span> <a href="#overview">to overview</a></td></tr>\n',linkStr,categoryDescr{iCat});
  rel_width=r_stddev./r_mean;
  for k=1:length(idx)
    ind1a=ind1(idx(k),:);
    ind2a=ind2(idx(k),:);
    adr1=mk_address(ind1a);
    adr2=mk_address(ind2a);
    resiName1=model.structures{ind1a(1)}(ind1a(2)).residues{ind1a(3)}.info(ind1a(4)).name;
    resiName1(2:end)=lower(resiName1(2:end));
    resiName2=model.structures{ind2a(1)}(ind2a(2)).residues{ind2a(3)}.info(ind2a(4)).name;
    resiName2(2:end)=lower(resiName2(2:end));
    chainName1 = makeHtmlChainName(adr1);
    chainName2 = makeHtmlChainName(adr2);
    %if strcmp(chainName1,chainName2), chainName2=''; end
    resiIdx1 = adr1(strfind(adr1,'}')+1:end);
    resiIdx2 = adr2(strfind(adr2,'}')+1:end);
    if (nargin<8)
      fprintf(wfile,'<tr><td>%s (%s)</td><td>%s (%s)</td><td>%4.2f</td><td>%4.2f</td><td>% 3.1f%%</td></tr>\n',...
        chainName1,[resiIdx1 resiName1],...
        chainName2,[resiIdx2 resiName2],...
        r_mean(idx(k)),r_stddev(idx(k)),100*rel_width(idx(k)));
      %fprintf(wfile,'%s (%s) --- %s (%s):   %4.2f nm, std.dev. %4.2f nm, rel. width %3.1f%%<br>',...
      %  adr1,name1,adr2,name2,r_mean(idx(k)),r_stddev(idx(k)),100*rel_width(idx(k)));
    else
      fprintf(wfile,'Distance between residues %s (%s) at rotation angle %4.1f° occurs %i times and is %4.2f nm, , std.dev. %4.2f nm, rel. width %3.1f\n',...
        adr1,resiName1,angles(idx(k)),apop(idx(k)),r_mean(idx(k)),r_stddev(idx(k)),100*rel_width(idx(k)));
    end
  end
end
fprintf(wfile,'</table>\n');

function [rm,sr] = analyze_distribution(NOpos1,NOpos2)

pop1=NOpos1(:,4);
n1=length(pop1);
xyz1=NOpos1(:,1:3)/10; % divided by 10 for Å -> nm

pop2=NOpos2(:,4);
n2=length(pop2);
xyz2=NOpos2(:,1:3)/10; % divided by 10 for Å -> nm

xyz1sq = repmat(sum(xyz1.^2,2),1,n2);
xyz2sq = repmat(sum(xyz2.^2,2),1,n1).';
r = sqrt(abs(xyz1sq + xyz2sq - 2*xyz1*xyz2.'));
r = r(:);
pop = pop1*pop2.';
pop = pop(:);

rm = sum(r.*pop)/sum(pop);
dr = r-rm;
sr = sqrt(0.01+n1*n2*sum(pop.*dr.^2)/(sum(pop)*(n1*n2-1)));

function s = makeHtmlChainName(chainstr)
idx1 = strfind(chainstr,'[');
idx2 = strfind(chainstr,']');
pdb = chainstr(idx1+1:idx2-1);
idx1 = strfind(chainstr,'(');
idx2 = strfind(chainstr,')');
chain = chainstr(idx1+1:idx2-1);
idx1 = strfind(chainstr,'{');
idx2 = strfind(chainstr,'}');
model = chainstr(idx1+1:idx2-1);
s = ['<span class="pdb">' pdb '</span>/'...
     '<span class="chain">' chain '</span>'...
     '<span class="model">' model '</span>'];
return
