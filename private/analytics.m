function analytics

global hMain

state1='2X79';
state2='2JLO';

template='2WSX';
target='';

dali_lite=true;

basis=50;

gap_mode=true; % true: gaps are neglected in fit mode, false: gap repair by Modeller is accepted
give_title=false;

% 'changes', 'covariance', 'modes', 'spectrum', 'core_covariance', 'fit',
% 'theta'


an_mode='constraints';

TMD_colors=[0,0,139;65,105,225;0,139,139;0,100,0;0,128,0;127,255,0;...
    255,255,0;255,165,0;255,0,0;139,0,0];
TMD_colors=TMD_colors/255;

switch state1
    case {'2JLN','2JLO','2X79'}
        protein='Mhp1';
    case {'2A65','3QS4','3F3A','3GJC','3TT1','3TT3'}
        protein='LeuT';
    case {'2WSX','3HFX'}
        protein='CaiT';
    case {'3LRB','3L1L','3NCY','3OB6'}
        protein='AdiC';
    case {'3DH4','2XQ2'}
        protein='vSGLT';
    case {'2WIT','3P03'}
        protein='BetP';
    case {'3GIA','3GI8''3GI9'}
        protein='ApcT';
end;

switch template
    case {'2JLN','2JLO','2X79'}
        protein2='Mhp1';
    case {'2A65','3QS4','3F3A','3GJC','3TT1','3TT3'}
        protein2='LeuT';
    case {'2WSX','3HFX'}
        protein2='CaiT';
    case {'3LRB','3L1L','3NCY','3OB6'}
        protein2='AdiC';
    case {'3DH4','2XQ2'}
        protein2='vSGLT';
    case {'2WIT','3P03'}
        protein2='BetP';
    case {'3GIA','3GI8''3GI9'}
        protein2='ApcT';
end;

mask1=[];
switch state1
    case '2A65'
        mask1=133:134;
    case '3QS4'
        mask1=132:134;
    case '3F3A'
        mask1=132:134;
    case '3LRB'
        mask1=253:272;
    case '3L1L'
        mask1=181:191;
    case '3NCY'
        mask1=[177:181 317:321];
    case '3OB6'
        mask1=345:350;
    case '3DH4'
        mask1=179:184;
    case '2XQ2'
        mask1=[182:184 323:327];
    case '3TT1'
        mask1=130:134;
end;

mask2=[];
switch state2
    case '2A65'
        mask2=133:134;
    case '3QS4'
        mask2=132:134;
    case '3F3A'
        mask2=132:134;
    case '3LRB'
        mask2=253:272;
    case '3L1L'
        mask2=181:191;
    case '3NCY'
        mask2=[177:181 317:321];
    case '3OB6'
        mask2=345:350;
    case '3DH4'
        mask1=179:184;
    case '2XQ2'
        mask1=[182:184 323:327];
    case '3TT1'
        mask1=130:134;
end;

maskt=[];
switch template
    case '2A65'
        maskt=133:134;
    case '3QS4'
        maskt=132:134;
    case '3F3A'
        maskt=132:134;
    case '3LRB'
        maskt=253:272;
    case '3L1L'
        maskt=181:191;
    case '3NCY'
        maskt=[177:181 317:321];
    case '3OB6'
        maskt=345:350;
    case '3DH4'
        maskt=179:184;
    case '2XQ2'
        maskt=[182:184 323:327];
    case '3TT1'
        maskt=130:134;
end;

selection=[];

if strcmp(state1,'2JLN') && strcmp(state2,'2JLO'),
    selection=29:379;
end;

if strcmp(state1,'2JLN') && strcmp(state2,'2X79'),
    selection=29:379;
end;

if strcmp(state1,'2A65') && strcmp(state2,'3F3A'),
    selection=[11:131 135:424];
end;

if strcmp(state1,'3F3A') && strcmp(state2,'3GJC'),
    selection=[11:131 135:424];
end;

if strcmp(state1,'2A65') && strcmp(state2,'3GJC'),
    selection=[11:132 135:424];
end;

if strcmp(state1,'2A65') && strcmp(state2,'3TT1'),
    selection=[11:129 135:424];
end;

if strcmp(state1,'3TT1') && strcmp(state2,'3TT3'),
    selection=[11:129 135:424];
end;

if strcmp(state1,'2A65') && strcmp(state2,'3TT3'),
    selection=[11:132 135:424];
end;

if strcmp(state1,'3L1L') && strcmp(state2,'3LRB'),
    selection=[11:180 192:252 273:374];
end;

if strcmp(state1,'3L1L') && strcmp(state2,'3NCY'),
    selection=[11:176 192:316 322:374];
end;

if strcmp(state1,'3L1L') && strcmp(state2,'3OB6'),
    selection=[11:180 192:344 351:374];
end;

if strcmp(state1,'3LRB') && strcmp(state2,'3OB6'),
    selection=[11:252 273:344 351:374];
end;

if strcmp(state1,'3DH4') && strcmp(state2,'2XQ2'),
    selection=[53:178 185:322 328:447];
end;

if strcmp(state1,'2WIT') && strcmp(state2,'3P03'),
    selection=[138:271 276:538];
end;

if dali_lite, 
    protein=[protein '_DL']; 
    protein2=[protein2 '_DL']; 
    template=[template '_DL'];
    state1=[state1 '_DL']; 
    state2=[state2 '_DL']; 
end;

switch protein
    case 'Mhp1'
        core_def=[29,53;59,84;102,133;142,157;164,187;206,230;249,274;294,325;340,355;356,379];
    case 'Mhp1_DL'
        core_def=[11,14;21,41;42,54;56,66;69,74;76,87;91,97;99,102;103,134;140,143;144,160;161,188;192,198;206,227;228,236;248,279;282,291;298,330;334,351;353,356;362,369;370,382;385,388];
    case 'CaiT'
        core_def=[88,112;133,158;188,219;232,247;255,278;309,333;347,372;403,434;449,464;469,492];
    case 'CaiT_DL'
        core_def=[9,12;80,100;104,116;133,143;144,149;150,161;169,175;188,191;193,224;225,228;230,246;252,279;280,286;310,331;333,341;343,374;393,402;403,435;447,464;467,470;474,481;483,495;496,499];
    case 'BetP'
        core_def=[138,162;180,205;234,265;280,295;302,325;359,383;397,422;449,480;492,507;515,538];
    case 'vSGLT'
        core_def=[53,77;82,107;126,157;162,177;185,208;249,273;282,307;346,377;399,414;424,447];
    case 'LeuT'
        core_def=[11,35;42,67;89,120;168,183;191,214;240,264;275,300;336,367;379,394;401,424];
    case 'LeuT_DL'
        core_def=[19,42;44,52;54,71;72,77;81,85;86,125;150,153;166,185;186,192;194,215;239,253;256,259;260,273;274,307;310,316;317,331;336,367;369,372;374,398;400,407;409,424;426,431];
    case 'AdiC'
        core_def=[11,35;41,66;81,112;125,140;144,167;190,214;224,249;274,305;327,342;351,374];
    case 'ApcT'
        core_def=[10,34;40,65;85,116;125,140;147,170;186,210;221,246;268,299;321,336;340,363];
end;

switch protein2
    case 'Mhp1'
        core_def2=[29,53;59,84;102,133;142,157;164,187;206,230;249,274;294,325;340,355;356,379];
    case 'Mhp1_DL'
        core_def2=[11,14;21,41;42,54;56,66;69,74;76,87;91,97;99,102;103,134;140,143;144,160;161,188;192,198;206,227;228,236;248,279;282,291;298,330;334,351;353,356;362,369;370,382;385,388];
    case 'CaiT'
        core_def2=[88,112;133,158;188,219;232,247;255,278;309,333;347,372;403,434;449,464;469,492];
    case 'CaiT_DL'
        core_def2=[9,12;80,100;104,116;133,143;144,149;150,161;169,175;188,191;193,224;225,228;230,246;252,279;280,286;310,331;333,341;343,374;393,402;403,435;447,464;467,470;474,481;483,495;496,499];
    case 'BetP'
        core_def2=[138,162;180,205;234,265;280,295;302,325;359,383;397,422;449,480;492,507;515,538];
    case 'vSGLT'
        core_def2=[53,77;82,107;126,157;162,177;185,208;249,273;282,307;346,377;399,414;424,447];
    case 'LeuT'
        core_def2=[11,35;42,67;89,120;168,183;191,214;240,264;275,300;336,367;379,394;401,424];
    case 'LeuT_DL'
        core_def2=[19,42;44,52;54,71;72,77;81,85;86,125;150,153;166,185;186,192;194,215;239,253;256,259;260,273;274,307;310,316;317,331;336,367;369,372;374,398;400,407;409,424;426,431];
    case 'AdiC'
        core_def2=[11,35;41,66;81,112;125,140;144,167;190,214;224,249;274,305;327,342;351,374];
    case 'ApcT'
        core_def2=[10,34;40,65;85,116;125,140;147,170;186,210;221,246;268,299;321,336;340,363];
end;

if ~isempty(selection),
    selection=selection-core_def(1,1)+1;
end;

switch an_mode
    case 'changes'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coarse0=coarse;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        outfile=sprintf('%s_to_%s_changes.txt',state1,state2);
        fid=fopen(outfile,'wt');
        fprintf(fid,'Coordinate changes of TMDs and loops for %s between states %s and %s\n\n',protein,state1,state2);
        rel_core=core_def-core_def(1,1)+1; % relative core definition
        rangeb=[rel_core(1,1):rel_core(1,2),rel_core(2,1):rel_core(2,2),rel_core(6,1):rel_core(6,2),rel_core(7,1):rel_core(7,2)];
        rangeh=[rel_core(3,1):rel_core(3,2),rel_core(4,1):rel_core(4,2),rel_core(8,1):rel_core(8,2),rel_core(9,1):rel_core(9,2)];
        rel_core=rel_core-1;
        rmsdb0=rmsd_superimpose(coarse0.Ca_coor(rangeb,:),coarse.Ca_coor(rangeb,:));
        rmsdh=rmsd_superimpose(coarse0.Ca_coor(rangeh,:),coarse.Ca_coor(rangeh,:));
        diffb=coarse0.Ca_coor(rangeb,:)-coarse.Ca_coor(rangeb,:);
        diffh=coarse0.Ca_coor(rangeh,:)-coarse.Ca_coor(rangeh,:);
        rmsdba=sqrt(sum(sum(diffb.^2))/length(diffb));
        rmsdha=sqrt(sum(sum(diffh.^2))/length(diffh));
        diff=coarse0.Ca_coor-coarse.Ca_coor;
        per_residue=sqrt(sum(diff.^2,2));
        mask1=mask1-core_def(1,1)+1;
        mask2=mask2-core_def(1,1)+1;
        per_residue(mask1)=0;
        per_residue(mask2)=0;
        [m,n]=size(rel_core);
        figure(1); clf;
        hold on;
        for k=1:m,
            res_ax=core_def(k,1):core_def(k,2);
            pr=per_residue(rel_core(k,1)+1:rel_core(k,2)+1);
            plot(res_ax,pr,'Color',TMD_colors(k,:));
            if k<m,
                res_ax=core_def(k,2)+1:core_def(k+1,1)-1;
                pr=per_residue(rel_core(k,2)+2:rel_core(k+1,1));
                plot(res_ax,pr,'k');
            end;
        end;
        set(gca,'FontSize',14);
        axis([core_def(1,1),core_def(m,2),0,16]);
        title('Coordinate change per residue');
        xlabel('Residue number');
        rmsda=sqrt(sum(sum(diff.^2))/length(diff));
        rmsdb=rmsd_superimpose(coarse0.Ca_coor,coarse.Ca_coor);
        fprintf(fid,'General: global %4.2f Å, internal %4.2f Å\n',rmsda,rmsdb);
        fprintf(fid,'Bundle : global %4.2f Å, internal %4.2f Å\n',rmsdba,rmsdb0);
        fprintf(fid,'Hash : global %4.2f Å, internal %4.2f Å\n',rmsdha,rmsdh);
        ivec=[];
        for k=1:m,
            ivec=[ivec rel_core(k,1)+1:rel_core(k,2)+1];
        end;
        coor0=coarse0.Ca_coor(ivec,:);
        coor=coarse.Ca_coor(ivec,:);
        diff=coor0-coor;
        rmsda=sqrt(sum(sum(diff.^2))/length(diff));
        rmsdb=rmsd_superimpose(coor0,coor);
        fprintf(fid,'TMDs  : global %4.2f Å, internal %4.2f Å\n',rmsda,rmsdb);
        TMD=true;
        k=1;
        while k<=m,
            if TMD,
                ivec=rel_core(k,1)+1:rel_core(k,2)+1;
            elseif k<m
                ivec=rel_core(k,2)+2:rel_core(k+1,1);
            end;
            coor0=coarse0.Ca_coor(ivec,:);
            coor=coarse.Ca_coor(ivec,:);
            if ~isempty(coor0) && ~isempty(coor),
                diff=coor0-coor;
                rmsda=sqrt(sum(sum(diff.^2))/length(diff));
                rmsdb=rmsd_superimpose(coor0,coor);
            else
                rmsda=-1;
                rmsdb=-1;
            end;
            if TMD,
                fprintf(fid,'TMD %2i: global %4.2f Å, internal %4.2f Å\n',k,rmsda,rmsdb);
            elseif k<m
                fprintf(fid,'L%2i-%2i: global %4.2f Å, internal %4.2f Å\n',k,k+1,rmsda,rmsdb);                
            end;
            TMD=~TMD;
            if TMD,
                k=k+1;
            end;
        end;
        fclose(fid);
        hMain.report_file=outfile;
        report_editor;
    case 'covariance'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        tcorr0=ENM.covariance;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        tcorr=ENM.covariance;
        sc=sum(sum(tcorr0.*tcorr))/sum(sum(tcorr0.*tcorr0));
        if isnan(sc), sc=1; end;
        diff_corr=sc*tcorr-tcorr0;
        chisq=sum(sum((tcorr0-sc*tcorr).^2));
        err_left=chisq/sum(sum(tcorr0.^2));
        fprintf(1,'Difference of covariance matrices: %6.4f\n',err_left);
        figure;
        tcorr0=255*tcorr0/max(max(tcorr0));
        image(tcorr0);
        hold on;
        rel_core=core_def-core_def(1,1)+1;
        aa=1;
        ee=rel_core(end,2);
        for k=1:10,
            h=plot([rel_core(k,1),rel_core(k,1)],[aa,ee],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([rel_core(k,2),rel_core(k,2)],[aa,ee],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([aa,ee],[rel_core(k,1),rel_core(k,1)],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([aa,ee],[rel_core(k,2),rel_core(k,2)],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
        end;
        set(gca,'YDir','normal');
        axis equal
        axis tight
        axis off
        if give_title,
            title(sprintf('Covariance matrix for state %s',state1));
        end;
        greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
        greyscale=flipud(greyscale');
        colormap(greyscale.^16); 
        figure;
%         v=max(max(tcorr))*[0.0025,0.005,0.01,0.02,0.04,0.08];
%         [C,h]=contourf(tcorr,v);
%         set(h,'LineColor',[1,1,1]);
%         axis equal
%         axis tight
        tcorr=255*tcorr/max(max(tcorr));
        image(tcorr);
        hold on;
        rel_core=core_def-core_def(1,1)+1;
        aa=1;
        ee=rel_core(end,2);
        for k=1:10,
            h=plot([rel_core(k,1),rel_core(k,1)],[aa,ee],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([rel_core(k,2),rel_core(k,2)],[aa,ee],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([aa,ee],[rel_core(k,1),rel_core(k,1)],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([aa,ee],[rel_core(k,2),rel_core(k,2)],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
        end;
        set(gca,'YDir','normal');
        axis equal
        axis tight
        axis off
        if give_title,
            title(sprintf('Covariance matrix for state %s',state2));
        end;
        greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
        greyscale=flipud(greyscale');
        colormap(greyscale.^16); 
        figure;
        tcorrd=tcorr0;
        [m,n]=size(tcorrd);
        for k=1:n,
            for kk=k+1:m,
                tcorrd(k,kk)=tcorr(k,kk);
            end;
        end;
        image(2*tcorrd);
        hold on;
        rel_core=core_def-core_def(1,1)+1;
        aa=1;
        ee=rel_core(end,2);
        for k=1:10,
            h=plot([rel_core(k,1),rel_core(k,1)],[aa,ee],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([rel_core(k,2),rel_core(k,2)],[aa,ee],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([aa,ee],[rel_core(k,1),rel_core(k,1)],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
            h=plot([aa,ee],[rel_core(k,2),rel_core(k,2)],'LineWidth',2.5);
            set(h,'Color',TMD_colors(k,:));
        end;
        set(gca,'YDir','normal');
        axis equal
        axis tight
        axis off
        if give_title,
            title(sprintf('Combined covariance matrix for states %s,%s',state1,state2));
        end;
        greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
        greyscale=flipud(greyscale');
        colormap(greyscale.^16); 
    case 'core_covariance'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        tcorr0=ENM.core_covariance;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        tcorr=ENM.core_covariance;
        sc=sum(sum(tcorr0.*tcorr))/sum(sum(tcorr0.*tcorr0));
        if isnan(sc), sc=1; end;
        chisq=sum(sum((tcorr0-sc*tcorr).^2));
        err_left=chisq/sum(sum(tcorr0.^2));
        fprintf(1,'Difference of core covariance matrices: %6.4f\n',err_left);
        figure;
        tcorr0=255*tcorr0/max(max(tcorr0));
        image(tcorr0);
        set(gca,'YDir','normal');
        axis equal
        axis tight
        title(sprintf('Core covariance matrix for state %s',state1));
        greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
        greyscale=flipud(greyscale');
        colormap(greyscale.^16); 
        figure;
        tcorr=255*tcorr/max(max(tcorr));
        image(tcorr);
        set(gca,'YDir','normal');
        axis equal
        axis tight
        title(sprintf('Core covariance matrix for state %s',state2));
        greyscale=[linspace(0,1,100);linspace(0,1,100);linspace(0,1,100)];
        greyscale=flipud(greyscale');
        colormap(greyscale.^16); 
    case 'modes'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        u0=ENM.u;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        u1=ENM.u;
        [m,m1]=size(u0);
        overlap=zeros(m,m);
        for k0=1:m,
            evec0=u0(:,k0);
            evec0=evec0/norm(evec0);
            for k1=1:m,
                evec1=u1(:,k1);
                evec1=evec1/norm(evec1);
                overlap(k0,k1)=sum(evec0.*evec1);
            end;
        end;
        assignment=zeros(1,m);
        quality=zeros(1,m);
        modax=1:m;
        for k=1:m,
            ol=abs(overlap(k,:));
            [ols,poi]=sort(ol,2,'descend');
            found=false;
            kk=1;
            while ~found && kk<=m,
                used=find(assignment==poi(kk));
                if isempty(used),
                    assignment(k)=poi(kk);
                    quality(k)=ols(kk);
                    found=true;
                else
                    kk=kk+1;
                end;
            end;
        end;
        figure(1); clf;
        plot(modax,assignment,'k.');
        title(sprintf('Mode correspondence between states %s and %s',state1,state2));
        figure(2); clf;
        plot(modax,quality,'k.');
        title(sprintf('Overlap of assigned modes between states %s and %s',state1,state2));
    case 'spectrum'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coor0=coarse.Ca_coor;
        u0=ENM.u;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        coor=coarse.Ca_coor;
        u1=ENM.u;
        [m,m1]=size(u0);
        diff=coor-coor0;
        [n,n1]=size(diff);
        change=reshape(diff',1,3*n);
        coeff0=u0\change';
        coeff0b=u0(:,1:basis)\change';
        coeff1=u1\(-change');
        coeff1b=u1(:,1:basis)\(-change');
        diff_test=zeros(n,3);
        diff_test_b=zeros(n,3);
        for k=1:m,
            evec=u0(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test=diff_test+mode1'*coeff0(k);
        end;
        for k=1:basis,
            evec=u0(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test_b=diff_test_b+mode1'*coeff0b(k);
        end;
        coor_b=coor0+diff_test_b;
        rmsd=rmsd_superimpose(coor0,coor);
        rmsdb=rmsd_superimpose(coor,coor_b);
        fprintf('%s to %s:\n',state1,state2);
        fprintf(1,'%i modes account for %4.2f%% of change. Deviation is %4.2f Å\n',basis,100*(rmsd-rmsdb)/rmsd,rmsdb);
        figure(1); clf;
        plot(diff(:,1),'r');
        hold on;
        plot(diff(:,2),'g');
        plot(diff(:,3),'b');
        plot(diff_test_b(:,1),'k:');
        plot(diff_test_b(:,2),'y:');
        plot(diff_test_b(:,3),'c:');
        title(sprintf('Test of forward coordinates %s to %s',state1,state2));
        diff_test2=zeros(n,3);
        diff_test2_b=zeros(n,3);
        for k=1:m,
            evec=u1(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test2=diff_test2+mode1'*coeff1(k);
        end;
        for k=1:basis,
            evec=u1(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test2_b=diff_test2_b+mode1'*coeff1b(k);
        end;
        coor2_b=coor+diff_test2_b;
        rmsd=rmsd_superimpose(coor0,coor);
        rmsdb=rmsd_superimpose(coor0,coor2_b);
        fprintf('%s to %s:\n',state2,state1);
        fprintf(1,'%i modes account for %4.2f%% of change. Deviation is %4.2f Å\n',basis,100*(rmsd-rmsdb)/rmsd,rmsdb);
        figure(2); clf;
        plot(-diff(:,1),'r');
        hold on;
        plot(-diff(:,2),'g');
        plot(-diff(:,3),'b');
        plot(diff_test2_b(:,1),'k:');
        plot(diff_test2_b(:,2),'y:');
        plot(diff_test2_b(:,3),'c:');
        title(sprintf('Test of backward coordinates %s to %s',state1,state2));
        figure(3); clf;
        plot(abs(coeff0),'k.');
        title(sprintf('Mode spectrum for structural change between states %s and %s',state1,state2));
        hold on;
        set(gca,'FontSize',12);
        axis([0,length(coeff0),0,1.1*max(abs(coeff0))]);
        figure(4); clf;
        plot(abs(coeff1),'k.');
        hold on;
        title(sprintf('Mode spectrum for structural change between states %s and %s',state2,state1));
        coeff0=abs(coeff0(7:end))';
        ax0=1:length(coeff0);
        mom1_0=sum(coeff0.*ax0)/sum(coeff0);
        coeff1=abs(coeff1(7:end))';
        ax1=1:length(coeff1);
        mom1_1=sum(coeff1.*ax1)/sum(coeff1);
        fprintf(1,'First moment of linear coefficients for state %s is %5.1f\n',state1,mom1_0);
        fprintf(1,'First moment of linear coefficients for state %s is %5.1f\n',state2,mom1_1);
    case 'check'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coor0=coarse.Ca_coor;
        u0=ENM.u;
        lambda0=ENM.lambda;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        coor=coarse.Ca_coor;
        u1=ENM.u;
        [m,m1]=size(u0);
        diff=coor-coor0;
        [n,n1]=size(diff);
        fprintf(1,'%i Calpha atoms in structure %s\n',n,state1);
        change=reshape(diff',1,3*n);
        change=change';
        for zz=1:100,
            coeff0=u0\change;
        end;
        for zz=1:100,
            coeff0a=zeros(m,1);
            for k=1:m,
                coeff0a(k)=dot(u0(:,k),change);
            end;
        end;
        for zz=1:100,
            coeff0a=zeros(m,1);
            for k=1:m,
                coeff0a(k)=sum(u0(:,k).*change);
            end;
        end;
        coeff0b=u0(:,1:basis)\change;
        coeff0c=coeff0(1:basis);
        diff0=sqrt(sum((coeff0-coeff0a).^2));
        fprintf(1,'Coefficient difference dot product: %5.3f\n',diff0);
        diff=sqrt(sum((coeff0(1:basis)-coeff0b).^2));
        fprintf(1,'Coefficient difference: %5.3f\n',diff);
        for zz=1:100,
            diff_test=zeros(n,3);
            for k=1:m,
                evec=u0(:,k);
                mode1=reshape(evec,3,m/3);
                diff_test=diff_test+mode1'*coeff0(k);
            end;
            coor_x=coor0+diff_test;
            rmsd=rmsd_superimpose(coor,coor_x);
        end;
        fprintf(1,'Full solution r.m.s.d.: %5.3f\n',rmsd);
        for zz=1:100,
            diff_test=zeros(3*n,1);
            for k=1:m,
                diff_test=diff_test+u0(:,k)*coeff0(k);
            end;
            diff_test=reshape(diff_test,3,m/3);
            coor_y=coor0+diff_test';
            rmsd=rmsd_superimpose(coor,coor_y);
        end;
        fprintf(1,'Full solution alternative r.m.s.d.: %5.3f\n',rmsd);
        for zz=1:100,
            diff_test=u0*coeff0;
            diff_test=reshape(diff_test,3,m/3);
            coor_z=coor0+diff_test';
            rmsd=rmsd_superimpose(coor,coor_z);
        end;
        fprintf(1,'Full solution alternative 2 r.m.s.d.: %5.3f\n',rmsd);
        diff_test_b=zeros(n,3);
        for k=1:basis,
            evec=u0(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test_b=diff_test_b+mode1'*coeff0b(k);
        end;
        coor_b=coor0+diff_test_b;
        rmsdb=rmsd_superimpose(coor,coor_b);
        fprintf(1,'Partial solution r.m.s.d.: %5.3f\n',rmsdb);
        diff_test=zeros(3*n,1);
        for k=1:basis,
            diff_test=diff_test+u0(:,k)*coeff0c(k);
        end;
        diff_test=reshape(diff_test,3,m/3);
        coor_y=coor0+diff_test';
        rmsdb=rmsd_superimpose(coor,coor_y);
        fprintf(1,'Partial solution alternative r.m.s.d.: %5.3f\n',rmsdb);
        diff_test=u0(:,1:basis)*coeff0c;
        diff_test=reshape(diff_test,3,m/3);
        coor_z=coor0+diff_test';
        rmsdb=rmsd_superimpose(coor,coor_z);
        fprintf(1,'Partial solution alternative 2 r.m.s.d.: %5.3f\n',rmsdb);
        for zz=1:100,
            kappa = mode_collectivity(u0);
        end;
        figure(1); clf;
        plot(kappa,'.');
        for zz=1:100,
            [kappas,poi]=sort(kappa,'descend');
        end;
        hold on;
        plot(kappas,'.r');
        figure(2); clf;
        plot(poi,'k.');
        ub=u0(:,poi(1:basis));
        coeff=coeff0a(poi(1:basis));
        diff_test=ub*coeff;
        diff_test=reshape(diff_test,3,m/3);
        coor_a=coor0+diff_test';
        rmsdc=rmsd_superimpose(coor,coor_a);
        fprintf(1,'Partial solution collectivity r.m.s.d.: %5.3f\n',rmsdc);
        u2=zeros(m,basis);
        for k=1:basis,
            u2(:,k)=u0(:,poi(k));
        end;
        coeff0c=u2\change;
        diff_test=u2*coeff0c;
        diff_test=reshape(diff_test,3,m/3);
        coor_a=coor0+diff_test';
        rmsdc=rmsd_superimpose(coor,coor_a);
        fprintf(1,'Partial solution collectivity alternative r.m.s.d.: %5.3f\n',rmsdc);  
        for zz=1:100,
            [coeff,poi]=sort(abs(coeff0a),'descend');
        end;
        figure(3); clf;
        plot(coeff,'k.');
        coeff=coeff0a(poi(1:basis));
        uc=u0(:,poi(1:basis));
        diff_test=uc*coeff;
        diff_test=reshape(diff_test,3,m/3);
        coor_b=coor0+diff_test';
        rmsdd=rmsd_superimpose(coor,coor_b);
        fprintf(1,'Partial solution overlap r.m.s.d.: %5.3f\n',rmsdd);
        figure(4); clf;
        plot(lambda0,'k.');
    case 'fit'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coor0=coarse.Ca_coor;
        u0=ENM.u;
        [m,m1]=size(u0);
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        target=coarse.Ca_coor;
        diff=target-coor0;
        [n,n1]=size(diff);
        change=reshape(diff',1,3*n);
        coeff0b=u0(:,1:basis)\change';
        diff_test_b=zeros(n,3);
        for k=1:basis,
            evec=u0(:,k);
            mode1=reshape(evec,3,m/3);
            diff_test_b=diff_test_b+mode1'*coeff0b(k);
        end;
        coor_b=coor0+diff_test_b;
        rmsd=rmsd_superimpose(target,coor0);
        rmsdb=rmsd_superimpose(target,coor_b);
        fprintf('%s to %s:\n',state1,state2);
        if gap_mode && ~isempty(selection),
            rmsd=rmsd_superimpose(target(selection,:),coor0(selection,:));
            rmsdb=rmsd_superimpose(target(selection,:),coor_b(selection,:));
            fprintf(1,'%i modes account for %4.2f%% of change. Deviation is %4.2f Å\n',basis,100*(rmsd-rmsdb)/rmsd,rmsdb);
            [rmsd,coor1a] = fit_by_ANM(coor0(selection,:),target(selection,:),basis);
            ermsd=rmsd_superimpose(target(selection,:),coor1a);
            coor1=coor0;
            coor1(selection,:)=coor1a;
        else
            fprintf(1,'%i modes account for %4.2f%% of change. Deviation is %4.2f Å\n',basis,100*(rmsd-rmsdb)/rmsd,rmsdb);
            [rmsd,coor1] = fit_by_ANM(coor0,target,basis);
            ermsd=rmsd_superimpose(target,coor1);
        end;
        
        fprintf(1,'Final r.m.s.d. after %i iterations: %4.2f Å\n',length(rmsd),ermsd);
        figure(1);
        hold on
        plot([0,length(rmsd)],[rmsdb,rmsdb],'r:');
        diff=coor1-target;
        per_residue=sqrt(sum(diff.^2,2));
        mask1=mask1-core_def(1,1)+1;
        mask2=mask2-core_def(1,1)+1;
        per_residue(mask1)=0;
        per_residue(mask2)=0;
        rel_core=core_def-core_def(1,1); % relative core definition
        [m,n]=size(rel_core);
        figure(2); clf;
        hold on;
        for k=1:m,
            res_ax=core_def(k,1):core_def(k,2);
            pr=per_residue(rel_core(k,1)+1:rel_core(k,2)+1);
            plot(res_ax,pr,'Color',TMD_colors(k,:));
            if k<m,
                res_ax=core_def(k,2)+1:core_def(k+1,1)-1;
                pr=per_residue(rel_core(k,2)+2:rel_core(k+1,1));
                plot(res_ax,pr,'k');
            end;
        end;
        set(gca,'FontSize',14);
        axis([core_def(1,1),core_def(m,2),0,16]);
        title('Remaining coordinate deviation per residue');
    case 'theta'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coor0=coarse.Ca_coor;
        u0=ENM.u;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        coor=coarse.Ca_coor;
        rel_core=core_def-core_def(1,1)+1;
        range4=rel_core(4,1):rel_core(4,2);
        rangeb=[rel_core(1,1):rel_core(1,2),rel_core(2,1):rel_core(2,2),rel_core(6,1):rel_core(6,2),rel_core(7,1):rel_core(7,2)];
        rangeb2=[rel_core(1,1)+15:rel_core(1,2),rel_core(2,1):rel_core(2,2),rel_core(6,1):rel_core(6,2),rel_core(7,1):rel_core(7,2)];
        range5=rel_core(5,1):rel_core(5,2);
        n5=length(range5);
        n5h=floor(n5/2)-1;
        range5c=range5(1):range5(1)+n5h; % cytoplasmic half
        range5p=range5(end)-n5h:range5(end); % periplasmic half
        range10=rel_core(10,1):rel_core(10,2);
        n10=length(range10);
        n10h=floor(n10/2)-1;
        range10p=range10(1):range10(1)+n10h; % periplasmic half
        range10c=range10(end)-n10h:range10(end); % cytoplasmic half
        mid5=sum(coor0(range5,:),1)/length(range5);
        mid10=sum(coor0(range10,:),1)/length(range10);
        yz=mid10-mid5;
        [p0,h4vec]=rmsd_line_3D(coor0(range4,:));
        [p0,h5vec]=rmsd_line_3D(coor0(range5,:));
        [p0,h10vec]=rmsd_line_3D(coor0(range10,:));
        [p0h,hbvec]=rmsd_line_3D(coor0(rangeb,:));
        [p0,hb2vec]=rmsd_line_3D(coor0(rangeb2,:));
        [p0,h5pvec]=rmsd_line_3D(coor0(range5p,:));
        [p0,h5cvec]=rmsd_line_3D(coor0(range5c,:));
        [p0,h10pvec]=rmsd_line_3D(coor0(range10p,:));
        [p0,h10cvec]=rmsd_line_3D(coor0(range10c,:));
        h4vec=h4vec/norm(h4vec);
        h5vec=h5vec/norm(h5vec);
        h5pvec=h5pvec/norm(h5pvec);
        h5cvec=h5cvec/norm(h5cvec);
        h10vec=h10vec/norm(h10vec);
        h10pvec=h10pvec/norm(h10pvec);
        h10cvec=h10cvec/norm(h10cvec);
        hbvec=hbvec/norm(hbvec);        
        hb2vec=hb2vec/norm(hb2vec);        
        yz=yz/norm(yz);
        x=cross(yz,hb2vec);
        x=x/norm(x);
        y=cross(hb2vec,x);
        y=y/norm(y);
        evec=[x;y;hb2vec];
        ang0=180*acos(sum(h4vec.*hbvec))/pi;
        ang02=180*acos(sum(h4vec.*hb2vec))/pi;
        eb=hb2vec*evec';
        ex=x*evec';
        ey=y*evec';
        e4=h4vec*evec';
        e5=-h5vec*evec';
        e10=h10vec*evec';
        e5c=-h5cvec*evec';
        e10c=h10cvec*evec';
        e5p=-h5pvec*evec';
        e10p=h10pvec*evec';
        figure(2); clf;
        plot3([0,eb(1)],[0,eb(2)],[0,eb(3)],'k');
        hold on;
        plot3([0,ex(1)],[0,ex(2)],[0,ex(3)],'k:');
        plot3([0,ey(1)],[0,ey(2)],[0,ey(3)],'g:');
        plot3([0,e4(1)],[0,e4(2)],[0,e4(3)],'b');
        plot3([0,e5(1)],[0,e5(2)],[0,e5(3)],'m');
        plot3([0,e5c(1)],[0,e5c(2)],[0,e5c(3)],'m:');
        plot3([0,e5p(1)],[0,e5p(2)],[0,e5p(3)],'m-.');
        plot3([0,e10(1)],[0,e10(2)],[0,e10(3)],'r');
        plot3([0,e10c(1)],[0,e10c(2)],[0,e10c(3)],'r:');
        plot3([0,e10p(1)],[0,e10p(2)],[0,e10p(3)],'r-.');
        e5xy=e5;
        e5xy(3)=0;
        e5xy=e5xy/norm(e5xy);
        ang5=180*acos(e5xy(1))/pi;
        e5pxy=e5p;
        e5pxy(3)=0;
        e5pxy=e5pxy/norm(e5pxy);
        ang5p=180*acos(e5pxy(1))/pi;
        e5cxy=e5c;
        e5cxy(3)=0;
        e5cxy=e5cxy/norm(e5cxy);
        ang5c=180*acos(e5cxy(1))/pi;
        e10xy=e10;
        e10xy(3)=0;
        e10xy=e10xy/norm(e10xy);
        ang10=180*acos(e10xy(1))/pi;
        e10pxy=e10p;
        e10pxy(3)=0;
        e10pxy=e10pxy/norm(e10pxy);
        ang10p=180*acos(e10pxy(1))/pi;
        e10cxy=e10c;
        e10cxy(3)=0;
        e10cxy=e10cxy/norm(e10cxy);
        ang10c=180*acos(e10cxy(1))/pi;
        range5=rel_core(5,1):rel_core(5,2);
        range10=rel_core(10,1):rel_core(10,2);
        mid5=sum(coor(range5,:),1)/length(range5);
        mid10=sum(coor(range10,:),1)/length(range10);
        yz=mid10-mid5;
        [p0,h4vec]=rmsd_line_3D(coor(range4,:));
        [p0,h5vec]=rmsd_line_3D(coor(range5,:));
        [p0,h5cvec]=rmsd_line_3D(coor(range5c,:));
        [p0,h5pvec]=rmsd_line_3D(coor(range5p,:));
        [p0,h10vec]=rmsd_line_3D(coor(range10,:));
        [p0,h10cvec]=rmsd_line_3D(coor(range10c,:));
        [p0,h10pvec]=rmsd_line_3D(coor(range10p,:));
        [p0,hbvec]=rmsd_line_3D(coor(rangeb,:));
        [p0,hb2vec]=rmsd_line_3D(coor(rangeb2,:));
        h4vec=h4vec/norm(h4vec);
        h5vec=h5vec/norm(h5vec);
        h5cvec=h5cvec/norm(h5cvec);
        h5pvec=h5pvec/norm(h5pvec);
        h10vec=h10vec/norm(h10vec);
        h10cvec=h10cvec/norm(h10cvec);
        h10pvec=h10pvec/norm(h10pvec);
        hbvec=hbvec/norm(hbvec);        
        hb2vec=hb2vec/norm(hb2vec);        
        yz=yz/norm(yz);
        x=cross(yz,hb2vec);
        x=x/norm(x);
        y=cross(hb2vec,x);
        y=y/norm(y);
        evec=[x;y;hb2vec];
        eb=hb2vec*evec';
        ex=x*evec';
        ey=y*evec';
        e4=h4vec*evec';
        e5=-h5vec*evec';
        e5c=-h5cvec*evec';
        e5p=-h5pvec*evec';
        e10=h10vec*evec';
        e10c=h10cvec*evec';
        e10p=h10pvec*evec';
        e5xy=e5;
        e5xy(3)=0;
        e5xy=e5xy/norm(e5xy);
        ang5_2=180*acos(e5xy(1))/pi;
        e5cxy=e5c;
        e5cxy(3)=0;
        e5cxy=e5cxy/norm(e5cxy);
        ang5c_2=180*acos(e5cxy(1))/pi;
        e5pxy=e5p;
        e5pxy(3)=0;
        e5pxy=e5pxy/norm(e5pxy);
        ang5p_2=180*acos(e5pxy(1))/pi;
        e10xy=e10;
        e10xy(3)=0;
        e10xy=e10xy/norm(e10xy);
        ang10_2=180*acos(e10xy(1))/pi;
        e10cxy=e10c;
        e10cxy(3)=0;
        e10cxy=e10cxy/norm(e10cxy);
        ang10c_2=180*acos(e10cxy(1))/pi;
        e10pxy=e10p;
        e10pxy(3)=0;
        e10pxy=e10pxy/norm(e10pxy);
        ang10p_2=180*acos(e10pxy(1))/pi;
        ang=180*acos(sum(h4vec.*hbvec))/pi;
        ang2=180*acos(sum(h4vec.*hb2vec))/pi;
        fprintf(1,'For state %s, theta is %4.1f degree\n',state1,ang0); 
        fprintf(1,'For state %s, theta(4;1b,2,6,7) is %4.1f degree\n',state1,ang02); 
        fprintf(1,'For state %s, fi5 is %4.1f degree\n',state1,ang5); 
        fprintf(1,'For state %s, fi5c is %4.1f degree\n',state1,ang5c); 
        fprintf(1,'For state %s, fi5p is %4.1f degree\n',state1,ang5p); 
        fprintf(1,'For state %s, fi10 is %4.1f degree\n',state1,ang10); 
        fprintf(1,'For state %s, fi10c is %4.1f degree\n',state1,ang10c); 
        fprintf(1,'For state %s, fi10p is %4.1f degree\n',state1,ang10p); 
        fprintf(1,'For state %s, theta is %4.1f degree\n',state2,ang); 
        fprintf(1,'For state %s, theta(4;1b,2,6,7) is %4.1f degree\n',state2,ang2); 
        fprintf(1,'For state %s, fi5 is %4.1f degree\n',state2,ang5_2); 
        fprintf(1,'For state %s, fi5c is %4.1f degree\n',state2,ang5c_2); 
        fprintf(1,'For state %s, fi5p is %4.1f degree\n',state2,ang5p_2); 
        fprintf(1,'For state %s, fi10 is %4.1f degree\n',state2,ang10_2); 
        fprintf(1,'For state %s, fi10c is %4.1f degree\n',state2,ang10c_2); 
        fprintf(1,'For state %s, fi10p is %4.1f degree\n',state2,ang10p_2); 
    case 'kinks'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coor0=coarse.Ca_coor;
        u0=ENM.u;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        coor=coarse.Ca_coor;
        rel_core=core_def-core_def(1,1)+1;
        kink1=get_kink(rel_core,coor0,1);
        kink5=get_kink(rel_core,coor0,5);
        kink6=get_kink(rel_core,coor0,6);
        kink10=get_kink(rel_core,coor0,10);
        kink1_2=get_kink(rel_core,coor,1);
        kink5_2=get_kink(rel_core,coor,5);
        kink6_2=get_kink(rel_core,coor,6);
        kink10_2=get_kink(rel_core,coor,10);
        fprintf(1,'For state %s, kink of TMD 1  is %4.1f degree\n',state1,kink1); 
        fprintf(1,'For state %s, kink of TMD 5  is %4.1f degree\n',state1,kink5); 
        fprintf(1,'For state %s, kink of TMD 6  is %4.1f degree\n',state1,kink6); 
        fprintf(1,'For state %s, kink of TMD 10 is %4.1f degree\n',state1,kink10); 
        fprintf(1,'For state %s, kink of TMD 1  is %4.1f degree\n',state2,kink1_2); 
        fprintf(1,'For state %s, kink of TMD 5  is %4.1f degree\n',state2,kink5_2); 
        fprintf(1,'For state %s, kink of TMD 6  is %4.1f degree\n',state2,kink6_2); 
        fprintf(1,'For state %s, kink of TMD 10 is %4.1f degree\n',state2,kink10_2); 
    case 'constraints'
        infile1=sprintf('%s_ANM.mat',state1);
        load(infile1);
        coarse0=coarse;
        infile2=sprintf('%s_ANM.mat',state2);
        load(infile2);
        rel_core=core_def-core_def(1,1)+1; % relative core definition
        core_ext=[];
        core2_ext=[];
        displacements=coarse.Ca_coor-coarse0.Ca_coor;
        [mc,nc]=size(core_def2);
        for k=1:mc,
            core2_ext=[core2_ext core_def2(k,1):core_def2(k,2)];
        end;
        [mc,nc]=size(rel_core);
        for k=1:mc,
            core_ext=[core_ext rel_core(k,1):rel_core(k,2)];
        end;
        infile3=sprintf('%s_ANM.mat',template);
        load(infile3);
        tcoarse=coarse;
        diff=tcoarse.Ca_coor-coarse.Ca_coor;
        weights=sum(diff.^2,2);
        weights=ones(size(weights))./weights;
        weights=weights/sum(weights);
        figure(1); clf;
        plot(weights,'k.');
        if strcmpi(template,state1),
            displace_mode='uniform';
        else
            displace_mode='uniform';
        end;
        if dali_lite,
            outfile=sprintf('%s_to_%s_displacements_for_%s_to_%s.dat',state1,state2,template,target);
        else
            outfile=sprintf('%s_to_%s_displacements_for_%s_to_%s.dat',state1,state2,template,target);
        end;
        fid=fopen(outfile,'wt');
        if ~isempty(target),
            fprintf(fid,'%% MMM simulated displacement constraints for target structure [%s]\n',target);
            fprintf(fid,'# TARGET %s\n',target);
        else
            fprintf(fid,'%% MMM simulated displacement constraints for template structure [%s]\n',template);
        end;
        fprintf(fid,'# PDB %s\n',template);
        fprintf(fid,'# TRANSITION %s %s\n',state1,state2);
        fprintf(fid,'# DISPLACEMENTS %s\n',displace_mode);
        for k=1:length(core_ext),
            cdisp=displacements(core_ext(k),:);
            resnum=core2_ext(k);
            fprintf(fid,'%10i  %6.2f %6.2f %6.2f\n',resnum,cdisp);
        end;
        fprintf(fid,'# END\n');
        fclose(fid);
end;

function [rmsd,coor1] = fit_by_ANM(coor0,target,basis)

coor=coor0;
iter=0;
rmsd=rmsd_superimpose(target,coor);

figure(1); clf;
set(gca,'FontSize',14);
title('Deviation (C^\alpha r.m.s.d.) from target structure');
h=line(iter,rmsd,'Color','k');

rmsd_diff=rmsd;
nit=0;

while nit<400 && rmsd_diff > 5e-6,
    nit=nit+1;
    u=get_ANM(coor);
    diff=target-coor;
    [m,m1]=size(u);
    [n,n1]=size(diff);
    change=reshape(diff',1,3*n);
    coeff=u(:,1:basis)\change';
    step=zeros(n,3);
    for k=1:basis,
        evec=u(:,k);
        mode1=reshape(evec,3,m/3);
        step=step+mode1'*coeff(k);
    end;
    per_residue=sqrt(sum(step.^2,2));
    sc=0.2/max(per_residue);
    if sc>1, sc=1; end;
    coor=coor+sc*step;
    crmsd=rmsd_superimpose(target,coor);
    iter=[iter nit];
    rmsd=[rmsd crmsd];
    if nit>10,
        rmsd_diff=(rmsd(end-10)-rmsd(end))/10;
    end;
    set(h,'XData',iter,'YData',rmsd);
    drawnow;
end;
axis([0,max(iter),0,1.1*max(rmsd)]);
coor1=coor;



function [u,lambda]=get_ANM(coor)

Hessian=setup_ANM_bonded(coor);
[u,d]=eig(Hessian);
clear Hessian
lambda=diag(d);

function Hessian=setup_ANM_bonded(Ca_coor)
% function Hessian=setup_ANM_bonded(Ca_coor)
%
% returns the Hessian for a fully connected anisotropic network model (ANM)
% when given the Calpha coordinates
% the force constants have inverse 6th power dependence on interresidue
% distance, except for next neighbors (i,i+1) and second neighbors (i,i+2)
% along the polypetide chain
%
% ### This version assumes a single-chain peptide without gaps in the sequence ###
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% and:
% L. Yang, G. Song, R L. Jernigan, Proc. Natl. Acad. Sci. 106 (2009)
% 12347–12352.
%
% the force constant at a distance of 1 Å p_ANM_gamma is defined in
% initialize_MMM.m  and stored in global variable ENM_param
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010

global ENM_param

bond_force_1=10000*3.8^(-6); % 1200
bond_force_2=10000*6^(-6); % 80
exponent=6;

[m,n]=size(Ca_coor); % m is number of Calpha atoms

Hessian=zeros(3*m);
r1=zeros(1,m);
p1=0;
r2=zeros(1,m);
p2=0;
for i=1:m-1,
    basi=3*(i-1);
    for j=i+1:m,
        basj=3*(j-1);
        rvec=Ca_coor(i,:)-Ca_coor(j,:); % distance vector of Calpha atoms
        r=norm(rvec);
        diff=i-j;
        switch diff
            case 1
                cgamma=-bond_force_1*ENM_param.p_ANM_gamma;
                p1=p1+1;
                r1(p1)=r;
            case 2
                cgamma=-bond_force_2*ENM_param.p_ANM_gamma;
                p2=p2+1;
                r2(p2)=r;
            otherwise
                cgamma=-ENM_param.p_ANM_gamma*r^(-exponent);
        end;
        if cgamma~=0,
            submat=kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
            Hessian(basi+1:basi+3,basj+1:basj+3)=cgamma*submat;
            Hessian(basj+1:basj+3,basi+1:basi+3)=cgamma*submat;
        end;
    end;
end;

conn=0;
for i=1:m,
    basi=3*(i-1);
    submat=zeros(3);
    for j=1:m,
        basj=3*(j-1);
        submat=submat+Hessian(basi+1:basi+3,basj+1:basj+3);
    end;
    conn=conn+trace(submat);
    Hessian(basi+1:basi+3,basi+1:basi+3)=-submat;
end;

if conn==0, % ANM is unconnected
    Hessian=[];
end;

function kink = get_kink(rel_core,coor0,n)

range=rel_core(n,1):rel_core(n,2);
ll=length(range);
lh=floor(ll/2)-1;
rangec=range(1):range(1)+lh; % cytoplasmic half
rangep=range(end)-lh:range(end); % periplasmic half
[p0,cvec]=rmsd_line_3D(coor0(rangec,:));
cvec=cvec/norm(cvec);
[p0,hvec]=rmsd_line_3D(coor0(rangep,:));
hvec=hvec/norm(hvec);

kink=180*acos(sum(cvec.*hvec))/pi;
if kink>90, kink=180-kink; end;
