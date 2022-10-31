%% Figure 5
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
h3cc=load('./H3CC_henikoff.mat');
h3smooth=conv(mean(h3cc.p_all,2),normpdf([-50:50]',0,25),'same');
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
load('SC_genome.mat')
GP=load('./group_imp.mat')
load('./promCorrSort.mat')

%% Figure 5D - Examples: signal on idividual promoters, sumProm scatters and binding signal corr matrix
examplePairs = {'Stp2','Stp1'; 'Pdr3', 'Pdr1'; 'Dal80', 'Gzf3'};
targets = {'BAP3'; 'SNQ2'; 'GAP1'};
motifList = {'ACGGC', 'CGGAA', 'GATAA'};
zscoreTH = 3.5;
quantileTH = 0.99;
Xlimit = [-4 4];
gb = true(sum(GP.chr_len(1:16)), 1);
gb(GP.chrIdx(12)+[451575-51:468958+51]) = false;
yLimitSignal = 0.5;
intBases = intBasesVec(1:6701);
intBasesMotifs = createIntBasesForMotif();
intoGene = 150;
promoterL = repmat(700,6701,1);

lineCol = lines(2);
col = brewermap(1,'Dark2');%cbrewer2('Dark2',1);
cMapCorr = flipud(bone);
cMapScatter = brewermap(256,'OrRd');%cbrewer2('OrRd');
markerCol = cMapScatter(128,:);
colSignal = repmat([0 0 0; 0.5 0.5 0.5], 2,1);
colMotifs = lines(6);

yspacer = 0.05;
xspacer = 0.06;
Wsignal = 0.25;
xPos = [0.05: Wsignal+xspacer*1.3: 1];
Hsignal = 0.07;
Wsignal = 0.25;
Wscatter = 0.62*Wsignal-0.03;
Wcorrmat = Wsignal-Wscatter-xspacer;
Hscatter =Wscatter*1.8;
yPos = 0.5;
figure('Units','pixels','Position',[141,19,1732.00000000000,962], 'Color','w')
for i = 1:size(examplePairs,1)
    targetGeneidx = GP.gene_table.(targets{i,1});
    selTargets = summaryTable.nHighTargetsIdx{contains(summaryTable.label, examplePairs{i,1})};
    clear logChange significantGenes zScoreTF
    for tf = 1:2
        currF1 = examplePairs{i,tf};
        m1 = [currF1,'_d',upper(examplePairs{i,3-tf})];
        zscore_WT1 = nanZscore(checWTdelLactisSwap.sumProm.(currF1));
        fitGenes =  zscore_WT1 > min(zscoreTH, quantile(zscore_WT1, quantileTH));
        bestFit = robustfit(checWTdelLactisSwap.sumProm.(currF1)(fitGenes), checWTdelLactisSwap.sumProm.(m1)(fitGenes), [], [],'off' );
        normFactor(tf) = log2(bestFit);
        logChange(tf,:) = log2(checWTdelLactisSwap.sumProm.(m1)+700)- log2(checWTdelLactisSwap.sumProm.(currF1)+700)- normFactor(tf);
        significantGenes(tf,:) =  abs((logChange(tf,:)-median(logChange(tf,fitGenes)))/std(logChange(tf,fitGenes))) >= 1;
    end
    
    % scatter paralog that changed
    axes('Position', [xPos(i) yPos Wscatter Hscatter])
    maxX = max(checWTdelLactisSwap.sumProm.(examplePairs{i,1}));
    maxY = max(checWTdelLactisSwap.sumProm.([examplePairs{i,1}, '_d', upper(examplePairs{i,2})]));
    maxCol = max(checWTdelLactisSwap.sumProm.(examplePairs{i,2}));

    scatter(checWTdelLactisSwap.sumProm.(examplePairs{i,1})/maxX,...
        checWTdelLactisSwap.sumProm.([examplePairs{i,1}, '_d', upper(examplePairs{i,2})])/maxY,[],...
        checWTdelLactisSwap.sumProm.(examplePairs{i,2})/maxCol, 'filled','MarkerEdgeColor', 'k')
    hold on
    plot(xlim, xlim*2^(normFactor(1))*maxX/maxY, 'k:')
    text(checWTdelLactisSwap.sumProm.(examplePairs{i,1})(targetGeneidx)/maxX-0.02*max(xlim),...
        checWTdelLactisSwap.sumProm.([examplePairs{i,1}, '_d', upper(examplePairs{i,2})])(targetGeneidx)/maxY, targets{i}, 'HorizontalAlignment','right')
    
    xlabel(examplePairs{i,1}, 'fontSize',15)
    ylabel([examplePairs{i,1}, ' Delta', upper(examplePairs{i,2})], 'fontSize',15)
    text(0.015, 1,...
        sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(examplePairs{i,1}),...
        checWTdelLactisSwap.sumProm.([examplePairs{i,1}, '_d', upper(examplePairs{i,2})]), 'rows','pairwise')), 'fontSize',15, 'HorizontalAlignment', 'left', 'VerticalAlignment','top')
    ylim([0 1])
    xlim([0 1])
    set(gca, 'XTick', [0:0.5:1], 'YTick', [0:0.5:1])
    cbr = colorbar()
    ylabel(cbr, examplePairs{i,2}, 'fontSize',12)
    %cbr.Position = cbr.Position.*[1 1 0.5 1];
    caxis([0 0.8*max(caxis)])
    colormap(gca, cMapScatter)
    cbr.Position = [0.1797-xPos(1)+xPos(i) yPos 0.005 Hscatter]
    set(cbr, 'Ticks', [0:0.4:0.8])
%     setAxisExponent('axes',cbr,'axis','c')
%     if numel(cbr.Ticks)>5
%         cbr.Ticks = cbr.Ticks(1:2:end);
%     end
    
    % scatter other paralog
    axes('Position', [xPos(i)+Wscatter+xspacer yPos Wcorrmat Wcorrmat*1.3])
    maxX = max(checWTdelLactisSwap.sumProm.(examplePairs{i,2}));
    maxY = max(checWTdelLactisSwap.sumProm.([examplePairs{i,2}, '_d', upper(examplePairs{i,1})]));
    maxCol = max(checWTdelLactisSwap.sumProm.(examplePairs{i,1}));

    scatter(checWTdelLactisSwap.sumProm.(examplePairs{i,2})/maxX,...
        checWTdelLactisSwap.sumProm.([examplePairs{i,2}, '_d', upper(examplePairs{i,1})])/maxY,[],...
        checWTdelLactisSwap.sumProm.(examplePairs{i,1})/maxCol, 'filled','MarkerEdgeColor', 'k')
    hold on
    plot(xlim, xlim*2^(normFactor(2))*maxX/maxY, 'k:')
    xlabel(examplePairs{i,2}, 'fontSize',12)
    ylabel([examplePairs{i,2}, ' Delta', upper(examplePairs{i,1})], 'fontSize',12)
    text(0.02,1,...
        sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(examplePairs{i,2}),...
        checWTdelLactisSwap.sumProm.([examplePairs{i,2}, '_d', upper(examplePairs{i,1})]), 'rows','pairwise')), 'fontSize',10,  'HorizontalAlignment', 'left', 'VerticalAlignment','top')
    ylim([0 1])
    xlim([0 1])
    cbr = colorbar()
    ylabel(cbr, examplePairs{i,1}, 'fontSize',12)
    caxis([0 0.8*max(caxis)])
    colormap(gca, cMapScatter)
    set(cbr, 'Ticks', [0:0.4:0.8])
    
    % corr matrix
    axes('Position', [xPos(i)+Wscatter+xspacer yPos+yspacer+Wcorrmat*1.3 Wcorrmat Wcorrmat*1.3])
    intStrains = {examplePairs{i,2}, [examplePairs{i,2}, '_d', upper(examplePairs{i,1})], examplePairs{i,1}, [examplePairs{i,1}, '_d', upper(examplePairs{i,2})]};
    % WITH REPEATS
    %     [corrSumProm, idxSumProm] = getRepeatsCorr(intStrains, 'dataType','sumProm');
    %     [corr7mer, ~] = getRepeatsCorr(intStrains, 'dataType','7mer');
    %
    %     combineMat = tril(corrSumProm) + triu(corr7mer,1);
    %     imagesc(combineMat)
    %     borders = [0; cumsum(accumarray(idxSumProm,1))]+0.5;
    %     tickPos = movmean(borders,2,'Endpoints','discard');
    %
    %     hold on
    %     plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
    %     plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
    
    % WITHOUT REPEATS
    corrSumProm = plotSumPromCorr(intStrains, checWTdelLactisSwap,0);
    normProfile = chromosomes2fullProfile(checWTdelLactisSwap, intStrains);
    mers = mer_occupancy(normProfile, 7, 'method','else','intBases',intBasesMotifs);
    corr7mer = corr(mers.score, 'rows','pairwise');
    combinedMat = tril(corrSumProm) + triu(corr7mer,1);
    imagesc(combinedMat)
    plotgrid(combinedMat)
    plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)

    colormap(gca, cMapCorr)
    caxis([0.5 1])
    set(gca, 'YTick', [1:numel(intStrains)], 'YTickLabel',{[examplePairs{i,2},' (P1)'],'P1DeltaP2',[examplePairs{i,1},' (P2)'],'P2DeltaP1'}, 'XTick',[],'fontSize',9)
    hold on
    plot(xlim, ylim, 'color',[0.7 0.7 0.7])
    cbr = colorbar()
    ylabel(cbr,'correlation', 'fontSize',15)
    set(cbr, 'Ticks', [0.6:0.2:1])
            
    % log2 change
    axes('Position', [xPos(i) yPos-2*yspacer-Hsignal/2 Wsignal Hsignal/2])
    for tf =1:2
        wtVec(tf,:) = checWTdelLactisSwap.sumProm.(examplePairs{i,tf});
    end
    for tf = 1:2
        plot([Xlimit], [tf tf], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
        hold on
        grayIdx = selTargets(~significantGenes(tf,selTargets));
        significantIdx = selTargets(significantGenes(tf,selTargets));
        sizeVec = rescale(checWTdelLactisSwap.sumProm.(examplePairs{i,tf})(selTargets),20,100);
        % color rescaled over all genes
        fullColVec =  rescale(wtVec(3-tf,:),0,1, 'InputMax', 0.8*max(wtVec(3-tf,:)), 'InputMin',0);
        colVec = fullColVec(selTargets);
        % color rescaled only to selected targets
        %colVec = rescale(checWTdelLactisSwap.sumProm.(examplePairs{i,3-tf})(selTargets),0,1, 'InputMax', 0.8*max(wtVec(3-tf,:)), 'InputMin',0);
        scatter(logChange(tf,grayIdx), tf*ones(numel(grayIdx),1),...
            sizeVec(ismember(selTargets,grayIdx)), [0.7 0.7 0.7],'filled',...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0.7 0.7 0.7]);
        scatter(logChange(tf, significantIdx), tf*ones(numel(significantIdx),1), sizeVec(ismember(selTargets,significantIdx)),...
           colVec(ismember(selTargets,significantIdx)),...
            'filled','MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cMapScatter(end,:));
    end
    %set(gca, 'YDir', 'reverse')
    ylim([0.5 2.5])
    plot([0 0],ylim, 'k--')
    colormap(gca, cMapScatter)
    caxis([0 1])
    title('Binding to top targets after paralog deletion', 'fontSize',12, 'FontWeight','normal')
    set(gca, 'YTick', [1,2], 'YTickLabel', examplePairs(i,:), 'fontSize',12, 'XTick', [-4:2:4], 'TickLength',[0 0])
    xlabel('Fold change (log2)', 'fontSize',11)
    
    target = targets{i};
    geneIdx = GP.gene_table.(target);
    plot(logChange(1, geneIdx)*[1,1], [1.5 2.2], 'color','k') 
    scatter(logChange(1, geneIdx), 1.65,20, 'vk','filled')
    if i == 2
        text(logChange(1, geneIdx)+0.05*logChange(1, geneIdx), 2.3, target,'fontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','top')
    else
        text(logChange(1, geneIdx)+0.05*logChange(1, geneIdx), 2.3, target,'fontSize',8, 'HorizontalAlignment','right', 'VerticalAlignment','top')
    end
    
    % colorbar log2
    axes('Position', [xPos(i)+0.007 yPos-2*yspacer-Hsignal/2-0.05 0.05 0.01])
    colormap(gca, cMapScatter)
    cbr = colorbar('Location','south');
    cbr.AxisLocation =  'out';
%     set(cbr, 'Ticks', [min(caxis), max(caxis)], 'TickLabels', {'min','max'}, 'TickLength', [0 0], 'FontSize',7)
    set(cbr, 'Ticks', [min(caxis), max(caxis)], 'TickLabels', {'0','0.8'}, 'TickLength', [0 0], 'FontSize',8)
    title(cbr, 'Paralogs signal', 'fontSize',10)
    axis off
    
    % promoter
    target = targets{i};
    geneIdx = GP.gene_table.(target);
    pos = GP.gene_infoR64.position(geneIdx,:);
    intBasesGene = GP.chr_idx(pos(1))+pos(2)+[-promoterL(geneIdx):intoGene].*GP.gene_infoR64.dir(geneIdx);
    PromSeq = upper(findPromoterSeq(geneIdx, GP, promoterL,'tss','orf', 'intoGene', intoGene));
    chr = pos(1);
    if contains(motifList{i}, {'|','[','.'})
        regExpPattern = motifList{i};
    else
        regExpPattern = [motifList{i}, '|' seqrcomplement(motifList{i})];
    end
    [S,E] = regexp(PromSeq, regExpPattern,'start','end');
    clear motifPos
    motifPos.S = S-promoterL(geneIdx);
    motifPos.E = E-promoterL(geneIdx);

    axes('Position', [xPos(i) yPos-3.8*yspacer-Hsignal*1.8 Wsignal Hsignal*1.5])
    intStrainsS = {examplePairs{i,2}, [examplePairs{i,2},'_d',upper(examplePairs{i,1})],...
        examplePairs{i,1}, [examplePairs{i,1},'_d',upper(examplePairs{i,2})]};
    for z = 1:numel(intStrainsS)
        currNorm = chromosomes2fullProfile(checWTdelLactisSwap, intStrainsS(z));
        maxNorm = max(currNorm(gb));
        currProfile = currNorm(intBasesGene)./maxNorm;
        plot([-promoterL(geneIdx):intoGene], rescale(currProfile, numel(intStrainsS)-z+1.5, numel(intStrainsS)-z+2.3,'InputMin',0, 'InputMax',yLimitSignal),...
            'color', colSignal(z,:))
        hold on
    end
    plot(-promoterL(geneIdx):intoGene, rescale(h3smooth(intBasesGene), 0.5, 1.3),'k', 'LineWidth',1)
    area(-promoterL(geneIdx):intoGene, rescale(h3smooth(intBasesGene), 0.5, 1.3),'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'LineStyle', 'none')
    ylim([0.5 numel(intStrainsS)+1.5])
    xlim([-promoterL(geneIdx) intoGene])
    TssLoc = (GP.gene_infoR64.tamarTss(geneIdx,2)-pos(2))*GP.gene_infoR64.dir(geneIdx);
    %yLabels = strrep(intStrainsS([numel(intStrainsS):-1:1]), '_d', ' \Delta');
    yLabels = {'P2DeltaP1', [examplePairs{i,1},' (P2)'],'P1P2',[examplePairs{i,2},' (P1)']};
    set(gca, 'YTick', [1:numel(intStrainsS)+1], 'YTickLabel', {'NucOcc', yLabels{:}},...
        'TickLength',[0 0], 'XTick',[], 'fontSize',11)
    plot(TssLoc*[1 1] ,ylim, ':r', 'LineWidth',1.5)
    scatter(0 ,0.5,[], '>r','filled')
    if numel(motifPos.S) > 0
        for p = 1:numel(motifPos.S)
            fill([motifPos.S(p)*[1,1], motifPos.E(p)*[1,1]], [0.5, numel(intStrainsS)+1.5.*[1,1], 0.5], colMotifs(1,:),...
                'LineStyle','none', 'FaceAlpha', 0.5)
        end
    end
    title(sprintf('%s promter', target),'fontSize',9)
    box on
    text(min(xlim), 0.5, ['in-vitro motif: ', motifList{i}],'Color', colMotifs(1,:),'VerticalAlignment','top', 'fontSize',12)
    text(TssLoc,0.5, 'TSS', 'Color', 'r','VerticalAlignment','top', 'HorizontalAlignment', 'center', 'fontSize',10)

    targetGenePos = GP.gene_infoR64.position(targetGeneidx,:);
    targetGeneTss = GP.gene_infoR64.tamarTss(targetGeneidx,:);
    PromSeq = upper(findPromoterSeq(targetGeneidx, GP, promoterL,'tss','orf', 'intoGene', intoGene));
    chr = targetGenePos(1);
    relPromoterPos =  targetGenePos(2)+[-promoterL(targetGeneidx):intoGene-1]*GP.gene_infoR64.dir(targetGeneidx);    
end
%save_gf(gcf,'Fig3Exam','paper','tamar21','type','svg','size',[])


 %% Figure 5C - summary plot (scatters by groups)
 groups = fieldnames(promCorrSort);
 familyNames = {'Zinc finger', 'Zinc cluster, bZIP', 'others'};
 groupsIdx = {1, [2,3], 4};
 TFsChanged.ZFs(:,1) = {'Rph1', 'Stp4', 'Stp2','Mig2'};
 TFsChanged.ZCs_bZs(:,1) = {'Pdr3', 'Ecm22', 'Yap4','Yrr1','Tbs1'};
 TFsChanged.others(:,1) = {'Fkh1', 'Smp1', 'Dal80'};
 ZCidx = [1:6];
 bZIPidx = [7:10];
 cMap = brighten(flipud(bone),0.4);
 c = 1;
 figure
 for g =  1: length(groupsIdx)
     axes('Position', [-0.06+(g-1)*0.305 0.07 0.45 0.45]);
     nParalogs =  sum(ismember(summaryTable.familyId, groupsIdx{g}));
     c = g;
     corrWTs = nan(nParalogs,1);
     corrM1W1= nan(nParalogs,1);
     corrM2W2= nan(nParalogs,1);
     corrM1W2= nan(nParalogs,1);
     corrM2W1= nan(nParalogs,1);
     x1neg = [];
     x2neg=[];
     y1neg = [];
     y2neg = [];
     
     t = {};
     for i = 1:nParalogs
         currF1 = promCorrSort.(groups{g}).TFsSort{i,1};
         currF2 = promCorrSort.(groups{g}).TFsSort{i,2};
         m1 = [currF1, '_d', upper(currF2)];
         m2 = [currF2, '_d', upper(currF1)];
         if isfield(checWTdelLactisSwap.sumProm, m1)
             %corrM1W1(i,1) = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(currF1), 'rows','pairwise');
             %corrM1W2(i,1) = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(currF2), 'rows','pairwise');
             [corrM1W1(i,1), ~] = corrBetweenRepeatsOf2Samples({currF1,m1},'dataType', 'sumProm');
             [corrM1W2(i,1), y1neg(i,1)] = corrBetweenRepeatsOf2Samples({currF2,m1},'dataType', 'sumProm');
         end
         if isfield(checWTdelLactisSwap.sumProm, m2)
%              corrM2W1(i,1) = corr(checWTdelLactisSwap.sumProm.(m2), checWTdelLactisSwap.sumProm.(currF1), 'rows','pairwise');
%              corrM2W2(i,1) = corr(checWTdelLactisSwap.sumProm.(m2), checWTdelLactisSwap.sumProm.(currF2), 'rows','pairwise');
                [corrM2W1(i,1), y2neg(i,1)] = corrBetweenRepeatsOf2Samples({currF1,m2},'dataType', 'sumProm');
                [corrM2W2(i,1), ~] = corrBetweenRepeatsOf2Samples({currF2,m2},'dataType', 'sumProm');
         end         
         corrWTs(i,1) = corr(checWTdelLactisSwap.sumProm.(currF1), checWTdelLactisSwap.sumProm.(currF2), 'rows','pairwise');
         t{i,1} = currF1;
         t{i,2} = currF2;
     end
     hold on
     for p =1:nParalogs
         plot([corrWTs(p), corrWTs(p)], [corrM1W2(p), corrM2W1(p)], 'Color', [0.7 0.7 0.7])
         hold on
     end
     y1neg(y1neg==0) = nan;
     y2neg(y2neg==0) = nan;
     
     if length(groupsIdx{g}) > 1
              scatter([corrWTs(ZCidx);corrWTs(ZCidx)], [corrM1W2(ZCidx); corrM2W1(ZCidx)], 85, [corrM1W1(ZCidx);corrM2W2(ZCidx)],...
                  'filled', 'MarkerEdgeColor','r', 'LineWidth',1.5);
              hold on
               scatter([corrWTs(bZIPidx);corrWTs(bZIPidx)], [corrM1W2(bZIPidx); corrM2W1(bZIPidx)], 85, [corrM1W1(bZIPidx);corrM2W2(bZIPidx)],...
                  'filled', 'MarkerEdgeColor','k', 'LineWidth',1.5);
     else        
        scatter([corrWTs;corrWTs], [corrM1W2; corrM2W1], 85, [corrM1W1;corrM2W2], 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.5);
        hold on
     end
     errorbar([corrWTs;corrWTs], [corrM1W2; corrM2W1], [y1neg; y2neg], [y1neg; y2neg], 'LineStyle','none', 'color','k', 'lineWidth', 1.5)
     axis square
     box on
     xlim([0 1]);
     ylim([0 1]);
     caxis([0.6 1])
     xlabel('WT-paralog correlation',  'FontSize',18);
     ylabel('mutant-paralog correlation',  'FontSize',18);
     hold on
     plot([1 0],[1 0], '--k','LineWidth',1);
     
     colormap(gca, cMap);
     set(gca, 'XTick', [0:0.2:1], 'YTick', [0:0.2:1],'fontSize',15)
     
     if c == 3;
         cbr = colorbar('Location', 'east')
         currCbrPos = cbr.Position;
         set(cbr, 'position', [0.91    0.0687    0.0125    0.4524]);
         set(cbr, 'AxisLocation', 'out')
         ylabel(cbr, 'mutant-WT correlation' , 'FontSize',18)
         set(cbr, 'YTick',[0.6:0.1:1])
         caxis([0.6 1])
     end
     textTFs(1:length(t)) = t(:,1);
     textTFs(length(t)+1:length(t)*2) = t(:,2);
     for tf = 1: length(TFsChanged.(groups{g}))
         tfIdx = find(strcmp(textTFs,TFsChanged.(groups{g}){tf}));
         if tfIdx < size(t,1)+1
             text([corrWTs(tfIdx)]+0.015, corrM1W2(tfIdx), t(tfIdx,1), 'FontSize', 15, 'HorizontalAlignment', 'left');
         else
             text([corrWTs(tfIdx-size(t,1))]+0.015, corrM2W1(tfIdx-size(t,1)), t(tfIdx-size(t,1),2), 'FontSize', 15, 'HorizontalAlignment', 'left');
         end
     end
     title(familyNames{g}, 'FontSize',20)
 end
 set(gcf, 'Color','w');
% save_gf(gcf,'/home/labs/barkailab/tamarge/Master/paperIllustrator/Rev/Fig3Sum','paper','tamar21','type','svg')
% saveas(gcf, 'Fig3SummaryPlotWithSTD.svg')


%% Figure 5A - scheme explaining summary figure
colVec = lines(2);
figure('Units','normalized','Position', [0 0 0.5 1])
w = 0.13;
axes('Position', [0.05 0.5 w w*1.996])

plot([0 1], [0 1], 'k--', 'lineWidth',2)
hold on
xlabel('WT-paralog correlation',  'FontSize',18);
ylabel('mutant-paralog correlation',  'FontSize',18);
set(gca, 'XTick', [0:0.2:1], 'YTick', [0:0.2:1], 'fontSize',12)
set(gcf, 'color','w')
fill([0 1 1], [0 0 1], colVec(2,:), 'FaceAlpha', 0.05, 'LineStyle','none')
fill([0 1 0], [0 1 1], colVec(1,:), 'FaceAlpha', 0.05, 'LineStyle','none')
text(0.8, 0.8, 'No change', 'Rotation',45,'HorizontalAlignment', 'center','VerticalAlignment','bottom', 'fontSize',12)
plot([0.55 0.55], [0.55 0.75],'color',[0.7 0.7 0.7], 'LineWidth',1)
%scatter(0.3, 0.55, 'k^','filled')
scatter(0.55, 0.75, 50,'k', 'filled')
text(0.45, 0.775, 'Gained paralog''s\newline       targets', 'HorizontalAlignment','center','VerticalAlignment','bottom')

plot([0.55 0.55], [0.4 0.55],'color',[0.7 0.7 0.7], 'LineWidth',1)
%scatter(0.8, 0.4, 'kv','filled')
scatter(0.55, 0.4,50, 'k', 'filled')
text(0.55, 0.35, 'Lost common\newline     targets', 'HorizontalAlignment','center','VerticalAlignment','top')

text(0.05, 0.95, 'Robustness', 'color',colVec(1,:), 'FontWeight','bold', 'fontSize',12)
text(0.8, 0.05, 'Fragility', 'color',colVec(2,:), 'FontWeight','bold', 'fontSize',12, 'HorizontalAlignment','center')

scatter(0.65, 0.65, 'k', 'filled')
%save_gf(gcf,'Fig3Scheme','paper','tamar21','type','svg')
%% bars correlations with WTs - with Names
[~,maxCorrIdx] = max(summaryTable.WTmutantCorr, [], 2);
 familyNames = {'Zinc finger', 'Zinc cluster, bZIP', 'others'};

% yData_old = summaryTable.WTmutantCorr(repmat(1,30,1).*[1:30]'+(maxCorrIdx-1)*30);
% xData_old = summaryTable.WTmutantCorr(repmat(1,30,1).*[1:30]'+(3-maxCorrIdx-1)*30);

yData = nan(size(summaryTable,1),1);
xData = nan(size(summaryTable,1),1);
yErr = nan(size(summaryTable,1),1);
xErr = nan(size(summaryTable,1),1);
for i = 1:size(summaryTable,1)
    p1 = sprintf('p%d',maxCorrIdx(i));
    p2 = sprintf('p%d',3-maxCorrIdx(i));
    currF1 = summaryTable.(p1){i};
    currF2 = summaryTable.(p2){i};
    m1 = [currF1, '_d', upper(currF2)];
    m2 = [currF2, '_d', upper(currF1)];
    if isfield(checWTdelLactisSwap.sumProm, m1)
        [yData(i,1), yErr(i,1)] = corrBetweenRepeatsOf2Samples({currF1,m1},'dataType', 'sumProm');
    end
    if isfield(checWTdelLactisSwap.sumProm, m2)
        [xData(i,1), xErr(i,1)] = corrBetweenRepeatsOf2Samples({currF2,m2},'dataType', 'sumProm');
    end
end

col = [0.1 0.1 0.1; 0.7 0.7 0.7];

figure('Units','normalized','Position', [0 0 0.5 1], 'color','w')
for i = 1:3
    axes('Position', [0.1+(i-1)*0.11 0.2 0.1 0.2])
    barh([1:10], -(1-yData(i*10-9:i*10)), 'FaceColor', col(1,:))
    hold on
    er = errorbar(-(1-yData(i*10-9:i*10)), [1:10], yErr(i*10-9:i*10), yErr(i*10-9:i*10), 'horizontal', 'lineStyle', 'none')
    er.Color = 'b';
    barh([1:10], 1-xData(i*10-9:i*10), 'FaceColor', col(2,:))
    er = errorbar(1-xData(i*10-9:i*10), [1:10], xErr(i*10-9:i*10), xErr(i*10-9:i*10), 'horizontal', 'lineStyle', 'none')
    er.Color = 'b';
    xlim([-0.45 0.45])
    set(gca, 'YDir','reverse', 'XTick', [-0.4:0.4:0.4], 'XTickLabel', [0.4,0,0.4], 'YTick',[], 'fontSize',10)
    ylim([0.5 10.5])
    Lab = [summaryTable.p1(i*10-9:i*10),summaryTable.p2(i*10-9:i*10)];
    text(repmat(-0.4, 10,1), [1:10], Lab([1:10]'+(maxCorrIdx(i*10-9:i*10)-1)*10), 'HorizontalAlignment', 'left','VerticalAlignment','middle', 'fontSize',9, 'color',[0.6 0.6 0.6])
    text(repmat(0.4, 10,1), [1:10], Lab([1:10]'+(3-maxCorrIdx(i*10-9:i*10)-1)*10), 'HorizontalAlignment', 'right','VerticalAlignment','middle', 'fontSize',9, 'color',[0.6 0.6 0.6])
    title(familyNames{i}, 'fontSize',12)
        if i==2
        xlabel('Effect of paralog deletion (1-correlation WT-mutant)', 'fontSize',12)
    end
end
% save_gf(gcf,'Fig3Bars','paper','tamar21','type','svg')
