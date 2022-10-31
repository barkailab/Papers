%% Figure 9
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
h3cc=load('./H3CC_henikoff.mat');
h3smooth=conv(mean(h3cc.p_all,2),normpdf([-50:50]',0,25),'same');
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
load('SC_genome.mat')
SC_genomeCat = upper(cat(2,SC_genome.Sequence));
load('./paraSeqs.mat');
if ~exist('checWTdelLactisSwap','var')
    load('checWTdelLactisSwap.mat')
end
GP=load('./group_imp.mat')


%% Figure 9A-D
examplePairs = {'Gis1','Rph1'};
targets = {'GSY2','HAP5'};
TargetsSortForFig5 = readtable('./TargetsSortForFig8.xlsx');
allSamples = fieldnames(checWTdelLactisSwap.sumProm);

patternsColl{1} = {'[ACG]A[AG]GGG[AT]|[AT]CCC[CT]T[CGT]','TA[GA]GGG[AT]|[AT]CCC[CT]TA'};
compType = {'_d','_DBD'};
promoterLengthPlot = 700;
promoterLengthPlotVec = repmat(promoterLengthPlot,6701,1);
quantileTH = 0.98; %for fit log2 fold change
zscoreTH = 2; %for fit log2 fold change
patternPos=[-272,-134]
xspacer = 0.03;
yspacer = 0.07;
xPos = 0.1;
yPos = 0.5;
Hsignal = 0.12;
Wsignal =  0.13;
Wscatter = 0.15;
Hscatter = Wscatter*1.6;
smallYspacer = 0.01;
Wline=Wscatter/8.5;
Hline = Hscatter;
zscoreTHclus = 4;

%colSignal = repmat([0 0 0; 0.5 0.5 0.5], 2,1);
colSignal = repmat([0 0 0],7,1);
colMotifs = lines(6);
caxisLim = [0 0.2];
cMapLac = cbrewer2('OrRd'); %brewermap(256,'OrRd');
cMapZscores = flipud(bone);
cMapMotifs = brighten(cbrewer2('Oranges',8),-0.2); %brighten(brewermap(8,'Oranges'),-0.2);
cMapMotifs = cMapMotifs([1:4],:);
cMapMotifs(1,:) = [1,1,1];
cMapClusters = [0.7 0.7 0.7; 0.2 0.2 0.2];

% cMapMotifs = brighten(cbrewer2('Purples',4),-0.2);
% cMapMotifs = brighten(cbrewer2('Oranges',8),-0.2);
% cMapMotifs = cMapMotifs([1:4],:);
% cMapMotifs(1,:) = [1,1,1];
cMapLogChange = flipud(cbrewer2('RdBu')); %flipud(brewermap(256,'RdBu'));
% cMapClusters = [[141 211 199];[61 122 143];[190 186 218]]/255;
% cMapClusters = brighten(cMapClusters,-0.4);
cMapMotifs2 = brighten(cbrewer2('PuRd'),-0.4);
cMapMotifs2(1,:) = [1,1,1];
intoGene = 20;
N = 30;
lineCol = [cbrewer2('Blues',1); cMapLac(150,:)] %[brewermap(1,'Blues'); cMapLac(150,:)];
kLacNtargets = 100;

for tf = 1:size(examplePairs,1)
    figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
    targetList = targets(tf,:);
    TF1 = examplePairs{tf,1};
    TF2 = examplePairs{tf,2};
    m1 = [TF1, '_d',upper(TF2)];
    m2 = [TF2, '_d',upper(TF1)];
    s1 = [TF1,'_',TF2,'_DBD'];
    s2 = [TF2,'_',TF1,'_DBD'];
    lac = 'Rph1_lactis';
    intStrains = {TF1,m1,s1,TF2,m2,s2,lac};
    intStrains = intStrains(ismember(intStrains,allSamples));
    pattern = patternsColl{tf};
    % scatters
    clear S
    for p = 1:numel(pattern)
        S{p} = regexp(SC_genomeCat, pattern{p});
    end
    occPattern = [];
    for i = 1:6701
        if ~isnan(promoterLengthsORF(i))
            intBases = GP.chr_idx(GP.gene_infoR64.position(i,1)) + GP.gene_infoR64.position(i,2)+[-promoterLengthsORF(i):0]*(GP.gene_infoR64.dir(i));
            for p = 1:numel(pattern)
                occPattern(i,p) = sum(ismember(S{p}, intBases));
            end
        end
    end
    
    xdataOriginal = checWTdelLactisSwap.sumProm.(TF1);
    ydataOriginal = checWTdelLactisSwap.sumProm.(TF2);
    lacIdx = find(contains(allSamples, 'lactis') & contains(allSamples, {TF1,TF2}));
    
    if isempty(lacIdx)
        colDataZ= repmat(1, 6701,1);
        colDataS= repmat(1, 6701,1);
    else
        lac = allSamples{lacIdx};
        colDataZ= nanZscore(checWTdelLactisSwap.sumProm.(lac));
        colDataS= checWTdelLactisSwap.sumProm.(lac);
    end
    
    % scatter WTs
    axes('Position', [xPos+0.1*Wscatter yPos+0.1*Hscatter 0.9*Wscatter 0.9*Hscatter])
    scatter(xdataOriginal/max(xdataOriginal), ydataOriginal/max(ydataOriginal), 20, [0.7 0.7 0.7], '.','MarkerEdgeColor',[.3 .3 .3], 'DisplayName','all genes')
    hold on
    [~,maxIdx] = maxk(colDataS,150);
    scatter(xdataOriginal(maxIdx)/max(xdataOriginal), ydataOriginal(maxIdx)/max(ydataOriginal), 20, colDataS(maxIdx)/max(colDataS),'filled','MarkerEdgeColor',[.3 .3 .3],'DisplayName','{\itK.lactis} top targets')
    colormap(gca, cMapLac)
    cbr = colorbar()
    %title(cbr, '\itK.lac', 'fontSize',12)
    title(cbr, 'K.lac', 'fontSize',12)
    axisPos = gca;
    set(gca, 'YTick',[],'XTick',[],'fontSize',10)
    title(sprintf('r=%.2f', corr(xdataOriginal, ydataOriginal, 'rows','pairwise')), 'fontSize',11,'FontWeight','normal')
    caxis([0.4,1])
    set(cbr, 'Ticks', [0.4:0.2:1])

    % histogram yaxis
    axes('Position', [xPos yPos+0.1*Hscatter 0.1*Wscatter 0.9*Hscatter])
    bins = [0:0.05:1];
    binM = movmean(bins,2,'Endpoints','discard');
    countsH = histcounts(ydataOriginal(maxIdx)/max(ydataOriginal), bins,'Normalization','pdf');
    barh(movmean(bins, 2,'Endpoints','discard'),movmean(countsH,1), 'FaceColor', cMapLac(180,:),'LineStyle','none')
    hold on
    countsAll = histcounts(ydataOriginal/max(ydataOriginal), bins,'Normalization','pdf');
    fill([countsAll([1,1:end]),0,0], [0,binM([1:end, end]),0], [.7 .7 .7], 'LineStyle','none', 'FaceAlpha',0.75)
    set(gca, 'XDir','reverse','YLim', axisPos.YLim,'XTick',[], 'fontSize',10)
    xlim([0 max(countsH)])
    set(gca,'YTick', [0:0.5:1])
    ylabel(TF2, 'fontSize',12)

    % histogram xaxis
    axes('Position', [xPos+0.1*Wscatter yPos axisPos.Position(3) 0.1*Hscatter])
    countsH = histcounts(xdataOriginal(maxIdx)/max(xdataOriginal), bins,'Normalization','pdf');
    bar(movmean(bins,2,'Endpoints','discard'), movmean(countsH,1), 'FaceColor', cMapLac(180,:),'LineStyle','none','DisplayName','K.lactis top targets')
    hold on
    countsAll = histcounts(xdataOriginal/max(xdataOriginal), bins,'Normalization','pdf');
    fill([0,binM([1:end, end]),0], [countsAll([1,1:end]),0,0], [.7 .7 .7], 'LineStyle','none', 'FaceAlpha',0.75,'DisplayName','all genes')
    set(gca, 'YDir','reverse','XLim', axisPos.XLim, 'YTick',[],'fontSize',10)
    ylim([0 max(countsH)])
    set(gca,'XTick', [0:0.5:1])
    xlabel(TF1,'fontSize',12)
    legend()
    
    % scatters with movement
    k=10;
    [~, xTopIdx] = maxk(xdataOriginal,k);
    [~, yTopIdx] = maxk(ydataOriginal,k);
    topBoth = unique([xTopIdx; yTopIdx]);
    for ct = 1:numel(compType)
        comp1Idx = find(contains(allSamples, compType{ct}) & startsWith(allSamples,TF1));
        comp2Idx = find(contains(allSamples,compType{ct}) & startsWith(allSamples,TF2));
        if ~isempty(comp1Idx)
            comp1 = allSamples{comp1Idx};
        else
            comp1 = TF1;
        end
        if ~isempty(comp2Idx)
            comp2 = allSamples{comp2Idx};
        else
            comp2 = TF2;
        end
        if ~isempty(comp1Idx) | ~isempty(comp2Idx)
            
            ydataNew = checWTdelLactisSwap.sumProm.(comp2)/max(checWTdelLactisSwap.sumProm.(comp2));
            xdataNew = checWTdelLactisSwap.sumProm.(comp1)/max(checWTdelLactisSwap.sumProm.(comp1));
            
            normF1 =mean(xdataNew(xTopIdx))/mean(xdataOriginal(xTopIdx)/max(xdataOriginal));
            normF2 = mean(ydataNew(yTopIdx))/mean(ydataOriginal(yTopIdx)/max(ydataOriginal));
            
            fitGenes1=getFitGenes(xdataOriginal,3.5,0.99);
            normF1=robustfit(xdataOriginal(fitGenes1)/max(xdataOriginal),xdataNew(fitGenes1),[],[],'off');
            
            fitGenes2=getFitGenes(ydataOriginal,3.5,0.99);
            normF2=robustfit(ydataOriginal(fitGenes2)/max(ydataOriginal),ydataNew(fitGenes2),[],[],'off');
            
            titleStr = sprintf('%s vs %s (%.2f, %.2f)', TF1, TF2, normF1,normF2);
            
            xOriginalAdjusted = xdataOriginal*normF1/max(xdataOriginal);
            yOriginalAdjusted = ydataOriginal*normF2/max(ydataOriginal);
            
            axes('position', [xPos+(Wscatter+xspacer)*ct yPos Wscatter Hscatter])
            scatter(xdataNew, ydataNew, [], [0.4 0.4 0.4], '.', 'MarkerEdgeAlpha', 0.5)
            hold on
            for g = xTopIdx'
                if ydataNew(g)>yOriginalAdjusted(g)
                    pUp = plot(xdataNew(g)*[1,1], [ydataNew(g), yOriginalAdjusted(g)], 'color', lineCol(1+(ydataNew(g)>yOriginalAdjusted(g)),:), 'lineWidth',2, 'DisplayName','mutant>WT')
                else
                    pDown = plot(xdataNew(g)*[1,1]/max(xdataNew), [ydataNew(g)/max(ydataNew), yOriginalAdjusted(g)/max(yOriginalAdjusted)], 'color', lineCol(1+(ydataNew(g)>yOriginalAdjusted(g)),:), 'lineWidth',2,'DisplayName','mutant<WT')
                end
            end
            
            for g = yTopIdx'
                plot( [xdataNew(g), xOriginalAdjusted(g)], ydataNew(g)*[1,1], 'color', [0.7 0.7 0.7], 'lineWidth',1)
            end
            scatter(xdataNew(xTopIdx), ydataNew(xTopIdx), [],occPattern(xTopIdx,3-ct), 'filled', 'MarkerEdgeColor','k')
            scatter(xdataNew(yTopIdx), ydataNew(yTopIdx), [],occPattern(yTopIdx,3-ct), 'filled', 'MarkerEdgeColor','k')
            colormap(gca,cMapMotifs)
            cbr = colorbar()
            caxis([-0.5 3.5])
                title(sprintf('r=%.2f', corr(xdataNew, ydataNew, 'rows','pairwise')), 'fontSize',11,'FontWeight','normal')
            set(cbr, 'Ticks', [0:3], 'TickLabels', {'0', '1','2','>2'})
            %             title(cbr, ['#motif ',char(64+ct)], 'fontSize',12)
            title(cbr, ['#motif',char(64+3-ct)], 'fontSize',12)
            labelSamples = regexprep({comp1,comp2}, {'_d', '_lactis', '_DBD','_'}, {' \\Delta', ' \\itK.lac', ' DBD',' '});
            set(gca,'fontSize',10)
%             xlabel(labelSamples{1},'fontSize',12)
%             ylabel(labelSamples{2},'fontSize',12)
            xlabel(strrep(comp1,'_',' '),'fontSize',12)
            ylabel(strrep(comp2,'_',' '),'fontSize',12)
            set(gca, 'XTick', [0:0.5:1], 'YTick', [0:0.5:1])
            legend(pUp,'Location','northeast')
        end
    end    
    
    % unqiue targets
    intStrains = {TF1,TF2,lac};
    clear patternCorrMat
    k=7;
    uTargetIdx = [];
    while numel(uTargetIdx)<N
        [~,idxMax{1}] = maxk(checWTdelLactisSwap.sumProm.(TF1),k);
        [~,idxMax{2}] = maxk(checWTdelLactisSwap.sumProm.(TF2),k);
        uTargetIdx = unique(cat(1, idxMax{:})', 'stable');
        k=k+1;
    end
    
    % pattern:
    sumSignalPattern=nan(numel(uTargetIdx),numel(pattern),2);
    for p = 1:numel(pattern)
        for g = 1:numel(uTargetIdx)
            currPromoterSeq = findPromoterSeq(uTargetIdx(g), GP, promoterLengthsORF, 'tss', 'ORF');
            currPromoterSeq = upper(currPromoterSeq);
            [currMatches, S, E] = regexp(currPromoterSeq, pattern{p},'match', 'start','end');
            currPromoterSeq(S) = 'N';
            [currMatches2, S2, E2] = regexp(currPromoterSeq, pattern{p},'match', 'start','end');
            currMatches = [currMatches, currMatches2];
            S = [S,S2];
            E = [E,E2];
            nPatternInPromoter(p,g)  = numel(currMatches);
            M = round((S+E)/2);
            if numel(M) == 0
                M=nan;
                E=nan;
                S=nan;
            end
            promPos = GP.gene_infoR64.position(uTargetIdx(g),2)+[-promoterLengthsORF(uTargetIdx(g)):-1]*GP.gene_infoR64.dir(uTargetIdx(g));
            intPromPos = promPos(min(abs([1:promoterLengthsORF(uTargetIdx(g))]-M'),[],1)<(25+max((E-S)/2)));
            sumSignalPattern(g,p,1) = sum(checWTdelLactisSwap.norm.(TF1){GP.gene_infoR64.position(uTargetIdx(g),1)}(intPromPos))./(sum(checWTdelLactisSwap.norm.(TF1){GP.gene_infoR64.position(uTargetIdx(g),1)}(promPos)));
            sumSignalPattern(g,p,2) = sum(checWTdelLactisSwap.norm.(TF2){GP.gene_infoR64.position(uTargetIdx(g),1)}(intPromPos))./(sum(checWTdelLactisSwap.norm.(TF2){GP.gene_infoR64.position(uTargetIdx(g),1)}(promPos)));
            forProm(g).pos{p}=GP.chrIdx(GP.gene_infoR64.position(uTargetIdx(g),1))+promPos;
            forProm(g).intPos{p}=GP.chrIdx(GP.gene_infoR64.position(uTargetIdx(g),1))+intPromPos;
        end
    end
    nPatternInPromoter = occPattern(uTargetIdx,:)';
    
    % zscore:
    intStrains = intStrains(~contains(intStrains,'_d'));
    kLacIdx=find(contains(allSamples,intStrains)&endsWith(allSamples,'lactis'));
    if numel(kLacIdx)
        [~, idxVec, sumPromRep, repIdx] = getRepeatsCorr([intStrains allSamples{kLacIdx}],'dataType','sumProm')
        [~,uRpt]=unique(repIdx,'stable')
        sumPromRep=sumPromRep(:,uRpt)
        idxVec=idxVec(uRpt)
        nRpt=accumarray(idxVec,1)
        maxSize=lcm(lcm(sum(idxVec==1),sum(idxVec==2)),2*sum(idxVec==3));    
    else
        [~, idxVec, sumPromRep, repIdx] = getRepeatsCorr(intStrains,'dataType','sumProm')
        [~,uRpt]=unique(repIdx,'stable')
        sumPromRep=sumPromRep(:,uRpt)
        idxVec=idxVec(uRpt)
        nRpt=accumarray(idxVec,1)
        maxSize=lcm(sum(idxVec==1),sum(idxVec==2));            
    end
    zScoresMat=cell(1,sum(nRpt));
    for j=1:sum(nRpt)
        zScoresMat{j} = repmat(nanZscore(sumPromRep(:,j)),1,maxSize/nRpt(idxVec(j))/(1+(idxVec(j)==3)));
    end
    zScoresMat=cell2mat(zScoresMat);
    zScoreSelectedTargets =  zScoresMat(uTargetIdx,:);
    clearvars zScoresMat
    zScoresMat = nan(6701,numel(intStrains));
    for j = 1: numel(intStrains)
        zScoresMat(:,j) = nanZscore(checWTdelLactisSwap.sumProm.(intStrains{j}));
    end 
    zScoreSelectedTargetsMean =  zScoresMat(uTargetIdx,:);
    clearvars zScoresMat
    
    % log2 fold change
    clear normFactor logChange
    for par = 1:2
        currF1 = eval(sprintf('TF%d',par));
        currF2 = eval(sprintf('TF%d',3-par));
        zscore_WT1 = nanZscore(checWTdelLactisSwap.sumProm.(currF1));
        fitGenes =  zscore_WT1 > min(zscoreTH, quantile(zscore_WT1, quantileTH));
        for ct = 1:numel(compType)
            m = allSamples{startsWith(allSamples, currF1) & contains(allSamples, compType{ct})};
            if numel(m) > 0
                bestFit = robustfit(checWTdelLactisSwap.sumProm.(currF1)(fitGenes), checWTdelLactisSwap.sumProm.(m)(fitGenes), [], [],'off' );
                normFactor(par) = log2(bestFit);
                logChange{par}(ct,:) = log2(checWTdelLactisSwap.sumProm.(m)+700)- log2(checWTdelLactisSwap.sumProm.(currF1)+700)- normFactor(par);
                %                 significantGenes(par,:) =  abs((logChange(par,:)-median(logChange(par,fitGenes)))/std(logChange(par,fitGenes))) >= 1;
            else
                logChange{par}(ct,:) = nan(1,6701);
            end
        end
    end
    logChangeSelectedTargets =  logChange{2}(:,uTargetIdx);
    
    % cluster by targets
    [~,idx] = ismember(TargetsSortForFig5.([examplePairs{tf,1},'_',examplePairs{tf,2}]), uTargetIdx);
    idx(idx==0) = [];
   
    % imagesc Zscores
    axes('Position', [xPos yPos-1.5*yspacer-Hline Wline*3 Hline])
    imagesc(zScoreSelectedTargets(idx,:))
    hold on
    plot([maxSize 2*maxSize]+[0.5;.5],ylim,'k-')
    plot(repmat(xlim',1,numel(idx)-1),repmat(1.5:numel(idx),2,1),'k-')
    caxis([0 10])
    set(gca,'YTick', [],  'XTick',[(maxSize+1)/2 maxSize+(maxSize+1)/2 maxSize*2+(maxSize/2+1)/2], 'XTickLabel', {TF1,TF2,'klac'}, 'fontSize',11, 'TickLength',[0 0]);
    colormap(gca,cMapZscores)
    if tf==1
        text(1.5,-3,'Binding', 'fontSize',12, 'HorizontalAlignment', 'center')
        text(1.5,-1, '(wt)', 'fontSize',12, 'HorizontalAlignment', 'center')
    end
    xlim([.5 size(zScoreSelectedTargets,2)+.5])
    
    colormap(gca,cMapZscores)
    [~, targetIdx]=ismember(targetList,GP.gene_infoR64.nameNew)
    [~, targetPos]=ismember(targetIdx, uTargetIdx(idx))
    text(2,-3,'Binding', 'fontSize',11, 'HorizontalAlignment', 'center')
    text(2,-1, '(wt)', 'fontSize',11, 'HorizontalAlignment', 'center')
  
    axes('Position', [xPos-Wline/2 yPos-1.5*yspacer-Hline Wline/2 Hline])
    scatter(repmat(0,2,1),targetPos,15, cMapClusters, 'filled');
    ylim([0 numel(idx)]+.5);
    set(gca,'YDir','reverse')
    axis off
    
    % imagesc log2 fold change
    axes('Position', [xPos+Wline*3.5 yPos-1.5*yspacer-Hline Wline*2 Hline])
    imagesc(logChangeSelectedTargets(:,idx)')
    plotgrid(logChangeSelectedTargets(:,idx)')
    caxis([-2 2])
%     set(gca, 'XTick', [1,2], 'XTickLabel', {'\DeltaP', 'swap'},'YTick',[], 'TickLength',[0 0], 'fontSize',11,  'XTickLabelRotation', 30)
    set(gca, 'XTick', [1,2], 'XTickLabel', {'P', 'swap'},'YTick',[], 'TickLength',[0 0], 'fontSize',11,  'XTickLabelRotation', 30)

    colormap(gca,cMapLogChange)
    text(mean(xlim()),-3,'Rph1 binding', 'fontSize',11, 'HorizontalAlignment', 'center')
    text(mean(xlim()),-1,'change', 'fontSize',11, 'HorizontalAlignment', 'center')

    % motifs 
    temp=getMotifScoresRepeats(forProm(idx),intStrains);
    axes('Position', [xPos+Wline*6 yPos-1.5*yspacer-Hline Wline*2 Hline])

    [~,idxForPattern] = max( zScoreSelectedTargetsMean(idx,:),[],2);
    idxForPattern=sub2ind(size(temp),[1:size(temp,1)]',idxForPattern)
    %idxForPattern=idxForPattern';
    sumSignalPatternPlot = cell2mat(cat(1,temp{idxForPattern}))%[sumSignalPattern(sub2ind(size(sumSignalPattern), 1:size(sumSignalPattern,1), repmat(1,1, size(sumSignalPattern,1)), idxForPattern)); sumSignalPattern(sub2ind(size(sumSignalPattern), 1:size(sumSignalPattern,1), repmat(2,1, size(sumSignalPattern,1)), idxForPattern))]';
    imagesc(sumSignalPatternPlot,[0 .8])
    hold on
    plot(repmat(xlim',1,numel(idx)-1),repmat(1.5:numel(idx),2,1),'k-')
    plot(numel(temp{1}{1})+[0.5;0.5],ylim,'k-')
    set(gca, 'XTick', [1:numel(pattern)], 'XTickLabel', {'mA','mB'},'YTick',[], 'TickLength',[0 0 ],'fontSize',11,  'XTickLabelRotation', 30)
    colormap(gca,cMapMotifs2)
    text(mean(xlim),-3,'Motif', 'fontSize',11, 'HorizontalAlignment', 'center')
    text(mean(xlim),-1,'score', 'fontSize',11, 'HorizontalAlignment', 'center')
    
    intStrains = {TF1,m1,s1,TF2,m2,s2,lac};
    intStrains = intStrains(ismember(intStrains,allSamples));
    
    % signal on example promoters
    for g = 1: length(targetList)
        if ~isempty(targetList{g})
            target = targetList{g};
            geneIdx = GP.gene_table.(target);
            pos = GP.gene_infoR64.position(geneIdx,:);
            intBasesGene = pos(2)+[-promoterLengthPlotVec(geneIdx):intoGene].*GP.gene_infoR64.dir(geneIdx);
            %intBases = pos(2)+[-promoterLengthPlot:100].*GP.gene_infoR64.dir(geneIdx);
            
            clear signalMat
            for p = 1:length(intStrains)
                normProfile = chromosomes2fullProfile(checWTdelLactisSwap, intStrains(p));
                maxNorm = max(normProfile(logical(promoterIDXvec)));
                signalMat(p,:) = checWTdelLactisSwap.norm.(intStrains{p}){pos(1)}(intBasesGene)./maxNorm;
            end
            patternIdx = []; patternIdxRC = [];motifPatternLocSAll = []; motifPatternLocEAll = [];motifPatternLocSRCAll=[]; motifPatternLocERCAll=[];
            currPromSeq = upper(findPromoterSeq(geneIdx, GP, promoterLengthPlotVec, 'tss', 'ORF'));
            for p = 1:length(pattern)
                motifPattern = pattern{p};
                [motifPatternLocS,motifPatternLocE] = regexp(currPromSeq, motifPattern, 'start','end');
                motifPatternLocS = motifPatternLocS-promoterLengthPlotVec(geneIdx);
                motifPatternLocSAbs{p} = motifPatternLocS*GP.gene_infoR64.dir(geneIdx)+pos(2)+GP.chr_idx(pos(1));
                motifPatternLocEAbs{p} = motifPatternLocE*GP.gene_infoR64.dir(geneIdx)+pos(2)+GP.chr_idx(pos(1));
                motifPatternLocE = motifPatternLocE-promoterLengthPlotVec(geneIdx);
                patternIdx = [patternIdx; repmat(p, length(motifPatternLocS),1)];
                motifPatternLocSAll = [motifPatternLocSAll,motifPatternLocS];
                motifPatternLocEAll = [motifPatternLocEAll,motifPatternLocE];
            end
            
            labelSamples = regexprep(intStrains, {'_d', '_lactis', '_DBD','_'}, {' \\Delta', ' \\itK.lac', ' DBD',' '});
            if g == 1
                axes('Position', [xPos+Wsignal+3.1*xspacer+(g-1)*(Wsignal+xspacer) yPos-1.5*yspacer-Hsignal Wsignal Hsignal])
            else
                axes('Position', [xPos+Wsignal+2.45*xspacer+(g-1)*(Wsignal+xspacer) yPos-1.5*yspacer-Hsignal Wsignal Hsignal])
            end
            hold on
            for z = 1:numel(intStrains)
                plot([-promoterLengthPlotVec(geneIdx):intoGene], rescale(signalMat(z,:),numel(intStrains)-z+1.5, numel(intStrains)-z+2.3,'InputMin',0, 'InputMax',0.4),...
                    'color', colSignal(z,:))
            end
            plot(-promoterLengthPlotVec(geneIdx):intoGene, rescale(h3smooth(GP.chr_idx(pos(1))+intBasesGene), 0.5, 1.3),'k', 'LineWidth',1)
            area(-promoterLengthPlotVec(geneIdx):intoGene, rescale(h3smooth(GP.chr_idx(pos(1))+intBasesGene), 0.5, 1.3),'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'LineStyle', 'none')
            ylim([0.5 numel(intStrains)+1.5])
            xlim([-promoterLengthPlotVec(geneIdx) intoGene])

            TssLoc = (GP.gene_infoR64.tamarTss(geneIdx,2)-pos(2))*GP.gene_infoR64.dir(geneIdx);
            if g==1
                %set(gca, 'YTick', [1:numel(labelSamples)], 'YTickLabel', labelSamples(end:-1:1), 'fontSize',12, 'XTick', [-800:200:-400, TssLoc,100], 'XTickLabel', {'-800','-600','-400','TSS', '100'})
                %set(gca, 'YTick', [1:numel(labelSamples)+1], 'YTickLabel', {'NucOcc',labelSamples{end:-1:1}},'XTick',[], 'fontSize',11, 'TickLength',[0 0])
                set(gca, 'YTick', [1:numel(labelSamples)+1], 'YTickLabel', ['NucOcc',strrep(labelSamples(end:-1:1), '\','')],'XTick',[], 'fontSize',11, 'TickLength',[0 0])

            else
                set(gca, 'YTick', [],'XTick',[], 'fontSize',11)
            end
            plot(TssLoc*[1 1] ,ylim, ':r', 'LineWidth',1.5)
            scatter(0 ,0.5,[], '>r','filled')
            if ~isempty(motifPatternLocSAll)
                for p = 1:length(patternIdx)
                    %scatter(motifPatternLocSAll(p), 0.5,50, [colMotifs(patternIdx(p),:)],'filled', 'MarkerEdgeColor','k');
                    fill([motifPatternLocSAll(p)*[1,1], motifPatternLocEAll(p)*[1,1]], [0.5, numel(intStrains)+1.5.*[1,1], 0.5], colMotifs(patternIdx(p),:),...
                        'LineStyle','none', 'FaceAlpha', 0.3)
                end
            end
            if g==1
                title(sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s',cMapClusters(1,:),target),'fontSize',12)
            else
                title(sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s',cMapClusters(2,:),target),'fontSize',12)
            end
            box on
            
            if g == 1
                axes('Position', [xPos+Wsignal+3.1*xspacer+(g-1)*(Wsignal+xspacer) yPos-1.5*yspacer-Hline Wsignal Hline-Hsignal-smallYspacer])
            else
                axes('Position', [xPos+Wsignal+2.45*xspacer+(g-1)*(Wsignal+xspacer) yPos-1.5*yspacer-Hline Wsignal Hline-Hsignal-smallYspacer])
            end
            hold on
            for z = 1:numel(intStrains)
                plot([-promoterLengthPlotVec(geneIdx):intoGene], rescale(signalMat(z,:),numel(intStrains)-z+1.5, numel(intStrains)-z+2.3,'InputMin',0, 'InputMax',0.4),...
                    'color', colSignal(z,:))
            end
            plot(-promoterLengthPlotVec(geneIdx):intoGene, rescale(h3smooth(GP.chr_idx(pos(1))+intBasesGene), 0.5, 1.3),'k', 'LineWidth',1)
            area(-promoterLengthPlotVec(geneIdx):intoGene, rescale(h3smooth(GP.chr_idx(pos(1))+intBasesGene), 0.5, 1.3),'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'LineStyle', 'none')
            ylim([0.5 numel(intStrains)+1.5])
            xlim([-promoterLengthPlotVec(geneIdx) intoGene])

            ylim([1.5 numel(intStrains)+1.5])
            xlim([-promoterLengthPlotVec(geneIdx) intoGene])
            TssLoc = (GP.gene_infoR64.tamarTss(geneIdx,2)-pos(2))*GP.gene_infoR64.dir(geneIdx);
            xlim(patternPos(g)+[-50 50])
            if g==1
                %set(gca, 'YTick', [1:numel(labelSamples)], 'YTickLabel', labelSamples(end:-1:1), 'fontSize',12, 'XTick', [-800:200:-400, TssLoc,100], 'XTickLabel', {'-800','-600','-400','TSS', '100'})
%                 set(gca, 'YTick', [1:numel(labelSamples)+1], 'YTickLabel', {'NucOcc',labelSamples{end:-1:1}},'XTick',patternPos(g)+[-45 45], 'fontSize',11, 'TickLength',[0 0],...
%                     'XTickLabel',{})
                set(gca, 'YTick', [1:numel(labelSamples)+1], 'YTickLabel', ['NucOcc',strrep(labelSamples(end:-1:1),'\', '')],'XTick',patternPos(g)+[-45 45], 'fontSize',11, 'TickLength',[0 0],...
                    'XTickLabel',{})
            else
                set(gca, 'YTick', [],'XTick',patternPos(g)+[-45 45], 'fontSize',11,'XTickLabel',{})
            end
            text(patternPos(g)+[-45 45],repmat(1.5,2,1),sprintfc('%d',patternPos(g)+[-45 45]),'VerticalAlignment','top','HorizontalAlignment','center')
            plot(TssLoc*[1 1] ,ylim, ':r', 'LineWidth',1.5)
            scatter(0 ,0.5,[], '>r','filled')
            if ~isempty(motifPatternLocSAll)
                for p = 1:length(patternIdx)
                    %scatter(motifPatternLocSAll(p), 0.5,50, [colMotifs(patternIdx(p),:)],'filled', 'MarkerEdgeColor','k');
                    fill([motifPatternLocSAll(p)*[1,1], motifPatternLocEAll(p)*[1,1]], [0.5, numel(intStrains)+1.5.*[1,1], 0.5], colMotifs(patternIdx(p),:),...
                        'LineStyle','none', 'FaceAlpha', 0.3)
                end
            end
            box on
        end
    end
end

% colorbars:
axes('Position', [xPos yPos-2.5*yspacer-Hline Wline*3 0.05])
caxis([0 10])
colormap(gca,cMapZscores)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'z-score', 'fontSize',10)
set(cbr, 'Ticks', [0.1+min(caxis),max(caxis)-1] , 'TickLabels', caxis)

axes('Position', [xPos+Wline*3.5 yPos-2.5*yspacer-Hline Wline*2 0.05])
caxis([-2 2])
colormap(gca,cMapLogChange)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'Deltalog2', 'fontSize',10)
set(cbr, 'Ticks', [0.5+min(caxis),0,max(caxis)-0.1] , 'TickLabels', [min(caxis),0,max(caxis)],'TickLength',[0 0])

axes('Position', [xPos+Wline*6 yPos-2.5*yspacer-Hline Wline*2 0.05])
caxis([0 .8])
colormap(gca,cMapMotifs2)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'motif signal', 'fontSize',10)
set(cbr, 'Ticks', [0, 0.8], 'TickLabels', {'0','80%'})

set(gcf, 'Renderer','painters')
% saveas(gcf, 'Fig6Gis1Rph1Scatters.svg')


%% Figure 9E - sequence alignment of JmjC and DBD of Gis1/Rph1
examplePairs = {'Gis1','Rph1'};
allDomain=readtable('./allDomainsAdj.txt');
chosenDomains = {'PF00096', 'PF02373', 'PF02375'}; % ZF, JmjC, JmjN
domainsLabel = {'ZF', 'JmjC','JmjN'};
domains = allDomain(contains(allDomain.queryName, 'Scer') & contains(allDomain.queryName, {'GIS1','RPH1'})& contains(allDomain.accession, chosenDomains) & allDomain.score_1>0,:);
domainsLac = allDomain(contains(allDomain.queryName, 'Klac') & contains(allDomain.queryName, {'GIS1','RPH1'})& contains(allDomain.accession, chosenDomains) & allDomain.score_1>0,:)
rIdxTable = find(contains(summaryTable.p1, examplePairs));
maxL = max(summaryTable.proteinLength(rIdxTable,:));

% schemes
domainsCol = cbrewer2('Pastel1',3);
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
axes('Position', [0.1 0.5 0.15 0.1])
for i = 1:size(examplePairs,2)
    TF = examplePairs{i};
    currDomainsIdx = find(contains(domains.queryName,upper(TF)));
    currL = domains.qlen(currDomainsIdx(1));
    plot([1 currL], [i i], 'color', 'k','LineWidth',2)
    hold on
    for d = currDomainsIdx'
        domainType = find(contains(chosenDomains, extractBefore(domains.accession(d),'.')));
        rectangle('Position', [domains.from_2(d), i-0.2, diff([domains.from_2(d), domains.to_2(d)]), 0.4],...
            'Curvature',0.2, 'FaceColor',domainsCol(domainType,:), 'EdgeColor', brighten(domainsCol(domainType,:),-0.8))
        if i==1
            text(mean([domains.from_2(d), domains.to_2(d)]), i-0.3 ,domainsLabel{domainType},...
                'fontSize',10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','Color',brighten(domainsCol(domainType,:),-0.6))
        end
    end
    text(-5,i,TF,'fontSize',12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
    text(currL+5, i ,num2str(currL),'fontSize',10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
end

currDomainsIdx = [1:4];
currL = domainsLac.qlen(1);
plot([1 currL], [3 3], 'color', 'k','LineWidth',2)
hold on
for d = currDomainsIdx
    domainType = find(contains(chosenDomains, extractBefore(domainsLac.accession(d),'.')));
    rectangle('Position', [domainsLac.from_2(d), 3-0.2, diff([domainsLac.from_2(d), domainsLac.to_2(d)]), 0.4],...
        'Curvature',0.2, 'FaceColor',domainsCol(domainType,:), 'EdgeColor', brighten(domainsCol(domainType,:),-0.8))
end
text(-5,3,'K.lac','fontSize',12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
text(currL+5, 3 ,num2str(currL),'fontSize',10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
xlim([1 maxL])
ylim([0.5 3.5])
set(gca, 'YDir', 'reverse')
axis off

% jmjc mutation
anno = [domains.from_2(contains(domains.accession, chosenDomains{2})), domains.to_2(contains(domains.accession, chosenDomains{2}))]
anno(3,:) = [domainsLac.from_2(contains(domainsLac.accession, chosenDomains{2})), domainsLac.to_2(contains(domainsLac.accession, chosenDomains{2}))];
mutationPos = -16;
paraIdx = find(contains({paraSeqs.name}, upper(examplePairs)));
if strcmp(paraSeqs(paraIdx).name, upper(examplePairs{1}))
    seq{2} = paraSeqs(paraIdx).seq{22};
    seq{1} = paraSeqs(paraIdx).seq{12};
else
    seq{1} = paraSeqs(paraIdx).seq{22};
    seq{2} = paraSeqs(paraIdx).seq{12};
end
seq{3} = paraSeqs(paraIdx).seq{18};
%[~,seqAlignPar] = nwalign(seq{1}(anno(1,1):anno(1,2)), seq{2}(anno(2,1):anno(2,2)));
for i =1:3
    seqsForAlign{i} = seq{i}(anno(i,1):anno(i,2));
end
seqAlign = multialign(seqsForAlign);
seqAlignPlot = seqAlign(:,end+mutationPos+[-9:9]);
alignStart = anno(:,2)+mutationPos-9;
alignEnd = anno(:,2)+mutationPos+9;

text(-300, 3, {'Gis1','Rph1','K.lac'}, 'Color', 'k', 'FontWeight', 'normal', 'FontName','Courier', 'HorizontalAlignment','left', 'VerticalAlignment','top')
text(-100, 3, strcat(num2str(alignStart([2,1,3])),repmat(char(8230),3,1)), 'Color', [.5 .5 .5], 'FontWeight', 'normal', 'FontName','Courier', 'HorizontalAlignment','center', 'VerticalAlignment','top')
for z = 1:length(seqAlignPlot)
    if z == 11
        colMut = [1 0 0];
        fWeight = 'bold';
    elseif ~strcmp(seqAlignPlot(1,z), seqAlignPlot(2,z))
        colMut = [0 0 0];
        fWeight = 'normal';
    else
        colMut = [.5 .5 .5];
        fWeight = 'normal';
    end
    text(-50+z*25, 3, seqAlignPlot([2,1,3],z), 'Color', colMut, 'FontWeight', fWeight, 'FontName','Courier', 'HorizontalAlignment','center', 'VerticalAlignment','top')
end
text((z+1)*25, 3, strcat(repmat(char(8230),3,1),num2str(alignEnd([2,1,3]))), 'Color', [.5 .5 .5], 'FontWeight', 'normal', 'FontName','Courier', 'HorizontalAlignment','center', 'VerticalAlignment', 'top')

% DBD mutation
anno = [domains.from_2(contains(domains.accession, chosenDomains{1})), domains.to_2(contains(domains.accession, chosenDomains{1}))]
anno([5,6],:) = [domainsLac.from_2(contains(domainsLac.accession, chosenDomains{1})), domainsLac.to_2(contains(domainsLac.accession, chosenDomains{1}))];
mutationPos = 9;
% [~,seqAlign] = nwalign(seq{1}(anno(1,1):anno(2,2)), seq{2}(anno(3,1):anno(4,2)));
for i =1:3
    seqsForAlign{i} = seq{i}(anno(i*2-1,1):anno(i*2,2));
end
seqAlign = multialign(seqsForAlign);
seqAlignPlot = seqAlign(:,mutationPos+[-5:5]);
alignStart = anno([1:2:end],1)+mutationPos-5;
alignEnd = anno([1:2:end],1)+mutationPos+5;

text(625, 3, strcat(num2str(alignStart([2,1,3])),repmat(char(8230),3,1)), 'Color', [.5 .5 .5], 'FontWeight', 'normal', 'FontName','Courier', 'HorizontalAlignment','center', 'VerticalAlignment','top')
for z = 1:length(seqAlignPlot)
    if z == 6
        colMut = [1 0 0];
        fWeight = 'bold';
    elseif strcmp(seqAlignPlot(2,z), ' ')|strcmp(seqAlignPlot(2,z), ':')
        colMut = [0 0 0];
        fWeight = 'normal';
    else
        colMut = [.5 .5 .5];
        fWeight = 'normal';
    end
    text(675+z*25, 3, seqAlignPlot([2,1,3],z), 'Color', colMut, 'FontWeight', fWeight, 'FontName','Courier', 'HorizontalAlignment','center','VerticalAlignment','top')
end
text(725+(z+1)*25, 3, strcat(repmat(char(8230),3,1),num2str(alignEnd([2,1,3]))), 'Color', [.5 .5 .5], 'FontWeight', 'normal', 'FontName','Courier', 'HorizontalAlignment','center','VerticalAlignment','top')
axis off

%set(gcf,'Renderer', 'painters')
% saveas(gcf, 'Fig6sequenceAlignment.svg')
