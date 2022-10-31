%% Figure 9—figure supplement 1.
clearvars -except checWTdelLactisSwap GP
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
if ~exist('GP','var')
    GP=load('./group_imp.mat')
end


%% Figure 9—figure supplement 1B
clearvars -except GP checWTdelLactisSwap SC_genomeCat allSamples
load('./promoterOL.mat','perOverlap','totalOL')
TargetsSortForFig5=readtable('./TargetsSortForFig8.xlsx');
zscoreTHclus = 4;
examplePairs = {'Vhr2','Vhr1'; 'Yap2','Yap1'; 'Swi5','Ace2'; 'Rlm1','Smp1';'Fkh2','Fkh1'};
allSamples = fieldnames(checWTdelLactisSwap.sumProm);

patternsColl{1} = {'TGACT[CA]|[TG]AGTCA'}
patternsColl{2} = {'T[TG]AC[TA]AA|TT[AT]GT[CA]A'};
patternsColl{3} = {'GCTGG|CCAGC'};
patternsColl{4} = {'TTAATAAA|TTTATTAA'};
patternsColl{5} = {'[GA]T[CA]AACA|TGTT[TG]A[CT]'}
    
compType = {'_d','_DBD'};
promoterLengthPlot = 700;
promoterLengthPlotVec = repmat(promoterLengthPlot,6701,1);
quantileTH = 0.98; %for fit log2 fold change
zscoreTH = 2; %for fit log2 fold change

Wscatter = 0.15;
Hscatter = Wscatter*1.6;
xspacer = 0.03;
yspacer = 0.07;
xPos = repmat([0.05, 0.05],3,1);
yPos = repmat([0.7,0.7-Hscatter-yspacer]',1,3);
Hsignal = 0.15;
Wsignal =  0.13;
smallYspacer = 0.01;
Hline = Hscatter% (Wscatter-3*(smallYspacer))/7.5;
Wline=Wscatter/10;
%colSignal = repmat([0 0 0; 0.5 0.5 0.5], 2,1);
colSignal = repmat([0 0 0],7,1);
colMotifs = lines(6);
caxisLim = [0 0.2];
cMapLac = brewermap(256,'OrRd');%cbrewer2('OrRd');
cMapZscores = flipud(bone);
cMapMotifs = brighten(brewermap(8,'Oranges'),-0.2);
cMapMotifs = cMapMotifs([1:4],:);
cMapMotifs(1,:) = [1,1,1];
% cMapMotifs = brighten(cbrewer2('Purples',4),-0.2);
% cMapMotifs = brighten(cbrewer2('Oranges',8),-0.2);
% cMapMotifs = cMapMotifs([1:4],:);
% cMapMotifs(1,:) = [1,1,1];
cMapLogChange = flipud(brewermap(128,'RdBu'));
cMapClusters = [[141 211 199];[61 122 143];[190 186 218]]/255;
cMapClusters = brighten(cMapClusters,-0.4);
intoGene = 20;
N = 30;

lineCol = [brewermap(1,'Blues'); cMapLac(150,:)];
kLacNtargets = 100;
figure('Units','pixels','Position',[1 41 1920 962], 'color','w')
for tf = 1:size(examplePairs,1)
    if ismember(tf,[3,5])        
        figure('Units','pixels','Position',[1 41 1920 962], 'color','w')
    end
    TF1 = examplePairs{tf,1};
    TF2 = examplePairs{tf,2};
    m1 = [TF1, '_d',upper(TF2)];
    m2 = [TF2, '_d',upper(TF1)];
    s1 = [TF1,'_',TF2,'_DBD'];
    s2 = [TF2,'_',TF1,'_DBD'];
    lac = allSamples{contains(allSamples, {TF1,TF2}) & contains(allSamples, '_lactis')};
    
    intStrains = {TF1,m1,s1,TF2,m2,s2,lac};
    intStrains = intStrains(ismember(intStrains,allSamples));
    
    % scatters
    clear S
    pattern = patternsColl{tf};
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
    xdataOriginal=xdataOriginal./max(xdataOriginal);
    ydataOriginal=ydataOriginal./max(ydataOriginal);

    lacIdx = find(contains(allSamples, 'lactis') & contains(allSamples, {TF1,TF2}));
    
    if isempty(lacIdx)
        colDataZ= repmat(1, 6701,1);
        colDataS= repmat(1, 6701,1);
    else
        lac = allSamples{lacIdx};
        colDataZ= nanZscore(checWTdelLactisSwap.sumProm.(lac));
        colDataS= checWTdelLactisSwap.sumProm.(lac);
        colDataS=colDataS./max(colDataS);
    end
    
    % scatter WTs
    axes('Position', [xPos(tf)+0.1*Wscatter yPos(tf)+0.1*Hscatter 0.9*Wscatter 0.9*Hscatter])
    scatter(xdataOriginal, ydataOriginal,20,[0.7 0.7 0.7],'.','MarkerEdgeColor',[.3 .3 .3], 'DisplayName','all genes')
    hold on
    [~,maxIdx] = maxk(colDataS,150);
    scatter(xdataOriginal(maxIdx), ydataOriginal(maxIdx),20,colDataS(maxIdx),'filled','MarkerEdgeColor',[.3 .3 .3],'DisplayName','K.lactis top targets')
    setAxisExponent()
    colormap(gca, cMapLac)
    cbr = colorbar()
    title(cbr, 'K.lac', 'fontSize',12)
    setAxisExponent('axes',cbr,'axis','c')
    axisPos = gca;
    set(gca, 'YTick',[],'XTick',[],'fontSize',10)
    
    % histogram yaxis
    axes('Position', [xPos(tf) yPos(tf)+0.1*Hscatter 0.1*Wscatter 0.9*Hscatter])
    bins = [0:0.05:1];
    binM = movmean(bins,2,'Endpoints','discard');
    countsH = histcounts(ydataOriginal(maxIdx), bins*max(ydataOriginal),'Normalization','pdf');
    barh(movmean(bins*max(ydataOriginal),2,'Endpoints','discard'),movmean(countsH,1), 'FaceColor', cMapLac(180,:),'LineStyle','none')
    hold on
    countsAll = histcounts(ydataOriginal, bins*max(ydataOriginal),'Normalization','pdf');
    fill([countsAll([1,1:end]),0,0], [0,binM([1:end, end]),0]*max(ydataOriginal),[.7 .7 .7], 'LineStyle','none', 'FaceAlpha',0.75)
    set(gca, 'XDir','reverse','YLim', axisPos.YLim,'XTick',[], 'fontSize',10)
    setAxisExponent()
    ylabel(TF2, 'fontSize',12)
    xlim([0 max(countsH)])
    
    % histogram xaxis
    axes('Position', [xPos(tf)+0.1*Wscatter yPos(tf) axisPos.Position(3) 0.1*Hscatter])
    countsH = histcounts(xdataOriginal(maxIdx), bins*max(xdataOriginal),'Normalization','pdf');
    bar(movmean(bins*max(xdataOriginal),2,'Endpoints','discard'), movmean(countsH,1), 'FaceColor', cMapLac(180,:),'LineStyle','none','DisplayName','K.lactis top targets')
    hold on
    countsAll = histcounts(xdataOriginal, bins*max(xdataOriginal),'Normalization','pdf');
    fill([0,binM([1:end, end]),0]*max(xdataOriginal), [countsAll([1,1:end]),0,0],[.7 .7 .7], 'LineStyle','none', 'FaceAlpha',0.75,'DisplayName','all genes')
    set(gca, 'YDir','reverse','XLim', axisPos.XLim, 'YTick',[],'fontSize',10)
    setAxisExponent()
    xlabel(TF1,'fontSize',12)
    ylim([0 max(countsH)])
    legend()
    
    % scatters with movement
    k=10;
    [~, xTopIdx] = maxk(xdataOriginal,k);
    [~, yTopIdx] = maxk(ydataOriginal,k);
    topBoth = unique([xTopIdx; yTopIdx]);
    if ismember(tf,[2,4])
        ctChoice =1;
    else
        ctChoice =[1,2];
    end
    for ct = ctChoice
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
            
            ydataNew = checWTdelLactisSwap.sumProm.(comp2);
            ydataNew=ydataNew./max(ydataNew);
            xdataNew = checWTdelLactisSwap.sumProm.(comp1);
            xdataNew=xdataNew./max(xdataNew);
            normF1 =mean(xdataNew(xTopIdx))/mean(xdataOriginal(xTopIdx));
            normF2 = mean(ydataNew(yTopIdx))/mean(ydataOriginal(yTopIdx));
            titleStr = sprintf('%s vs %s (%.2f, %.2f)', TF1, TF2, normF1,normF2);
            
            xOriginalAdjusted = xdataOriginal*normF1;
            yOriginalAdjusted = ydataOriginal*normF2;
            
            axes('position', [xPos(tf)+(Wscatter+xspacer)*ct yPos(tf) Wscatter Hscatter])
            scatter(xdataNew, ydataNew, [], [0.4 0.4 0.4], '.', 'MarkerEdgeAlpha', 0.5)
            hold on
            for g = xTopIdx'
                if ydataNew(g)>yOriginalAdjusted(g)
                    pUp = plot(xdataNew(g)*[1,1], [ydataNew(g), yOriginalAdjusted(g)], 'color', lineCol(1+(ydataNew(g)>yOriginalAdjusted(g)),:), 'lineWidth',2, 'DisplayName','mutant>WT')
                else
                    pDown = plot(xdataNew(g)*[1,1], [ydataNew(g), yOriginalAdjusted(g)], 'color', lineCol(1+(ydataNew(g)>yOriginalAdjusted(g)),:), 'lineWidth',2,'DisplayName','mutant<WT')
                end
            end
            
            for g = yTopIdx'
                plot( [xdataNew(g), xOriginalAdjusted(g)], ydataNew(g)*[1,1], 'color', [0.7 0.7 0.7], 'lineWidth',1)
            end
            scatter(xdataNew(xTopIdx), ydataNew(xTopIdx), [],occPattern(xTopIdx,1), 'filled', 'MarkerEdgeColor','k')
            scatter(xdataNew(yTopIdx), ydataNew(yTopIdx), [],occPattern(yTopIdx,1), 'filled', 'MarkerEdgeColor','k')
            colormap(gca,cMapMotifs)
            
%             scatter(xdataNew(xTopIdx), ydataNew(xTopIdx), 'filled', 'MarkerEdgeColor','k')
%             scatter(xdataNew(yTopIdx), ydataNew(yTopIdx), 'filled', 'MarkerEdgeColor','k')
            cbr = colorbar()
            caxis([-0.5 3.5])
            set(cbr, 'Ticks', [0:3], 'TickLabels', {'0', '1','2','>2'})
            %             title(cbr, ['#motif ',char(64+ct)], 'fontSize',12)
            title(cbr, '#motifs  ', 'fontSize',12)
            labelSamples = regexprep({comp1,comp2}, {'_d', '_lactis', '_DBD','_'}, {' Delta', ' K.lac', ' DBD',' '});
            set(gca,'fontSize',10)
            xlabel(labelSamples{1},'fontSize',12)
            ylabel(labelSamples{2},'fontSize',12)
            setAxisExponent()
            legend([pUp,pDown],'Location','northeast')          
        end
    end

    % zscores
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

    %[~,selProms]=ismember(targetList,GP.gene_infoR64.nameNew);
    %uTargetIdx=[uTargetIdx(ismember(uTargetIdx,selProms)),uTargetIdx(~ismember(uTargetIdx,selProms))];
    %uTargetIdx = uTargetIdx(ismember(uTargetIdx, TargetsSortForFig5.([zc{tf,1},'_',zc{tf,2}])));
    connMat=(perOverlap(uTargetIdx,uTargetIdx)>0.2) | totalOL(uTargetIdx,uTargetIdx);
    connGraph=graph(connMat);
    connBins=conncomp(connGraph);
    [~,keepIdx]=unique(connBins);
    uTargetIdx=uTargetIdx(keepIdx);
    
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
            sumSignalPattern(g,p,1) = sum(checWTdelLactisSwap.norm.(TF1){GP.gene_infoR64.position(uTargetIdx(g),1)}(intPromPos))./sum(checWTdelLactisSwap.norm.(TF1){GP.gene_infoR64.position(uTargetIdx(g),1)}(promPos));
            sumSignalPattern(g,p,2) = sum(checWTdelLactisSwap.norm.(TF2){GP.gene_infoR64.position(uTargetIdx(g),1)}(intPromPos))./sum(checWTdelLactisSwap.norm.(TF2){GP.gene_infoR64.position(uTargetIdx(g),1)}(promPos));
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
        for r=find(all(isnan(sumPromRep)))
            sumPromRep(:,r) =checWTdelLactisSwap.sumProm.(intStrains{idxVec(r)});
        end
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
        zscore_WT1 = nanZscore(checWTdelLactisSwap.sumProm.(currF1));
        fitGenes =  zscore_WT1 > min(zscoreTH, quantile(zscore_WT1, quantileTH));
        for ct = 1:numel(compType)
            m = allSamples(startsWith(allSamples,currF1)&contains(allSamples, compType{ct}));
            if numel(m) > 0
                m=m{1};
                bestFit = robustfit(checWTdelLactisSwap.sumProm.(currF1)(fitGenes), checWTdelLactisSwap.sumProm.(m)(fitGenes), [], [],'off' );
                normFactor(par) = log2(bestFit);
                logChange{par}(ct,:) = log2(checWTdelLactisSwap.sumProm.(m)+700)- log2(checWTdelLactisSwap.sumProm.(currF1)+700)- normFactor(par);
                %                 significantGenes(par,:) =  abs((logChange(par,:)-median(logChange(par,fitGenes)))/std(logChange(par,fitGenes))) >= 1;
            else
                logChange{par}(ct,:) = nan(1,6701);
            end
        end        
    end
    logChange=cat(1,logChange{:});
    logChangeSelectedTargets =  logChange([1,3,2,4],uTargetIdx);
   % cluster by targets
   if ismember([examplePairs{tf,1},'_',examplePairs{tf,2}],TargetsSortForFig5.Properties.VariableNames)
       [~,idx] = ismember(TargetsSortForFig5.([examplePairs{tf,1},'_',examplePairs{tf,2}]), uTargetIdx);
       idx(idx==0) = [];
   else
       % cluster by targets
       clusterIdx = sum((zScoreSelectedTargetsMean(:,[1,2])> zscoreTHclus).*[1,2],2)+1;
       clusterReorder = [2,1,3,2];
       clusterIdx = clusterReorder(clusterIdx);
       sortBy = [1,1,2,2];
       valForSort = nan(size(clusterIdx));
       for i = 1:max(clusterIdx)
           valForSort(clusterIdx==i) = -zScoreSelectedTargetsMean(clusterIdx==i,sortBy(i));
       end
       [~, idx] = sortrows(table(clusterIdx',valForSort'), [1,2]);
   end

    % imagesc Zscores    
    axes('Position', [xPos(tf)+0.54 yPos(tf) Wline*3 Hline])
    imagesc(zScoreSelectedTargets(idx,:))
    hold on
    plot([maxSize 2*maxSize]+[0.5;.5],ylim,'k-')
    plot(repmat(xlim',1,numel(idx)-1),repmat(1.5:numel(idx),2,1),'k-')
    caxis([0 10])
    set(gca,'YTick', 1:numel(idx),  'XTick',[(maxSize+1)/2 maxSize+(maxSize+1)/2 maxSize*2+(maxSize/2+1)/2], 'XTickLabel', {TF1,TF2,'klac'}, 'fontSize',11, 'TickLength',[0 0],'YTickLabel',sprintfc('%d',uTargetIdx(idx)));
    colormap(gca,cMapZscores)    
   
    % imagesc log2 fold change
    XL={'DeltaP2','DeltaP1','S1','S2'};
    axes('Position', [xPos(tf)+0.54+Wline*3.5 yPos(tf) Wline*2 Hline])
    imagesc(logChangeSelectedTargets(:,idx)')
    plotgrid(logChangeSelectedTargets(:,idx)')
    caxis([-2 2])
    set(gca, 'XTick', 1:numel(XL), 'XTickLabel', XL,'YTick',[])
    colormap(gca,cMapLogChange)
    if any(all(isnan(logChangeSelectedTargets([3,4],:)),2),1)
        xlim([.5 2.5])
    end
        if mod(tf,2)==0
%            save_gf(gcf,sprintf('FigS6_TF%d',tf),'paper','tamar21','type','svg','size',[])
            %close(gcf)
        end
end

% colorbars:
axes('Position', [xPos(tf)+3*(Wscatter+xspacer), yPos(tf)-yspacer*2, Wscatter/3, 0.02])
caxis([0 10])
colormap(gca,cMapZscores)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'z-score', 'fontSize',10)
set(cbr, 'Ticks', [0.1+min(caxis),max(caxis)-1] , 'TickLabels', caxis)

axes('Position', [xPos(tf)+3*(Wscatter+xspacer)+Wscatter/3, yPos(tf)-yspacer*2, Wscatter/3, 0.02])
caxis([-2 2])
colormap(gca,cMapLogChange)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'Fold change    (log2)', 'fontSize',10)
set(cbr, 'Ticks', [0.5+min(caxis),0,max(caxis)-0.1] , 'TickLabels', [min(caxis),0,max(caxis)],'TickLength',[0 0])
%save_gf(gcf,sprintf('FigS6_TF%d',tf),'paper','tamar21','type','svg')

axes('Position', [xPos(tf)+3*(Wscatter+xspacer)+2*Wscatter/3, yPos(tf)+Hscatter+yspacer-yspacer-4*smallYspacer-9*Hline, Wscatter/3, 0.02])
caxis([0 4]+0.5)
colormap(gca,cMapMotifs)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'#motifs', 'fontSize',10)
set(cbr, 'Ticks', [min(caxis)+0.5:1:max(caxis)-0.5] , 'TickLabels', {[0:max(caxis)-2.5], '>2'},'TickLength',[0 0])


%% Figure 9—figure supplement 1A: alignment of domains for Rph1 and Gis1
multiAlign = multialignread('./mCoffee/GIS1_fl.clustal')
jmjcPos = [510:538];
DBDPos = [1687:1745];

for i = 1:numel(multiAlign)
    alignForPlot(i).Sequence = [multiAlign(i).Sequence(jmjcPos), '---', multiAlign(i).Sequence(DBDPos)];
    alignForPlot(i).Header = multiAlign(i).Header;
end
alignForPlot=alignForPlot([19:28,11:18,1:10])

imageAlignDBD({'Gis1','Rph1'}, [],'alignment', alignForPlot)
pc = get(gca,'Children');
imMat=pc(end);
imMat.CData(:, numel(jmjcPos)+[1:3]) = 22;
currColormap = get(gca,'Colormap');
colormap([currColormap;1 1 1])
ax=gca;
cMap=ax.Colormap
cMap(14:20,:)=repmat(cMap(14,:),7,1);
cMap(6:13,:)=repmat(cMap(6,:),8,1);
colormap(cMap)
ax.Colormap=brighten(ax.Colormap,.9)
ax.Colormap(3:5,:)=repmat([212,233,255]/255,3,1)
intPos=[10,41] 
imMat.CData(:,intPos)=imMat.CData(:,intPos)+22
ax.Colormap=[ax.Colormap;cMap]
caxis([0.5 44.5])
save_gf(gcf,'FigS6Align','paper','tamar21','type','svg')
