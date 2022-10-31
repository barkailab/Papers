%% Figure 8 
clearvars -except checWTdelLactisSwap
load('paraSeqs.mat')
load('summaryTable.mat')
load('promoterIDXvec.mat')
allDomain = readtable('./allDomainsAdj.txt');
intDomains = {'Zn_clus','zf-C2H2','Forkhead','KilA-N','bZIP_1','Vhr1','Homeodomain','Myb_DNA-binding','SRF-TF',...
    'AFT','GATA','Copper-fist','HSF_DNA-bind','zf-MIZ','HMG_box','TIG','zf-BED'};
allDBD = allDomain(ismember(allDomain.x_TargetName,intDomains),:);
KlacDBDs = allDBD(contains(allDBD.queryName,'Klac'),:);
ScerDBDs = allDBD(contains(allDBD.queryName,'Scer'),:);
ScerDBDs = ScerDBDs(ScerDBDs.score_1>0,:);
ccPredictions = load('ccPredictions.mat');
CCpredictionsWJpred = load('./CCpredictionsWJpred.mat')
ccTH = 0.2; 
colRangeCC = [ccTH:0.01:1];
%ccColMap = cbrewer2('PuBu', length(colRangeCC)); 
ccColMap = brewermap( length(colRangeCC),'PuBu'); 
load('iupred2.mat')
h3cc=load('./H3CC_henikoff.mat');
h3smooth=conv(mean(h3cc.p_all,2),normpdf([-50:50]',0,25),'same');
load('SC_genome.mat')
SC_genomeCat=upper(cat(2,SC_genome.Sequence));
GP=load('./group_imp.mat');
if ~exist('checWTdelLactisSwap','var')
    load('checWTdelLactisSwap.mat')
end
load('promoterLengthsORF.mat')
% creating vector with the mid position of all promoters
promStartPos=GP.chrIdx(min(GP.gene_infoR64.position(:,1),18))+GP.gene_infoR64.position(:,2);
promStartPos(isnan(promoterLengthsORF))=NaN;
promMidVec = nan(6701,1);
promMidVec(~isnan(promoterLengthsORF)) = GP.chrIdx(GP.gene_infoR64.position(~isnan(promoterLengthsORF),1))+GP.gene_infoR64.position(~isnan(promoterLengthsORF),2)-promoterLengthsORF(~isnan(promoterLengthsORF)).*GP.gene_infoR64.dir(~isnan(promoterLengthsORF));
disMat = abs(promMidVec-promMidVec');
meanPromoterLengthMat = (promoterLengthsORF+promoterLengthsORF')/2;
overlapMat = meanPromoterLengthMat-disMat;
overlapMat(overlapMat<0) = 0;
perOverlap = overlapMat./min(promoterLengthsORF,promoterLengthsORF');
perOverlap(1:6702:end) = nan;
totalOL=(sign(promMidVec-promMidVec')+sign(promStartPos-promStartPos'))==0;
totalOL(1:6702:end)=0;

%% Figure 8A
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y','X', '-'};
B62 = blosum(62,'Order', cat(2,AAorder{1:end-1}));
B62(22,:) = -8;
B62(:,22) = -8;
colMapAlignment =brighten([1 1 1; flipud(bone)],0.6);
%DBDboxCol = [135 48 146]/255;
DBDboxCol = [0 230 230]/255;

% CCpredictionsWJpred = load('./CCpredictionsWJpred.mat')
% load('iupred2.mat')
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
zc = {'Tbs1', 'Hal9';'Oaf1','Pip2';'Pdr1','Pdr3'; 'Yrr1','Pdr8'}%; 'Skn7','Hms2'};
TargetsSortForFig5 = readtable('./TargetsSortForFig8.xlsx');
targets = { 'YCR102C','TPO1' ; 'ACS1','EEB1'; 'RPN4','PDR3'; 'SNQ2','RSB1'}
patternLoc = [-310, -612; -419, -164; -438, -219; -516,-811];

patternsColl{1} =  {'CGG.TCCG|CGGA.CCG',  'C[GC]G...CGG|CCG...C[GC]G|CGG..CGG|CCG..CCG'};
patternsColl{2} = {'CGG.{2,5}CGG|CCG.{2,5}CCG', 'CGG.{15,16}CCG'};
patternsColl{3} = {'TCCGTGG|CCACGGA','TCCGCGG[AG]|[TC]CCGCGGA'};
patternsColl{4} = {'CCG.{6,11}CC[AG]|[TC]GG.{6,11}CGG',...
    '[TC]CCGCGG[AC]|[GT]CCGCGG[GA]'};

maxL = max(summaryTable.proteinLength(contains(summaryTable.p1, zc),:),[],'all');
intoGene = 100;
promoterLengthPlot = 900;
promoterLengthPlotVec = repmat(promoterLengthPlot,6701,1);
intBases = intBasesVec(1:6701);

xspacer = 0.005;
yspacer = 0.04;
Wsignal = 0.05;
Hsignal = 0.2;
Wpromoter = 0.12;
Hpromoter = (Hsignal-yspacer)/2;
xPos = [0.2:xspacer+Wsignal:1-Wsignal];
yPos = [0.75:-(Hsignal+yspacer):0];
Wscheme = Wsignal*2.5;
Hscheme = Hsignal/3.5;
Wscatter = (Hsignal-Hscheme-yspacer)/1.8;
Hscatter = Wscatter*1.8;
Hzoom = Hpromoter/1.25;
Wzoom = Wpromoter;

colSignal = repmat([0 0 0; 0.5 0.5 0.5], 2,1);
colMotifs = lines(6);
caxisLim = [0 0.2];
cMapZscores = flipud(bone);
% cMapMotifs = brighten(cbrewer2('Purples',4),-0.2);
%cMapMotifs = brighten(cbrewer2('Oranges',8),-0.2);
cMapMotifs = brighten(brewermap(8,'Oranges'),-0.2);

cMapMotifs = cMapMotifs([1:4],:);
cMapMotifs(1,:) = [1,1,1];
%cMapMotifs2 = brighten(cbrewer2('PuRd'),-0.4);
cMapMotifs2 = brighten(brewermap(256,'PuRd'),-0.4);
cMapMotifs2(1,:) = [1,1,1];
%cMapLogChange = flipud(cbrewer2('RdBu'));
cMapLogChange = flipud(brewermap(256,'RdBu'));
%cMapClusters = cbrewer2('Set3',4);
cMapClusters = [0.7 0.7 0.7; 0.2 0.2 0.2];
N = 30;
zscoreTH = 2;
quantileTH = 0.98;
zscoreTHclus = 4;

figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
for tf = 1:size(zc,1)
    targetList = targets(tf,:);
    TF1 = zc{tf,1};
    TF2 = zc{tf,2};
    m1 = [TF1, '_d',upper(TF2)];
    m2 = [TF2, '_d',upper(TF1)];
    intStrains = {TF1,m1,TF2,m2};
    intStrains = intStrains(ismember(intStrains,allSamples));
    pattern = patternsColl{tf};
    
    % protein scheme: Disorder tendency and DBD annotations
    clear DBDpos proteinSeq disTenVec
    for j =1:2
        DBDpos{j} = [ScerDBDs.from_2(find(contains(ScerDBDs.queryName, upper(zc{tf,j})))), ScerDBDs.to_2(find(contains(ScerDBDs.queryName, upper(zc{tf,j}))))];
        disTenVec{j} = iupred2.score{strcmp(iupred2.commonName, ['Scer_', upper(zc{tf,j})])};
        paraIdx = find(contains({paraSeqs.name}, {upper(zc{tf,1}), upper(zc{tf,2})}));
        currProteinSeq =  paraSeqs(paraIdx).seq{strcmp(paraSeqs(paraIdx).commonName,  ['Scer_', upper(zc{tf,j})])};
        proteinSeq{j} =  strrep(currProteinSeq,'*','');
    end
    [~, alignment] = nwalign(proteinSeq{1}, proteinSeq{2}) ;
    paraLength = [size(proteinSeq{1},2), size(proteinSeq{2},2)];
    currMaxL = max(paraLength);
    
    clear posScore gaps
    [~, idxPara1] = ismember(alignment(1,:)', AAorder);
    [~, idxPara2] = ismember(alignment(3,:)', AAorder);
    gaps{1} = strfind(alignment(1,:), '-');
    gaps{2} = strfind(alignment(3,:), '-');

    posScore = B62(sub2ind(size(B62), idxPara1,idxPara2));
    posScoreNoGaps{1} = posScore(setdiff(1:length(posScore), gaps{1}));
    posScoreNoGaps{2} = posScore(setdiff(1:length(posScore), gaps{2}));

    imageMat = nan(2, currMaxL);
    for j =1:2
        imageMat(j,1:length(posScoreNoGaps{j})) = posScoreNoGaps{j};
    end
    imageMat = movmean(imageMat,20,2, 'omitnan').*((imageMat+50)./(imageMat+50));
    
    for j = 1:2
        axes('Position', [0.05, yPos(tf)+Hsignal-Hscheme/2-(j-1)*(Hscheme/2+yspacer/1.2), Wscheme,  Hscheme/2])
        imagesc(imageMat(j,:))
        hold on
        xVec = [DBDpos{j}, DBDpos{j}([2,1]), DBDpos{j}(1)];
        yVec = [1.5, 1.5, 0.5, 0.5, 1.5];
        plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
        plot([1:paraLength(j)], movmean(rescale(1-disTenVec{j}, 0.5, 1.5),5), 'color','k')
        
        set(gca, 'YTick', 1, 'YTickLabel', sprintf('%s_{P%d}', zc{tf,j}, j), 'fontSize',12, 'XTick',[], 'TickLength',[0 0],  'XColor', [1 1 1 0])
        text(paraLength(j)+5, 1.7, sprintf('%d aa', paraLength(j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'fontSize',11)
        plot(paraLength(j)*[1 1], [0.5, 1.5], 'k');
        plot([0,paraLength(j)]+0.5, [0.5 0.5], 'k');
        plot([1:paraLength(j)], repmat(1,paraLength(j),1), '--k', 'LineWidth', 0.5)
        plot([1:paraLength(j)], repmat(1.5,paraLength(j),1) , 'k')
        colormap(gca, colMapAlignment);
        caxis([-4 4]);
        box off;
        text(paraLength(j)+5, 0.5, '1', 'HorizontalAlignment', 'left', 'fontSize',10,'VerticalAlignment','middle')
        text(paraLength(j)+5, 1.5, '0', 'HorizontalAlignment', 'left', 'fontSize',10,'VerticalAlignment','middle')
    end
    
    %     axes('Position', [0.05+Wscheme+4*xspacer, yPos(tf)+Hsignal-Hscheme-(Hscheme/2+yspacer/1.5), Wscheme,  Hscheme+yspacer/1.5])
    %     text(0.1, 0.5, 'Disorder', 'Rotation',270, 'HorizontalAlignment', 'center', 'fontSize',11)
    %     text(0, 0.5, 'tendency', 'Rotation',270, 'HorizontalAlignment', 'center', 'fontSize',11)
    %     axis off
           
    % Top targets: pattern occurences +zscore mat + log2 fold change
    clear patternCorrMat
    k=7;
    uTargetIdx = [];
    idxMax={};
    while numel(uTargetIdx)<N
        [~,idxMax{1}] = maxk(checWTdelLactisSwap.sumProm.(TF1),k);
        [~,idxMax{2}] = maxk(checWTdelLactisSwap.sumProm.(TF2),k);
        uTargetIdx = unique(cat(1, idxMax{:})', 'stable');
        k=k+1;
    end
    % takong out overlapping promoters
    [~,selProms]=ismember(targetList,GP.gene_infoR64.nameNew);
    uTargetIdx=[uTargetIdx(ismember(uTargetIdx,selProms)),uTargetIdx(~ismember(uTargetIdx,selProms))];
    %uTargetIdx = uTargetIdx(ismember(uTargetIdx, TargetsSortForFig5.([zc{tf,1},'_',zc{tf,2}])));
    connMat=(perOverlap(uTargetIdx,uTargetIdx)>0.2) | totalOL(uTargetIdx,uTargetIdx);
    connGraph=graph(connMat);
    connBins=conncomp(connGraph);
    [~,keepIdx]=unique(connBins);
    uTargetIdx=uTargetIdx(keepIdx);
    %[ui, uj] = find(tril(perOverlap(uTargetIdx,uTargetIdx)>0.2));
    %uTargetIdx(uj) = [];
   
    % pattern:
    clear nPatternInPromoter sumSignalPattern
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
            if tf == 4 & p==1
                currMatches = regexp(currMatches, 'CCG...*?CC[AG]|[TC]GG...*?CGG', 'match','once');
                linker = cellfun(@(x)x(4:end-3), currMatches, 'UniformOutput', false);
                ATcontent = cellfun('prodofsize', regexp(linker, 'A|T'));
                keepPattern = ATcontent./(cellfun('prodofsize', currMatches)-6) >0.6;
                currMatches = currMatches(keepPattern);
                S = S(keepPattern);
                E = E(keepPattern);
            end
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
            %sumSignalPattern(g,p,1)=sumSignalPattern(g,p,1)/numel(M);
            %sumSignalPattern(g,p,2)=sumSignalPattern(g,p,2)/numel(M);
            forProm(g).pos{p}=GP.chrIdx(GP.gene_infoR64.position(uTargetIdx(g),1))+promPos;
            forProm(g).intPos{p}=GP.chrIdx(GP.gene_infoR64.position(uTargetIdx(g),1))+intPromPos;
        end
    end
    
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
        currF1 = zc{tf,par};
        m = [currF1,'_d',upper(zc{tf,3-par})];
        zscore_WT1 = nanZscore(checWTdelLactisSwap.sumProm.(currF1));
        fitGenes =  zscore_WT1 > min(zscoreTH, quantile(zscore_WT1, quantileTH));
        if any(ismember(allSamples,m))
            bestFit = robustfit(checWTdelLactisSwap.sumProm.(currF1)(fitGenes), checWTdelLactisSwap.sumProm.(m)(fitGenes), [], [],'off' );
            normFactor(par) = log2(bestFit);
            normFactorAllZC(tf,par) = log2(bestFit);
            logChange(par,:) = log2(checWTdelLactisSwap.sumProm.(m)+700)- log2(checWTdelLactisSwap.sumProm.(currF1)+700)- normFactor(par);
            significantGenes(par,:) =  abs((logChange(par,:)-median(logChange(par,fitGenes)))/std(logChange(par,fitGenes))) >= 1;
        else
            logChange(par,:) = nan(1,6701);
        end
    end
    logChangeSelectedTargets =  logChange(:,uTargetIdx);
    notNanFoldChange = ~isnan(logChangeSelectedTargets);
    
    % cluster by targets
    [~,idx] = ismember(TargetsSortForFig5.([zc{tf,1},'_',zc{tf,2}]), uTargetIdx);
    idx(idx==0) = []; 
    
    axes('Position', [xPos(1), yPos(tf), Wsignal, Hsignal])
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
    axes('Position', [xPos(1)-0.01, yPos(tf), 0.01, Hsignal])
    targetsPosPlot = find(contains(GP.gene_infoR64.name(uTargetIdx(idx)), targetList));
    scatter([0.5 0.5], targetsPosPlot, 15,cMapClusters, 'filled')
    ylim([0.5 numel(idx)+0.5])
    xlim([0 1])
    set(gca, 'YDir', 'reverse')
    axis off
    
    % motifs 
    temp=getMotifScoresRepeats(forProm(idx),intStrains);
    axes('Position', [xPos(1)+xspacer+Wsignal, yPos(tf), Wsignal, Hsignal])
    [~,idxForPattern] = max( zScoreSelectedTargetsMean(idx,:),[],2);
    idxForPattern=sub2ind(size(temp),[1:size(temp,1)]',idxForPattern)
    %idxForPattern=idxForPattern';
    sumSignalPatternPlot = cell2mat(cat(1,temp{idxForPattern}))%[sumSignalPattern(sub2ind(size(sumSignalPattern), 1:size(sumSignalPattern,1), repmat(1,1, size(sumSignalPattern,1)), idxForPattern)); sumSignalPattern(sub2ind(size(sumSignalPattern), 1:size(sumSignalPattern,1), repmat(2,1, size(sumSignalPattern,1)), idxForPattern))]';
    imagesc(sumSignalPatternPlot)
    hold on
    plot(repmat(xlim',1,numel(idx)-1),repmat(1.5:numel(idx),2,1),'k-')
    plot(numel(temp{1}{1})+[0.5;0.5],ylim,'k-')

    %plotgrid(sumSignalPatternPlot)
    
    set(gca,'YTick',[],'YtickLabels',uTargetIdx(idx),  'XTick', [1:2], 'XTickLabel', {'mA ', 'mB'} ,'fontSize',11,'TickLength',[0 0]);
    colormap(gca,cMapMotifs2(1:size(cMapMotifs2,1)/2, :))
    caxis([0 0.5])
    if tf==1
        text(1.5,-3,'Motif', 'fontSize',12, 'HorizontalAlignment', 'center')
        text(1.5,-1,'score', 'fontSize',12, 'HorizontalAlignment', 'center')
    end
    
    axes('Position', [xPos(1)+2*(xspacer+Wsignal), yPos(tf), Wsignal, Hsignal])
    imagesc(logChangeSelectedTargets(find(all(logChangeSelectedTargets,2)),idx)')
    plotgrid(logChangeSelectedTargets(find(all(logChangeSelectedTargets,2)),idx)')
    caxis([-2 2])
    if any(contains(intStrains(find(all(logChangeSelectedTargets,2))), 'Pip2'))
        set(gca, 'XTick', [1,2], 'XTickLabel', strcat(intStrains(find(all(logChangeSelectedTargets,2))),{'','*'}),  'fontSize',11, 'YTick',[], 'TickLength',[0 0])
    else
        set(gca, 'XTick', [1,2], 'XTickLabel', intStrains(find(all(logChangeSelectedTargets,2))), 'fontSize',11, 'YTick',[], 'TickLength',[0 0])
    end
    colormap(gca, cMapLogChange)
    if tf==1
        text(1.5,-3,'Binding change', 'fontSize',12, 'HorizontalAlignment', 'center')
            text(1.5,-1,'(DeltaParalog/wt)', 'fontSize',12, 'HorizontalAlignment', 'center')
    end    
  
% signal on example promoters
    intStrains = {TF1,m1,TF2,m2};
    intStrains = intStrains(ismember(intStrains,allSamples));
    for g = 1: length(targetList)
        if ~isempty(targetList{g})
            axes('Position', [xPos(1)+5*(Wsignal+xspacer)+(g-1)*(Wpromoter+xspacer) yPos(tf)+Hsignal-Hpromoter-yspacer/4 Wpromoter Hpromoter])
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
            hold on
            for z = 1:numel(intStrains)
                plot([-promoterLengthPlotVec(geneIdx):intoGene], rescale(signalMat(z,:),numel(intStrains)-z+1.5, numel(intStrains)-z+2.3,'InputMin',0, 'InputMax',0.2),...
                    'color', colSignal(z,:))
            end
            plot(-promoterLengthPlotVec(geneIdx):intoGene, rescale(h3smooth(GP.chr_idx(pos(1))+intBasesGene), 0.5, 1.3),'k', 'LineWidth',1)
            area(-promoterLengthPlotVec(geneIdx):intoGene, rescale(h3smooth(GP.chr_idx(pos(1))+intBasesGene), 0.5, 1.3),'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'LineStyle', 'none')
            ylim([0.5 numel(intStrains)+1.5])
            xlim([-promoterLengthPlotVec(geneIdx) intoGene])

            TssLoc = (GP.gene_infoR64.tamarTss(geneIdx,2)-pos(2))*GP.gene_infoR64.dir(geneIdx);
            labelsSignal = {'P1','P1DeltaP2', 'P2','P2DeltaP1'};
            %set(gca, 'YTick', [1:numel(labelSamples)], 'YTickLabel', labelSamples(end:-1:1), 'fontSize',12, 'XTick', [-800:200:-400, TssLoc,100], 'XTickLabel', {'-800','-600','-400','TSS', '100'})
            if g ==1
                set(gca, 'YTick', [1:numel(labelSamples)+1], 'YTickLabel', {'NucOcc', labelsSignal{numel(labelSamples):-1:1}},'XTick',[], 'fontSize',10, 'TickLength',[0 0])
            else
                set(gca, 'YTick', [], 'XTick',[], 'TickLength',[0 0])
            end
            
            plot(TssLoc*[1 1] ,ylim, ':r', 'LineWidth',1.5)
            scatter(0 ,0.5,[], '>r','filled')
            text(TssLoc, 0.5, 'TSS', 'color', 'r','HorizontalAlignment','center','VerticalAlignment','top','fontSize',8);
            if ~isempty(motifPatternLocSAll)
                for p = 1:length(patternIdx)
                    %scatter(motifPatternLocSAll(p), 0.5,50, [colMotifs(patternIdx(p),:)],'filled', 'MarkerEdgeColor','k');
                    fill([motifPatternLocSAll(p)*[1,1], motifPatternLocEAll(p)*[1,1]], [0.5, numel(intStrains)+1.5.*[1,1], 0.5], colMotifs(patternIdx(p),:),...
                        'LineStyle','none', 'FaceAlpha', 0.3)
                end
            end
            title(sprintf('%s',target),'fontSize',9)
            box on
         
            % zoom on a single pattern
            axes('Position', [xPos(1)+5*(Wsignal+xspacer)+(g-1)*(Wpromoter+xspacer) yPos(tf) Wzoom Hzoom])
            hold on
            for z = 1:numel(intStrains)
                plot([-promoterLengthPlotVec(geneIdx):intoGene], rescale(signalMat(z,:),numel(intStrains)-z+1.5, numel(intStrains)-z+2.3,'InputMin',0, 'InputMax',0.2),...
                    'color', colSignal(z,:))
            end
            ylim([1.5 numel(intStrains)+1.5])
            xlim([patternLoc(tf,g)+[-50 50]])
            if ~isempty(motifPatternLocSAll)
                for p = 1:length(patternIdx)
                    %scatter(motifPatternLocSAll(p), 0.5,50, [colMotifs(patternIdx(p),:)],'filled', 'MarkerEdgeColor','k');
                    fill([motifPatternLocSAll(p)*[1,1], motifPatternLocEAll(p)*[1,1]], [0.5, numel(intStrains)+1.5.*[1,1], 0.5], colMotifs(patternIdx(p),:),...
                        'LineStyle','none', 'FaceAlpha', 0.3)
                end
            end
            if g==1
                %set(gca, 'YTick', [1:numel(labelSamples)], 'YTickLabel', labelSamples(end:-1:1), 'fontSize',12, 'XTick', [-800:200:-400, TssLoc,100], 'XTickLabel', {'-800','-600','-400','TSS', '100'})
                set(gca, 'YTick', [1:numel(labelSamples)+1], 'YTickLabel', {'NucOcc', labelsSignal{numel(labelSamples):-1:1}},'XTick',[patternLoc(tf,g)+[-45 45]],'XTickLabel',[], 'fontSize',11)
            else
                set(gca,'YTick',[], 'XTick',[patternLoc(tf,g)+[-45 45]], 'XTickLabel', [])
            end
            text(patternLoc(tf,g)+[-45 45], [1.5,1.5], sprintfc('%d', patternLoc(tf,g)+[-45 45]'), 'fontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment','top')
            box on
        end
    end
end
%save_gf(gcf,sprintf('IndRepeats_Matrix_5'),'type',{'svg'},'paper','tamar21')


%% Figure 8B-D - Ecm22/Upc2 sequence alignment and WTs, deletion mutants and DBD-swap variants scatters (showing the change in binding relative to the wild-type)
examplePairs = {'Upc2','Ecm22'};
targets = {'ERG3','UPC2'};
patternsColl{1} = {'TA[TA]ACGA|TCGT[TA]TA'};
compType = {'_d','_DBD'};
promoterLengthPlot = 900;
promoterLengthPlotVec = repmat(promoterLengthPlot,6701,1);
kLacNtargets = 50;

colMapAlignment =brighten([1 1 1; flipud(bone)],0.6);
DBDboxCol = [0 230 230]/255;

xspacer = 0.04;
yspacer = 0.07;
xPos = 0.1;
yPos = 0.5;
Hsignal = 0.15;
Wsignal =  0.15;
Wscatter = Wsignal;
Hscatter = Wscatter*1.6;
smallYspacer = 0.01;
Hline = (Wsignal-3*(smallYspacer))/7.5;
Wscheme = 0.125;

%colSignal = repmat([0 0 0; 0.5 0.5 0.5], 2,1);
colSignal = repmat([0 0 0],7,1);
colMotifs = lines(6);
caxisLim = [0 0.2];
%cMapLac = cbrewer2('OrRd');
cMapLac = brewermap(256,'OrRd');
cMapZscores = flipud(bone);
% cMapMotifs = brighten(cbrewer2('Purples',4),-0.2);
% cMapMotifs(1,:) = [1,1,1];
%cMapMotifs = brighten(cbrewer2('Oranges',8),-0.2); %brighten(brewermap(8,'Oranges'),-0.2);
cMapMotifs = brighten(brewermap(8,'Oranges'),-0.2); %brighten(brewermap(8,'Oranges'),-0.2);
cMapMotifs = cMapMotifs([1:4],:);
cMapMotifs(1,:) = [1,1,1];
%lineCol = [cbrewer2('Blues',1);cMapLac(150,:)];
%ccColMap = cbrewer2('PuBu', length(colRangeCC)); 
lineCol = [brewermap(1,'Blues');cMapLac(150,:)];
ccColMap = brewermap(length(colRangeCC),'PuBu'); 

for tf = 1:size(examplePairs,1)
    figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
    targetList = targets(tf,:);
    TF1 = examplePairs{tf,1};
    TF2 = examplePairs{tf,2};
    m1 = [TF1, '_d',upper(TF2)];
    m2 = [TF2, '_d',upper(TF1)];
    s1 = [TF1,'_',TF2,'_DBD'];
    s2 = [TF2,'_',TF1,'_DBD'];
    lacIdx = find(contains(allSamples, 'lactis') & contains(allSamples, {TF1,TF2}));
    lac = allSamples{lacIdx};
    intStrains = {TF1,m1,s1,TF2,m2,s2,lac};
    intStrains = intStrains(ismember(intStrains,allSamples));
    pattern = patternsColl{tf};
    
     % protein scheme: DBD+ coiled coil
     paraIdx = find(contains({paraSeqs.name}, upper(examplePairs)));
     clear DBDpos proteinSeq disTenVec
     for j =1:2
         DBDpos{j} = [ScerDBDs.from_2(find(contains(ScerDBDs.queryName, upper(examplePairs{j})))), ScerDBDs.to_2(find(contains(ScerDBDs.queryName, upper(examplePairs{j}))))];
         disTenVec{j} = iupred2.score{strcmp(iupred2.commonName, ['Scer_', upper(examplePairs{j})])};
         currProteinSeq =  paraSeqs(paraIdx).seq{strcmp(paraSeqs(paraIdx).commonName,  ['Scer_', upper(examplePairs{j})])};
         proteinSeq{j} =  strrep(currProteinSeq,'*','');
     end
     kLacSeq=strrep(paraSeqs(paraIdx).seq{18},'*','');
     disTenVec{3} = iupred2.score{strcmp(iupred2.commonName, paraSeqs(paraIdx).commonName{18})};
     [~, alignment] = nwalign(proteinSeq{1}, proteinSeq{2}) ;
     paraLength = cellfun('prodofsize',[proteinSeq,kLacSeq])%[size(proteinSeq{1},2), size(proteinSeq{2},2)];
     currMaxL = max(paraLength);
     
     clear posScore idxPara
     [~, idxPara(1,:)] = ismember(alignment(1,:)', AAorder);
     [~, idxPara(2,:)] = ismember(alignment(3,:)', AAorder);
     posScore = B62(sub2ind(size(B62), idxPara(1,:),idxPara(2,:)));
     imageMat = nan(2, currMaxL);
     for j =1:2
         imageMat(j,1:sum(idxPara(j,:)~=22)) = posScore(idxPara(j,:)~=22);
     end
     
     % klac alignments
     DBDposLac = [KlacDBDs.from_2(find(contains(KlacDBDs.queryName, upper(examplePairs{tf,1})))), KlacDBDs.to_2(find(contains(KlacDBDs.queryName, upper(examplePairs{tf,1}))))];
     kLacMat= nan(2, currMaxL);
     for j=1:2
         [klacS(j), alignment] = nwalign(proteinSeq{j}, kLacSeq) ;
         clear posScore idxPara
         [~, idxPara(1,:)] = ismember(alignment(1,:)', AAorder);
         [~, idxPara(2,:)] = ismember(alignment(3,:)', AAorder);
         posScore = B62(sub2ind(size(B62), idxPara(1,:),idxPara(2,:)));
         kLacMat(j,1:sum(idxPara(2,:)~=22))=posScore(idxPara(2,:)~=22);
     end
     imageMat = movmean(imageMat,20,2, 'omitnan').*((imageMat+50)./(imageMat+50));
     kLacMat = movmean(kLacMat,20,2, 'omitnan').*((kLacMat+50)./(kLacMat+50));
     
     axes('Position', [xPos-0.02 yPos+Hscatter/2 Wscheme Hscatter/2])
     hold off
     imDis=2;
     for j=1:2
         imagesc(imageMat(j,1:paraLength(j)),'YData',j*imDis)
         hold on
         xVec = [DBDpos{j}, DBDpos{j}([2,1]), DBDpos{j}(1)];
         yVec = j*imDis+[.5, .5, -.5, -.5, .5];
         plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
         plot([1:paraLength(j)], movmean(rescale(1-disTenVec{j},j*imDis-.5, j*imDis+.5),5), 'color','k')         
         text(paraLength(j)+5, j*imDis+0.6, sprintf('%d aa', paraLength(j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'fontSize',11)
         plot(paraLength(j)+[.5 .5], j*imDis+[-.5, .5], 'k');
         plot(0+[.5 .5], j*imDis+[-.5, .5], 'k');
         plot([0.5,paraLength(j)]+0.5, j*imDis-[0.5 0.5], 'k');
         plot([0.5,paraLength(j)]+0.5, j*imDis+[0.5 0.5], 'k');         
         plot([0.5,paraLength(j)]+0.5, j*imDis.*[1 1], '--k', 'LineWidth', 0.5)
         text(paraLength(j)+5, j*imDis-0.5, '1', 'HorizontalAlignment', 'left', 'fontSize',10,'VerticalAlignment','middle')
         text(paraLength(j)+5, j*imDis+0.5, '0', 'HorizontalAlignment', 'left', 'fontSize',10,'VerticalAlignment','middle')         
     end
    imagesc(kLacMat(2,1:paraLength(3)),'YData',3*imDis)
    plot(paraLength(3)+[.5 .5], 3*imDis+[-.5, .5], 'k');
    plot(0+[.5 .5], 3*imDis+[-.5, .5], 'k');
    plot([0,paraLength(3)]+0.5, 3*imDis-[0.5 0.5], 'k');
    plot([0,paraLength(3)]+0.5, 3*imDis+[0.5 0.5], 'k');
    plot([0.5,paraLength(3)]+0.5, 3*imDis.*[1 1], '--k', 'LineWidth', 0.5)
    text(paraLength(3)+5, 3*imDis+0.6, sprintf('%d aa', paraLength(3)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'fontSize',11)
    text(paraLength(3)+5, 3*imDis-0.5, '1', 'HorizontalAlignment', 'left', 'fontSize',10,'VerticalAlignment','middle')
    text(paraLength(3)+5, 3*imDis+0.5, '0', 'HorizontalAlignment', 'left', 'fontSize',10,'VerticalAlignment','middle')
    xVec = DBDposLac([1,2,2,1,1]);
    yVec =3*imDis+[.5, .5, -.5, -.5, .5];
    plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
    plot([1:paraLength(3)], movmean(rescale(1-disTenVec{3},3*imDis-.5, 3*imDis+.5),5), 'color','k')         
    colormap(gca, colMapAlignment);
    caxis([-2 5]);
    ylim(imDis.*[1 3]+[-0.5 .5])
    box off
    axis off
    text(repmat(-5,3,1),[1,2,3]'.*imDis,[examplePairs,'\it{K.lac}'],'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',12)
        
    % scatters
    clear S
    for p = 1:numel(pattern)
        S{p} = regexp(SC_genomeCat, pattern{p});
    end
    occPattern = [];
    for i = 1:6701
        if ~isnan(promoterLengthsORF(i))
            intBasesGene = GP.chr_idx(GP.gene_infoR64.position(i,1)) + GP.gene_infoR64.position(i,2)+[-promoterLengthsORF(i):0]*(GP.gene_infoR64.dir(i));
            for p = 1:numel(pattern)
                occPattern(i,p) = sum(ismember(S{p}, intBasesGene));
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
    [~,maxIdx] = maxk(colDataS,kLacNtargets);
    
    % scatter WTs
    axes('Position', [xPos+Wscheme+xspacer+0.1*Wscatter yPos+0.1*Hscatter 0.9*Wscatter 0.9*Hscatter])
    scatter(xdataOriginal/max(xdataOriginal), ydataOriginal/max(ydataOriginal), 20,[0.7 0.7 0.7],'.','MarkerEdgeColor',[.3 .3 .3],'DisplayName','all genes')
    hold on
    scatter(xdataOriginal(maxIdx)/max(xdataOriginal), ydataOriginal(maxIdx)/max(ydataOriginal),20,colDataS(maxIdx)/max(colDataS),'filled','MarkerEdgeColor',[.3 .3 .3],'DisplayName','K.lactis top targets')
    colormap(gca, cMapLac)
    cbr = colorbar()
    ylabel(cbr, 'K.lac','fontSize',12)
    axisPos = gca;
    set(gca, 'YTick',[],'XTick',[])
    title(sprintf('r = %.2f',corr(xdataOriginal,ydataOriginal,'rows','pairwise')), 'fontSize',12, 'fontWeight','normal')
    caxis([0.2 0.6])
    set(cbr, 'Ticks', [0.2:0.2:0.6])
    ylim([0 1])
    xlim([0 1])
    
    % histogram yaxis
    axes('Position', [xPos+Wscheme+xspacer yPos+0.1*Hscatter 0.1*Wscatter 0.9*Hscatter])
    bins = [0:0.05:1];
    binM = movmean(bins,2,'Endpoints','discard');
    countsH = histcounts(ydataOriginal(maxIdx)/max(ydataOriginal), bins,'Normalization','pdf');
    barh(movmean(bins,2,'Endpoints','discard'),movmean(countsH,1), 'FaceColor', cMapLac(180,:),'LineStyle','none')
    hold on
    countsAll = histcounts(ydataOriginal/max(ydataOriginal), bins,'Normalization','pdf');
    fill([countsAll([1,1:end]),0,0], [0,binM([1:end, end]),0],[.7 .7 .7], 'LineStyle','none', 'FaceAlpha',0.75)
    set(gca, 'XDir','reverse','YLim', axisPos.YLim,'XTick',[], 'fontSize',10)
    ylabel(TF2,'fontSize',12)
    xlim([0 max(countsH)])
    currA = gca;
    set(gca, 'YTick', [0:0.5:1])
    
    % histogram xaxis
    axes('Position', [xPos+Wscheme+xspacer+0.1*Wscatter yPos axisPos.Position(3) 0.1*Hscatter])
    countsH = histcounts(xdataOriginal(maxIdx)/max(xdataOriginal), bins,'Normalization','pdf');
    bar(movmean(bins,2,'Endpoints','discard'), movmean(countsH,1), 'FaceColor', cMapLac(180,:),'LineStyle','none','DisplayName','K.lactis top targets')
    hold on
    countsAll = histcounts(xdataOriginal/max(xdataOriginal), bins,'Normalization','pdf');
    fill([0,binM([1:end, end]),0], [countsAll([1,1:end]),0,0],[.7 .7 .7], 'LineStyle','none', 'FaceAlpha',0.75,'DisplayName','all genes')
    set(gca, 'YDir','reverse','XLim', axisPos.XLim, 'YTick',[],'fontSize',10)
    xlabel(TF1,'fontSize',12)
    ylim([0 max(countsH)])
    set(gca, 'XTick', [0:0.5:1])
    legend('fontSize',10)
    
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
            fitGenes1=getFitGenes(xdataOriginal,3.5,0.99)
            normF1=robustfit(xdataOriginal(fitGenes1)/max(xdataOriginal),xdataNew(fitGenes1),[],[],'off');
            
            fitGenes2=getFitGenes(ydataOriginal,3.5,0.99)
            normF2=robustfit(ydataOriginal(fitGenes2)/max(ydataOriginal),ydataNew(fitGenes2),[],[],'off');
            
            %normF1 =mean(xdataNew(xTopIdx))/mean(xdataOriginal(xTopIdx)/max(xdataOriginal));
            %normF2 = mean(ydataNew(yTopIdx))/mean(ydataOriginal(yTopIdx)/max(ydataOriginal));
            titleStr = sprintf('%s vs %s (%.2f, %.2f)', TF1, TF2, normF1,normF2);
            
            xOriginalAdjusted = xdataOriginal*normF1/max(xdataOriginal);
            yOriginalAdjusted = ydataOriginal*normF2/max(ydataOriginal);
            
            %axes('position', [xPos+(Wscatter+xspacer)*ct yPos-Hscatter-yspacer Wscatter Hscatter])
            axes('Position', [xPos+Wscheme+xspacer+ct*(Wscatter+xspacer) yPos Wscatter Hscatter])
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
            setAxisExponent()
            set(gca,'fontSize',10)
            currA = gca;
            set(gca, 'XTIck', [0:0.5:1], 'YTIck', [0:0.5:1])
            cbr = colorbar()
            caxis([-0.5 3.5])
            set(cbr, 'Ticks', [0:3])
            title(cbr, '#motifs','fontSize',12)
            labelSamples = regexprep({comp1,comp2}, {'_d', '_lactis', '_DBD','_'}, {' Delta', ' K.lac', ' DBD',' '});
            xlabel(labelSamples{1},'fontSize',12)
            ylabel(labelSamples{2},'fontSize',12)
            legend([pUp,pDown], 'Location', 'southeast')
            %text(0.9*max(xlim), 0.05*max(ylim), sprintf('r = %.2f',corr(xdataNew,ydataNew,'rows','pairwise')), 'fontSize',12, 'HorizontalAlignment','center')
            title(sprintf('r = %.2f',corr(xdataNew,ydataNew,'rows','pairwise')), 'fontSize',12, 'fontWeight','normal')
        end
    end 
end
 
save_gf(gcf,'Fig5B','paper','tamar21','type','svg')




