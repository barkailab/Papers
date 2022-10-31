%% Figure 8—figure supplement 1
clearvars -except checWTdelLactisSwap
load('checWTdelLactisSwap.mat')
allDomain = readtable('./allDomainsAdj.txt');
intDomains = {'Zn_clus','zf-C2H2','Forkhead','KilA-N','bZIP_1','Vhr1','Homeodomain','Myb_DNA-binding','SRF-TF',...
    'AFT','GATA','Copper-fist','HSF_DNA-bind','zf-MIZ','HMG_box','TIG','zf-BED'};
allDBD = allDomain(ismember(allDomain.x_TargetName,intDomains),:);
KlacDBDs = allDBD(contains(allDBD.queryName,'Klac'),:);
ScerDBDs = allDBD(contains(allDBD.queryName,'Scer'),:);
ScerDBDs = ScerDBDs(ScerDBDs.score_1>0,:);
load('./promoterOL.mat','perOverlap','totalOL')
GP=load('group_imp.mat')
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
h3smooth=load('H3CC_henikoff.mat');
h3smooth=conv(mean(h3smooth.p_all,2),normpdf([-50:50]',0,25),'same');
clearvars -except summaryTable ScerDBDs KlacDBDs checWTdelLactisSwap perOverlap totalOL GP promoterIDXvec h3smooth promoterLengthsORF


%% Figure 8—figure supplement 1A,B: Hms2/Skn7, Ecm22/Upc2
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y','X', '-'};
B62 = blosum(62,'Order', cat(2,AAorder{1:end-1}));
B62(22,:) = -8;
B62(:,22) = -8;
colMapAlignment =brighten([1 1 1; flipud(bone)],0.6);
%DBDboxCol = [135 48 146]/255;
DBDboxCol = [0 230 230]/255;

CCpredictionsWJpred = load('./CCpredictionsWJpred.mat')
load('iupred2.mat')
load('./paraSeqs.mat');
%targets = {'YCR102C','MGA1','TPO1' ; 'ACS1','EEB1' ''; 'RPN4','','RSB1'; 'SNQ2','','SAM1'};
%targets = { 'YCR102C','TPO1','HXT6' ; 'ACS1','EEB1' 'INO1'; 'RPN4','RSB1','PDR3'; 'PLB1','AZR1','SAM1'}%; 'CCW12','YKL044W','PHD1'};

allSamples = fieldnames(checWTdelLactisSwap.sumProm);
zc = {'Skn7','Hms2';'Upc2','Ecm22'};
targets = {'NRG1','CWP2';'ERG25','ERG11'};
patternsColl{1} = {'GGCCA|TGGCC','GGCCG|CGGCC' };
patternsColl{2} = {'TATACGA|TCGTATA','TAAACGA|TCGTTTA' };

TargetsSortForFig5 = readtable('./TargetsSortForFig8.xlsx');
patternLoc = [-430, -543;-482,-849];
intoGene = 100;
promoterLengthPlot = 900;
promoterLengthPlotVec = repmat(promoterLengthPlot,6701,1);

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
save_gf(gcf,sprintf('IndRepeats_Matrix_Sup5'),'type',{'svg'},'paper','tamar21')

% colorbars:
axes('Units','normalized','Position',[xPos(1), yPos(2)-yspacer-0.01, Wsignal, 0.02])
caxis([0 10])
colormap(gca,cMapZscores)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'z-score', 'fontSize',10)
set(cbr, 'Ticks', [0.1+min(caxis),max(caxis)-1] , 'TickLabels', caxis)

axes('Position', [xPos(1)+xspacer+Wsignal, yPos(2)-yspacer-0.01, Wsignal, 0.02])
caxis([0 0.5])
colormap(gca,cMapMotifs2(1:size(cMapMotifs2,1)/2,:))
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'motif signal', 'fontSize',10)
set(cbr, 'Ticks', [min(caxis)+0.02,max(caxis)-0.02] , 'TickLabels', {'0', '50%'},'TickLength',[0 0])

axes('Position', [xPos(1)+2*(xspacer+Wsignal), yPos(2)-yspacer-0.01, Wsignal, 0.02])
caxis([-2 2])
colormap(gca,cMapLogChange)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'Deltalog2', 'fontSize',10)
set(cbr, 'Ticks', [0.5+min(caxis),0,max(caxis)-0.1] , 'TickLabels', [min(caxis),0,max(caxis)],'TickLength',[0 0])

% alignment colorbar
axes('Position', [0.05, yPos(4)+Hsignal-Hscheme/2-2.5*(Hscheme/2+yspacer/1.2), Wscheme,  Hscheme/2])
caxis([-4,4])
colormap(gca, colMapAlignment)
cbr = colorbar('Location', 'south')
cbr.AxisLocation = 'out'
axis off
ylabel(cbr,'aa conservation score', 'fontSize',12)
set(cbr,'fontSize',11)

%save_gf(gcf,'FigS5SU','paper','tamar21','type','svg','size',[])
