%% Figures 1 and 2
%%
load('summaryTable.mat')
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
h3cc=load('H3CC_henikoff.mat');
h3smooth=conv(mean(h3cc.p_all,2),normpdf([-50:50]',0,25),'same');
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
load('SC_genome.mat')
load('allZscoreMat.mat', 'allSamples', 'allZscoreMat','sumPromMat')
load('promCorrSort.mat')
load('checWTdelLactisSwap.mat')
GP=load('group_imp.mat');

%% Figure 2 - histograms
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
groups = fieldnames(promCorrSort);
examplePairs = {'Met31', 'Met32'; 'Gzf3','Dal80';'Ace2','Swi5'};
familyNames = {'Zinc finger', 'Zinc cluster, bZIP', 'others'};
groupsIdx = {1, [2,3], 4};
Blues = cbrewer2('Blues');
Bugn = cbrewer2('BuGn');
oranges = cbrewer2('Oranges');
purples = cbrewer2('Purples');
%colorsP = [77 148 255; 204 153 255; 255 173 51 ; 64 191 128]/255;
colorsP = [76 0 75; 130 18 126; 147 174 210 ; 138 139 192]/255;

familyNamesLabels = {sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Zinc finger', colorsP(1,:)),...
    sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Zinc cluster, \\color[rgb]{%.2f,%.2f,%.2f}bZIP', colorsP(2,:),colorsP(3,:)),...
    sprintf('\\color[rgb]{%.2f,%.2f,%.2f}Others', colorsP(4,:))};

% calculating the corr between all factors
promCorrMatall = plotSumPromCorr(geneName, checWTdelLactisSwap, 0); 
c = 1;
clear promCorrVecAll
for i = 1: length(geneName)
    for z = 1: length(geneName)
        if z > i
        promCorrVecAll(c) = promCorrMatall(i,z);
        c = c+1;
        end
    end
end
[counts, bins] = histcounts(promCorrVecAll ,20);

figure('Units','normalized','Position', [0 0 0.5 1],'color','w')
xPos = 0.1;
yPos = [0.7 0.58 0.46];
Hhist = 0.1;
Whist = 0.2;
for g =  1: length(groupsIdx)
    nPara =  sum(ismember(summaryTable.familyId, groupsIdx{g}));
    
    % plot histogram
    axes('Position',[xPos yPos(g) Whist Hhist])
    area(movmean(bins,2,'Endpoints','discard'), movmean(rescale(counts, 0,3),3), 'FaceColor',[0.7 0.7 0.7], 'FaceAlpha',0.5, 'LineStyle', 'none');
    hold on
    for f = groupsIdx{g}
        selTF = find(summaryTable.familyId == f);
        plot(summaryTable.sumPromCorr(selTF)'.*[1;1], repmat([0;1], 1, numel(selTF)),'color',colorsP(f,:), 'LineWidth',3)
    end
    ylim([0 3.5])
    set(gca, 'fontSize',15,'TickLength',[0 0])
    xlim([-0.0735 1])
    if g==3
        xlabel('Promoter signal correlation', 'FontSize',18);
        set(gca, 'TickLength',[0 0])
    else
        set(gca, 'XTick',[]);
    end
    if g==2
        ylabel('Number of pairs','FontSize',18);
    end
    text(sum(xlim.*[0.95 0.05]), max(ylim)*0.98, familyNamesLabels{g}, 'fontSize',15,'fontWeight','bold', 'VerticalAlignment','top')
    
    if any(contains(summaryTable.p1(selTF), examplePairs))
        selR = find(contains(summaryTable.p1(selTF), examplePairs));
        for r = selR'
            examplePairsIdxCorr = summaryTable.sumPromCorr(selTF(r));
            currExamplePairs = [summaryTable.p1{selTF(r)},...
                '\newline',summaryTable.p2{selTF(r)}];
            text(examplePairsIdxCorr, 1.9, currExamplePairs,'VerticalAlignment','bottom','HorizontalAlignment','center','fontSize',13)
            plot([examplePairsIdxCorr,examplePairsIdxCorr], [1.25,1.85],'k','LineWidth',1)
            scatter(examplePairsIdxCorr, 1.25, 25,'kv','filled')
        end
    end
end

%% Figure 2 - Zscores and barplot (bar colored by conservation group: strong, mid, weak)
% dividing into groups - ZF, zinc clusters+bZIP, other families, special cases

% load('corrBetweenRepeatsWTs.mat')
% for i = 1:size(summaryTable,1)
%     if isfield(corrBetweenRepeatsWTs, summaryTable.p1{i})
%     medianCorrRepeats(i) = median([corrBetweenRepeatsWTs.(summaryTable.p1{i}),corrBetweenRepeatsWTs.(summaryTable.p2{i})]);
%     else
%     end
% end
groups = fieldnames(promCorrSort);
familyNames = {'Zinc finger', 'Zinc cluster, bZIP', 'others'};
groupsIdx = {1, [2,3], 4};
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);

% calculating the corr between all factors
promCorrMatall = plotSumPromCorr(geneName, checWTdelLactisSwap, 0); 
c = 1;
clear promCorrVecAll
promCorrVecAll = 1-squareform(1-promCorrMatall);
xDis = 0.32;
xPos = [0.025:xDis:1];
Wbars = 0.05;
Wzscore = 0.2;
Hzscore = Wzscore*2;
yspacer = 0.06;
xspacer = 0.04;
Hhist = 0.1;
yPos = 0.8;
cMapZscore = brighten(flipud(bone),0.2);
%cMapBars = cbrewer2('BuPu',100);
cMapBars = [126,8,119; 138,139,192; 183,205,227]/255;

% groups classification by conservation
WTsCorr = summaryTable.sumPromCorr;
r = [0.81,0.4];
consGroupsLog = [WTsCorr>=r(1), WTsCorr<r(1)&WTsCorr>=r(2), WTsCorr<r(2)]; 

[counts, bins] = histcounts(promCorrVecAll ,20);
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
for g =  1: length(groupsIdx)
    nPara =  sum(ismember(summaryTable.familyId, groupsIdx{g}));   
    % plot bars
    axes('position', [xPos(g) yPos-Hzscore-yspacer Wbars Hzscore])
    c=1;
     for p = find(ismember(summaryTable.familyId, groupsIdx{g}))'
        hold on
        b = barh(c,summaryTable.sumPromCorr(p), 'FaceColor', cMapBars(consGroupsLog(p,:),:));
        c=c+1;
    end
    set(gca, 'YDir', 'reverse')
    set(gca, 'XDir', 'reverse')
    ylim([0, nPara]+0.5)    
    line([0.50 0.50], ylim, 'Color','k','LineStyle','--')
    box on;
    set(gca, 'ytick', [1:nPara], 'yticklabels',[], 'FontSize',12, 'XTick', [0:0.5:1], 'XAxisLocation', 'bottom');
    title('Correlation', 'FontSize',13,'FontWeight','bold');
    
    axes('position', [xPos(g) yPos-Hzscore-1.7*yspacer Wbars yspacer/4])
    imagesc([1:3])
    plotgrid([1:3])
    colormap(gca, cMapBars)
    set(gca, 'YTick',[],'XTick', [0.7,2,3.3], 'XTickLabel', {'strong','mid','weak'},...
        'TickLength',[0 0], 'fontSize',13,'XTickLabelRotation',20)

    
    % plot Zscores matrix
    axes('position', [xPos(g)+Wbars+xspacer yPos-Hzscore-yspacer Wzscore Hzscore])
    N = size(summaryTable.nHighTargetsIdx{1}, 1);
    c=1;
    selRows = find(ismember(summaryTable.familyId, groupsIdx{g}));
    currTFs = reshape([summaryTable.p1(selRows)';summaryTable.p2(selRows)'], 1,[]);
    allHighTargetsIdx = cat(1,summaryTable.nHighTargetsIdx{selRows});
    zscoreMat = nan(nPara*2, N*nPara);
    for p = selRows'
        nanVec1 = nanZscore(checWTdelLactisSwap.sumProm.(summaryTable.p1{p}));
        nanVec2 = nanZscore(checWTdelLactisSwap.sumProm.(summaryTable.p2{p}));
        zscoreMat(c,:) = nanVec1(allHighTargetsIdx)';
        zscoreMat(c+1,:) = nanVec2(allHighTargetsIdx)';
        c=c+2;
    end
    imagesc(zscoreMat)
    axis square
    caxis([0 10])
    colormap(gca, cMapZscore)
    hold on
    for z = 1: 2: size(zscoreMat,2)
        plot(xlim, [z+1+0.5 z+1+0.5], 'k', 'LineWidth', 2);
    end
    for z = 1:size(zscoreMat,2)/(N/2)
        plot([N*z N*z], [ylim], 'k', 'LineWidth', 2);
    end
    set(gca, 'YTick', [1:nPara*2], 'YTickLabel', currTFs ,'FontSize',12,'XTick',[], 'XAxisLocation','top');
    if g == 3
        cbr = colorbar('Location', 'east')
        set(cbr, 'position', [xPos(g)+xspacer/1.3+Wbars+Wzscore+0.015 yPos-Hzscore-yspacer 0.01 Hzscore]);
        set(cbr, 'AxisLocation', 'out','Ticks',[0:2:10],'FontSize',13)
        ylabel(cbr, 'Z score', 'fontSize',15,'fontweight','bold')
    end
    title(familyNames{g},'fontSize',18)
    
    axes('position', [xPos(g)+Wbars+xspacer yPos-Hzscore-yspacer-0.03 Wzscore Hzscore/20])
    plot(xlim/nPara-0.02, [0.5 0.5], 'color','k','LineWidth',1.5)
    xlim([0 1])
    hold on
    plot([0.102, (max(xlim)/nPara*2-0.02)], [0.5 0.5], 'color','k','LineWidth',1.5)
    plot([0.202, (max(xlim)/nPara*3-0.02)], [0.5 0.5], 'color','k','LineWidth',1.5)
    
    scatter((max(xlim)/nPara-0.02)+[0,0.1,0.2], repmat(0.5,3,1), 60, 'k>','filled')
    scatter((max(xlim)/nPara-0.02)+0.22+[0.025, 0.05, 0.075], repmat(0.5,3,1), 15, 'k', 'filled')
    text([0.1 0.2 0.3], repmat(2,3,1), {'40','40','40'}, 'HorizontalAlignment', 'center','VerticalAlignment','middle','fontSize',10)
    text(0, 2, '1', 'HorizontalAlignment', 'center','VerticalAlignment','middle','fontSize',10)

    axis off
    text(0.4, 0.5, 'Top targets (ordered)',...
        'HorizontalAlignment','left','VerticalAlignment','middle', 'fontSize',13)
end


%% Figure 1D - signal on promoter examples and sumProm scatter
colMotifs = lines(6);
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
examplePairs = {'Met31', 'Met32'; 'Gzf3','Dal80'; 'Ace2','Swi5'};
targets = {'MUP1', 'MET2'; 'GAT1','FUR4';'DSE1','PIR1'};
motifList = {'GCCAC', 'GATAA', 'GCTGG'};
intoGene = 150;
promoterL = repmat(700,6701,1);
lineCol = lines(2);
col = cbrewer2('Dark2',1);
scatterCol = [51 51 51]/255;%[148 175 211]/255;
markerCol = brighten(scatterCol, -0.5);
colSignal = [0 0 0]; 

% changing TSS position of GAT1
geneIdx=GP.gene_table.GAT1;
GP.gene_infoR64.tamarTss(geneIdx,2) = GP.gene_infoR64.stein_tss(geneIdx,2)

xPos = [0.1, 0.27, 0.44];
Hsignal = 0.04;
Wsignal = 0.12;
Hscatter =Wsignal*1.8;
yspacer = 0.03;
yPos = [0.9 0.9 0.9];

figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
for i = 1:size(examplePairs,1)
    % signal on example promoters
    for g = 1: size(targets,2)
        target = targets{i,g};
        geneIdx = GP.gene_table.(target);
        pos = GP.gene_infoR64.position(geneIdx,:);
        intBases = GP.chr_idx(pos(1))+pos(2)+[-promoterL(geneIdx):intoGene].*GP.gene_infoR64.dir(geneIdx);
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
         
         axes('Position', [xPos(i) yPos(i)-(g-1)*(Hsignal+yspacer) Wsignal Hsignal])
         for z = 1:2
             currNorm = chromosomes2fullProfile(checWTdelLactisSwap, examplePairs(i,z));
             maxNorm = max(currNorm(logical(promoterIDXvec)));
             currProfile = currNorm(intBases)./maxNorm;
             plot([-promoterL(geneIdx):intoGene], rescale(currProfile, 3-z+0.5, 3-z+1.3,'InputMin',0, 'InputMax',0.5),...
                 'color', colSignal)
             hold on
         end
         plot(-promoterL(geneIdx):intoGene, rescale(h3smooth(intBases), 0.5, 1.3),'k', 'LineWidth',1)
         area(-promoterL(geneIdx):intoGene, rescale(h3smooth(intBases), 0.5, 1.3),'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'LineStyle', 'none')
         ylim([0.5 2+1.5])
         xlim([-promoterL(geneIdx) intoGene])
         TssLoc = (GP.gene_infoR64.tamarTss(geneIdx,2)-pos(2))*GP.gene_infoR64.dir(geneIdx);
         set(gca, 'YTick', [1:3], 'YTickLabel', {'NucOcc', examplePairs{i,2},examplePairs{i,1}},...
             'TickLength',[0 0], 'XTick',[], 'fontSize',10)
         plot(TssLoc*[1 1] ,ylim, ':r', 'LineWidth',1.5)
         scatter(0 ,0.5,[], '>r','filled')
         if numel(motifPos.S) > 0
             for p = 1:numel(motifPos.S)
                 fill([motifPos.S(p)*[1,1], motifPos.E(p)*[1,1]], [0.5, 3+0.5.*[1,1], 0.5], colMotifs(1,:),...
                     'LineStyle','none', 'FaceAlpha', 0.5)
             end
         end
         text(TssLoc, 0.5, 'TSS','color','r','fontSize',8,'HorizontalAlignment','center','VerticalAlignment','top')

            title(sprintf('%s promoter', target),'fontSize',10)
            box on
            if g==2
                text(min(xlim), 0.5, ['{\itin-vitro} motif: ', motifList{i}],'Color', colMotifs(1,:),'VerticalAlignment','top','fontSize',10)
            end
    end
  
    %scatter
    axes('Position', [xPos(i) yPos(i)-Hsignal-2*yspacer-Hscatter-0.01 Wsignal Hscatter])
    maxX = max(checWTdelLactisSwap.sumProm.(examplePairs{i,1}));
    maxY = max(checWTdelLactisSwap.sumProm.(examplePairs{i,2}));
    scatter(checWTdelLactisSwap.sumProm.(examplePairs{i,1})/maxX,...
        checWTdelLactisSwap.sumProm.(examplePairs{i,2})/maxY, [],...
        scatterCol, 'filled','MarkerEdgeColor', markerCol)
    [~,intGenes] = ismember( targets(i,:), GP.gene_infoR64.nameNew);
    hold on
    scatter(checWTdelLactisSwap.sumProm.(examplePairs{i,1})(intGenes)/maxX,...
        checWTdelLactisSwap.sumProm.(examplePairs{i,2})(intGenes)/maxY, [], [140,152,199]/255, 'filled','MarkerEdgeColor', 'k')
    text(checWTdelLactisSwap.sumProm.(examplePairs{i,1})(intGenes)/maxX+0.02*max(xlim),...
        checWTdelLactisSwap.sumProm.(examplePairs{i,2})(intGenes)/maxY, targets(i,:), 'fontSize',12, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment','middle')
    xlabel(sprintf('%s', examplePairs{i,1}), 'fontSize',13)
    ylabel(sprintf('%s', examplePairs{i,2}), 'fontSize',13)
    set(gca, 'fontSize',12)
    title('   Promoter binding','FontWeight','normal', 'fontSize',13)
    set(gca,'XTick',[0:0.5:1], 'YTick',[0:0.5:1])
    setAxisExponent()
    yLim = 1.05;
    text(max(xlim)*0.05,max(yLim)*0.95,...
        sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(examplePairs{i,1}), checWTdelLactisSwap.sumProm.(examplePairs{i,2}), 'rows','pairwise')),...
        'fontSize',13, 'HorizontalAlignment', 'left')
    ylim([0 yLim])
    xlim([0 yLim])
end

%% Figure 1E - Auto and cross promoter binding 
% calculating parameter
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
WTs = geneName;
[~,WTsIdx] = ismember(geneName, allSamples);
[~,bIdx] = ismember(upper(WTs), GP.gene_infoR64.nameNew);
TFmat = allZscoreMat(bIdx, WTsIdx);
WTallMat = allZscoreMat(:,WTsIdx);
% defining zScore
nTargets = 50;
zTH1 = 3.5;
zTH = min(quantile(WTallMat, 0.99), zTH1);

% create a target matrix
allTargetsMat = WTallMat>=zTH;
TFtargetsMat = TFmat>=zTH;

clear newOrder motifType paraCorr TFlabel nTargets iTargets nInputs iInputs
for tf = 1:2:numel(WTs)
    if sum(TFtargetsMat([tf,tf+1], tf).*[1.5;1]) < sum(TFtargetsMat([tf+1,tf], tf+1).*[1.5;1])
        newOrder{(tf+1)/2} = [tf+1,tf];
    else
        newOrder{(tf+1)/2} = [tf,tf+1];
    end
end

newOrder = cat(2,newOrder{:});
allTargetsMatOrdered = allTargetsMat(:,newOrder);
TFtargetsMatOrdered = TFtargetsMat(newOrder,newOrder);
WTsOrdered = WTs(newOrder);
WTallMatOrdered = WTallMat(:,newOrder);
TFmatOrdered = TFmat(newOrder, newOrder);

for tf = 1:2:numel(WTsOrdered)
    motifType((tf+1)/2) = sum(TFtargetsMatOrdered(sub2ind(size(TFtargetsMatOrdered), [tf,tf+1,tf+1,tf],[tf,tf,tf+1,tf+1])).*[1 2 4 8])+1 ;
    imageMatB((tf+1)/2,:) =  TFtargetsMatOrdered(sub2ind(size(TFtargetsMatOrdered), [tf,tf+1,tf+1,tf],[tf,tf,tf+1,tf+1]));
    imageMatZ((tf+1)/2,:) = TFmatOrdered(sub2ind(size(TFmatOrdered), [tf,tf+1,tf+1,tf],[tf,tf,tf+1,tf+1]));
    paraCorr((tf+1)/2) = corr(WTallMatOrdered(:,tf), WTallMatOrdered(:,tf+1), 'rows','pairwise');
    TFlabel((tf+1)/2) = strcat(WTsOrdered(tf),',',WTsOrdered(tf+1));
    nTargets((tf+1)/2) = sum(any(allTargetsMatOrdered(:, [tf,tf+1]), 2));
    iTargets((tf+1)/2,:) = sum(allTargetsMatOrdered(:, [tf,tf+1]), 1);
    nInputs((tf+1)/2) = sum(any(TFtargetsMatOrdered([tf,tf+1],:),1));
    iInputs((tf+1)/2,:) = sum(TFtargetsMatOrdered([tf,tf+1],:),2)';
end

%% plot circuits
[matifsOrdered,imageOrder] = sort(motifType);
WTsCorr = summaryTable.sumPromCorr;
r = [0.81,0.4];
consGroupsLog = [WTsCorr>=r(1), WTsCorr<r(1)&WTsCorr>=r(2), WTsCorr<r(2)]; 
cMap = cbrewer2('BuPu') ;
cMapBars = [126,8,119; 138,139,192; 183,205,227]/255;

figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
axes('Position', [0.1 0.5 0.4 0.15])
imagesc(imageMatZ(imageOrder(matifsOrdered~=1),:)')
plotgrid(imageMatZ(imageOrder(matifsOrdered~=1),:)', 'color',[0.7 0.7 0.7])
caxis([0 10])
colormap(gca,cMap)
xTickLabels = TFlabel(imageOrder(matifsOrdered~=1));
for i = 1:numel(xTickLabels)
    GIdx(i) = 4-find(consGroupsLog(contains(summaryTable.p1,strsplit(xTickLabels{i},',')),:));
end
set(gca, 'XTick', [],'TickLength', [0 0 ],...
    'YTick', [1:4], 'YTickLabel', {['P1',char(8594),'P1'],['P1',char(8594),'P2'],['P2',char(8594),'P2'], ['P2',char(8594),'P1']}, 'fontSize',14)
cbr = colorbar()
ylabel(cbr, 'Z score', 'fontSize',12)
hold on
linePos = find(motifType(imageOrder(1:end-1)) ~= motifType(imageOrder(2:end)));
linePos = linePos(2:end) - linePos(1);
%plot(linePos+0.5*[1;1], repmat(ylim',1, numel(linePos)),'color','k', 'LineWidth',2)
currAxis = gca;

xTickLabelsCol = strcat('\color[rgb]{0,0,0}', xTickLabels);
axes('Position', currAxis.Position+[0 -0.02 0 -0.9*currAxis.Position(4)])
scatter([1:numel(xTickLabels)], repmat(0.5,1,numel(xTickLabels)), rescale(GIdx, 15,80), 'k', 'filled')
xlim([0.5, numel(xTickLabels)+0.5])
set(gca, 'XTick', [1:sum(matifsOrdered~=1,2)], 'XTickLabel', xTickLabelsCol,'TickLength', [0 0 ], 'XAxisLocation','bottom','XTickLabelRotation', 45,...
    'YTick', 0.5, 'YTickLabel', '\color[rgb]{0,0,0}Conservation', 'fontSize',12)
set(gcf, 'color','w')
set(gca, 'YColor',[1 1 1], 'XColor',[1 1 1])
