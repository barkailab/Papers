%% Figure 1-figure supplement 1
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
load('SGDTargetsRegulatorsResults.mat')
load('allSGDtargets.mat');
GP=load('./group_imp.mat')

%% Figure 1-figure supplement 1B: violin plots, SGD targets, sumProm scatters
familyNames = {'Zinc finger', 'Zinc cluster, bZIP', 'others'};
sumPromWTs = nan(6701, size(summaryTable,1)*2);
c=1;
for i = 1:size(summaryTable,1)
   sumPromWTs(:,c) = checWTdelLactisSwap.sumProm.(summaryTable.p1{i});
   sumPromWTs(:,c+1) = checWTdelLactisSwap.sumProm.(summaryTable.p2{i});
   c=c+2;
end
sumPromWTs(sumPromWTs == 0) = nan;
TFnames = [summaryTable.p1,summaryTable.p2]';

figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
W = 0.7;
H = 0.12;
yStart = 0.82;
xspacer = 2*W/19.8;
Wscatter = xspacer*2/3;
Hscatter = 1.8*Wscatter;
yspacer = Hscatter+4*H/5;

yPosScatter = repmat([yStart-yspacer+2*H/5: -(yspacer+H): 0],10,1);
xPosScatter = repmat([0.06: xspacer: 0.06+9*xspacer]',1,3);

ax(1) = axes('Position', [0.05 yStart W H]);
ax(2) = axes('Position', [0.05 yStart-yspacer-H W H]);
ax(3) = axes('Position', [0.05 yStart-2*(yspacer+H) W H]);
colM = cbrewer2('Paired',60);
scatterCol = [0.7 0.7 0.7];
pos = [1:2:30*2]+[0;0.8];

for i = 1:2*size(summaryTable,1)
    if mod(i-1,20) == 0
        axes(ax(ceil(i/20)));
    end
    violins(i) = Violin(log10(sumPromWTs(:,i)+700), pos(i));
    delete(violins(i).MedianPlot)
    violins(i).ViolinColor = colM(i,:);
    delete(violins(i).ScatterPlot)
    hold on
    bIdx = GP.gene_table.(upper(TFnames{i}));
    isTargetManual = ismember(allSGDtargets.locus1, GP.gene_infoR64.orf(bIdx))& contains(allSGDtargets.annotation_type, 'manually');
    isTargetChip = ismember(allSGDtargets.locus1, GP.gene_infoR64.orf(bIdx))& contains(allSGDtargets.evidence, 'chromatin immunoprecipitation-chip evidence');
    [~,targetIdxManual] = ismember(allSGDtargets.locus2(isTargetManual),GP.gene_infoR64.orf);
    [~,targetIdxChip] = ismember(allSGDtargets.locus2(isTargetChip),GP.gene_infoR64.orf);
    targetIdxManual = targetIdxManual(targetIdxManual>0);
    targetIdxChip = targetIdxChip(targetIdxChip>0);
    
    s1(i) = scatter(repmat(pos(i), numel(targetIdxManual),1), log10(sumPromWTs(targetIdxManual, i)+700), '.k','DisplayName', 'manually curated targets');
    s2(i) = scatter(pos(i), mean(log10(sumPromWTs(targetIdxChip, i)+700),'omitnan'), rescale(numel(targetIdxChip),20, 100, 'InputMin',5, 'InputMax', 50),...
        'r', 'LineWidth',2, 'DisplayName', 'ChIP-chip targets (mean)');
end

for i = 1:3
    axes(ax(i))
    set(ax(i), 'XTick', pos(i*20+[-19:0]), 'XTickLabel', TFnames(i*20+[-19:0]),'TickLength',[0 0], 'YLim', [log10(700),6],'Xlim',quantile(pos(i*20+[-19:0]),[0 1])+[-0.5 0.5], 'FontSize',12)
    text(sum([0.99 0.01].*xlim), 5.5, familyNames{i}, 'fontSize',15, 'HorizontalAlignment', 'left', 'VerticalAlignment','middle')
    if i==1
        legend([s1(1),s2(1)])
    end
    if i == 2
        ylabel('Promoter signal (log10)', 'fontSize',15)
    end
end

% scatters
for i = 1:size(summaryTable,1)
    axes('Position', [xPosScatter(i), yPosScatter(i), Wscatter,Hscatter]);
    scatter(sumPromWTs(:,i*2-1), sumPromWTs(:,i*2), 15, scatterCol, 'filled')
    setAxisExponent()
    text(0, max(ylim), sprintf(' %.2f', corr(sumPromWTs(:,i*2-1), sumPromWTs(:,i*2),'rows','pairwise')), 'HorizontalAlignment','left','VerticalAlignment','top')
end
%saveas(gcf, 'FigS1violinAndScatters.svg')


%% Figure 1-figure supplement 1C: binding signal corr matrix including repeats (promoters, 7mers)
goodPos = createIntBasesForMotif();
TFnames = [summaryTable.p1,summaryTable.p2]';
[~, idxVec, sumPromRep, repIdx] = getRepeatsCorr(TFnames(:),'dataType','sumProm');
[~, ~, mer7Rep, ~] = getRepeatsCorr(TFnames(:),'dataType','7mer');
samplesWOrepeats = unique(idxVec(all(isnan(sumPromRep),1)));
for i = samplesWOrepeats'
    sumPromRep(:,idxVec==i) = repmat(checWTdelLactisSwap.sumProm.(TFnames{i}),1,sum(idxVec==i));
    currNorm = chromosomes2fullProfile(checWTdelLactisSwap, TFnames(i));
    currMer = mer_occupancy(currNorm,7,'intBases', goodPos,'method','normal');
    mer7Rep(:,idxVec==i) = repmat(currMer.score,1,sum(idxVec==i));
end

sumPromCorrMat = corr(sumPromRep,'rows','pairwise');
mer7CorrMat = corr(mer7Rep,'rows','pairwise');
combineMat = tril(sumPromCorrMat) + triu(mer7CorrMat,1);
colMap = brighten(flipud(bone),0.2);

figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
imagesc(combineMat)
borders = [0; cumsum(accumarray(idxVec,1))]+0.5;
tickPos = movmean(borders,2,'Endpoints','discard');
nRepeats = accumarray(idxVec, repIdx, [], @(x)numel(unique(x(~isnan(x)))));
yLabels = strcat(TFnames(:), ' (',num2str(nRepeats), ')')
yLabels = strrep(yLabels,'(0)','(*)');

set(gca, 'YTick', tickPos, 'YTicklabel', yLabels,'XTick',tickPos,'XTickLabel',TFnames,'XTickLabelRotation',90, 'TickLength',[0 0])
hold on
plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)
colormap(gca, colMap)
caxis([0.2 1])
axis square

cbr = colorbar()
ylabel(cbr, 'correlation','fontSize',20)

