%% Figures 6 and 7
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
load('SC_genome.mat')
load('Abf2Ixr1.mat')
load('./treeFiles.mat');
load('./paraSeqs.mat');

%% Figure 6D - bar plot of paralog's distance based on nonDBD trees (values normalized to the distance of lactis)
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
geneName = [{'Ixr1'},{'Abf2'},geneName];
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
lacSamples = allSamples(contains(allSamples,'lactis') & ~contains(allSamples,'Sut1'));
YL=[{'Ixr1','Abf2'};summaryTable{[1:10,21:30,17:20,11:16],{'p1','p2'}}];

treeFiles = treeFiles(contains(treeFiles.name, '_non')&contains(treeFiles.name, geneName, 'IgnoreCase', true),:);
% 3-bars graph
[~,treeIdx]=ismember(upper(YL),regexp(treeFiles.name,'.{3}\d+','match','once'))
treeIdx=sum(treeIdx,2);
[~,sumidx]=ismember(YL(:,1),summaryTable.p1)

for i = 1: size(YL,1)    
    try
    [crMat, idxVec, ~, repIdx] = getRepeatsCorr(YL(i,[1,2]),'dataType','sumProm');
    [repIdx,uRpt]=unique(repIdx,'stable');
    idxVec=idxVec(uRpt);
    crMat=crMat(uRpt,uRpt);
    meanCrProm(i,1)=mean(crMat(idxVec==1,idxVec==2),'all');
    stdProm(i,1)=std(crMat(idxVec==1,idxVec==2),[],'all');
    catch
        try
            meanCrProm(i,1)=summaryTable.sumPromCorr(contains(summaryTable.p1,YL(i,:)));
        catch
            meanCrProm(i,1)=NaN;
        end                
        stdProm(i,1)=NaN;
    end
end

%resort treeFiles
paralogsDis = abs(treeFiles.dis1(treeIdx) - treeFiles.dis2(treeIdx));

familyNames = {'Zinc finger', 'others','bZIP', 'Zinc cluster'};
familyDiv=[1,11.5,22,26.7]+0.75;
axes('Position', [0.2 0.2 0.1 0.7])
barCol = [cbrewer2('Pastel1',2) ;[1,1,1]];
barCol(1,:) = [1  1 1]*0.7;
barPos = [1,1.5+[1:10, 11.5:20.5, 22:25, 26.5:31.5]];

% vertical bars
close all
figure('Units','normalized','Position',[0.5698 0.2324 0.3344 0.5019], 'color','w')

subplot(2,1,1)
bar(barPos-0.1, paralogsDis','FaceColor',barCol(2,:))
set(gca,'Ylim',[0 5],'Xlim',quantile(barPos,[0 1])+[-.5 .5])
hold on
plot( [1;1].*familyDiv,repmat(ylim',1,4), 'color',[1 1 1]*0.3, 'LineStyle', '--','LineWidt',1.5)
text(barPos([1,11,21,25]+1),repmat(0.99.*max(ylim),4,1), familyNames, 'fontSize',12, 'HorizontalAlignment','left','VerticalAlignment','top')
ylabel('D2')
set(gca, 'XTick', barPos, 'XTicklabel',strcat(YL(:,1),', ',YL(:,2)), 'fontSize',10,'XAxisLocation','top','XTickLabelRotation',90)

subplot(6,1,4)
bar(barPos,treeFiles.dis1(treeIdx),'FaceColor',barCol(1,:))
set(gca,'Ylim',[0 5],'Xlim',quantile(barPos,[0 1])+[-.5 .5])
hold on
plot( [1;1].*familyDiv,repmat(ylim',1,numel(familyDiv)), 'color',[1 1 1]*0.3, 'LineStyle', '--','LineWidt',1.5)
set(gca, 'XTick', [], 'fontSize',10)
ylabel('D1')
subplot(9,1,[7]+1)
crType=sum(meanCrProm>=[0,.4,.81],2)
cMapBars = flipud([126,8,119; 138,139,192; 183,205,227]/255);
hold off
for i=1:max(crType)
    bar(barPos,meanCrProm.*(crType==i),'BarWidth',.8,'FaceColor',cMapBars(i,:),'EdgeColor',[0 0 0])
    hold on
end
set(gca,'Ylim',[0 1],'Xlim',quantile(barPos,[0 1])+[-.5 .5])
plot( [1;1].*familyDiv,repmat(ylim',1,4), 'color',[1 1 1]*0.3, 'LineStyle', '--','LineWidt',1.5)
set(gca, 'XTick', [], 'TickLength',[0 0], 'fontSize',10,'XAxisLocation','bottom')
ylabel('binding correlation')
errorbar(barPos,meanCrProm,stdProm,'vertical','LineStyle','none','Color','k')

%% Ixr1 and Abf2: signal on example chr + mito + Ixr1 scatter with lactis
chr_idx = [0 cumsum(GP.chr_len)];
chr_id = sum([1:chr_idx(18)]'> chr_idx, 2); 
chrExample = 13;
chrMito = 17;
IxrCol = [153 102 153]/255;
AbfCol = [71 71 107]/255;
cMap = flipud(bone);
colSpot = [[140 113 180]/255;[148 175 211]/255];

figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
% signal on chr +mito
ax = axes('Position', [0.03 0.61 0.05 0.12*2+0.08])
set(ax,'xcolor','none', 'ycolor','none')
ylabel('Normalized signal','fontSize',20, 'Color','k')
ax.YAxis.Label.Visible='on';

axes('Position', [0.05 0.8 0.19 0.08])
plot(-Abf2Ixr1.meanNormProfiles(chr_id == chrExample,1),'Color',colSpot(2,:), 'LineWidth', 1.5);
hold on
plot(Abf2Ixr1.meanNormProfiles(chr_id == chrExample,2),'Color',colSpot(1,:), 'LineWidth', 1.5);
axis tight
ylim([-max(max(Abf2Ixr1.meanNormProfiles(chr_id == chrExample,2),Abf2Ixr1.meanNormProfiles(chr_id ==chrExample,1))), max(max(Abf2Ixr1.meanNormProfiles(chr_id == chrExample,2),Abf2Ixr1.meanNormProfiles(chr_id ==chrExample,1)))])
set(gca,'XTick', [1,GP.chr_len(chrExample)], 'XTicklabel', {'1', sprintf('%.0f kb',GP.chr_len(chrExample)/1000)} ,...
    'YTick', [-2000,0,2000], 'YTicklabel', abs([-2000,0,2000]), 'FontSize',12)
title('Chromosome XIII', 'fontSize',18)

ax = axes('Position', [0.25 0.8 0.05 0.08/2])
text(0.05, 0.5, 'Abf2','VerticalAlignment', 'middle', 'fontSize',13, 'Color',colSpot(2,:),'FontWeight','bold', 'Rotation',0, 'HorizontalAlignment','center')
axis off
ax = axes('Position', [0.25 0.8+0.08/2 0.05 0.08/2])
text(0.05, 0.5, 'Ixr1','VerticalAlignment', 'middle', 'fontSize',13, 'Color',colSpot(1,:),'FontWeight','bold', 'Rotation',0, 'HorizontalAlignment','center')
axis off

axes('Position', [0.05 0.65 0.19 0.08])
plot(-Abf2Ixr1.meanNormProfiles(chr_id == chrMito,1),'Color',colSpot(2,:), 'LineWidth', 1.5);
hold on
plot(Abf2Ixr1.meanNormProfiles(chr_id == chrMito,2),'Color',colSpot(1,:), 'LineWidth', 1.5);
axis tight
ylim([-max(max(Abf2Ixr1.meanNormProfiles(chr_id == chrMito,2),Abf2Ixr1.meanNormProfiles(chr_id == chrMito,1))), max(max(Abf2Ixr1.meanNormProfiles(chr_id == chrMito,2),Abf2Ixr1.meanNormProfiles(chr_id == chrMito, 1)))])
set(gca,'XTick', [1,GP.chr_len(chrMito)], 'XTicklabel', {'0', sprintf('%.0f kb',GP.chr_len(chrMito)/1000)},'YTick', [-10^6,0,10^6], 'YTicklabel', abs([-10^6,0,10^6]), 'FontSize',12)
title('mitochondrial genome', 'fontSize',18)
%legend({'Ixr1','Abf2'}, 'FontSize',12, 'Location', 'southwestoutside')
set(gcf,'Color','w')

ax = axes('Position', [0.25 0.65 0.05 0.08/2])
text(0.05, 0.5, 'Abf2','VerticalAlignment', 'middle', 'fontSize',13, 'Color',colSpot(2,:),'FontWeight','bold', 'Rotation',0, 'HorizontalAlignment','center')
axis off
ax = axes('Position', [0.25 0.65+0.08/2 0.05 0.08/2])
text(0.05, 0.5, 'Ixr1','VerticalAlignment', 'middle', 'fontSize',13, 'Color',colSpot(1,:),'FontWeight','bold', 'Rotation',0, 'HorizontalAlignment','center')
axis off

% Figure 7A - Ixr1 scatter with K. lactis
axes('Position', [0.05 0.21 0.2 0.2*1.8])
scatter(checWTdelLactisSwap.sumProm.Ixr1, checWTdelLactisSwap.sumProm.Ixr1_lactis, [],...
    checWTdelLactisSwap.sumProm.Ixr1_dABF2,'filled','MarkerEdgeColor','k')
linFit = linortfit2(checWTdelLactisSwap.sumProm.Ixr1(~isnan(checWTdelLactisSwap.sumProm.Ixr1)),...
    checWTdelLactisSwap.sumProm.Ixr1_lactis(~isnan(checWTdelLactisSwap.sumProm.Ixr1_lactis)));
xlabel('Ixr1', 'fontweight','bold', 'FontSize', 20);
ylabel('\itK.lac', 'fontweight','bold', 'FontSize',20);
set(gcf,'color','w');
hold on
plot(xlim, xlim*linFit(1),':')
caxis([0 0.5*max(caxis)])
cbr = colorbar()
ylabel(cbr, 'Ixr1 \DeltaABF2', 'fontSize',20, 'fontweight','bold')
colormap(gca, cMap)
set(cbr, 'fontSize',15)
set(gca, 'fontSize',15)
text(0.02*max(xlim),0.96*max(ylim),...
    sprintf('r=%.2f',corr(checWTdelLactisSwap.sumProm.Ixr1, checWTdelLactisSwap.sumProm.Ixr1_lactis,'rows','pairwise')),...
    'FontSize',15, 'HorizontalAlignment','left');

%% all protein alignment with lactis and z.rouxii GAPS REMOVED
% scoring lactis sequnce by multiple alignment sequence with ancestors 
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y','X', '-'};
B62 = blosum(62,'Order', cat(2,AAorder{1:end-1}));
B62(22,:) = -8;
B62(:,22) = -8;

ancestorMSA = dir('./mCoffee/*_fl.clustal*');
for a = 1:numel(ancestorMSA)
    msa = multialignread([ancestorMSA(a).folder,'/',ancestorMSA(a).name]);
    msa = msa(contains({msa.Header}, 'x'));
    [conSeq, ~] = seqconsensus(msa,'Gaps','all');
    if any(contains({msa.Header}, 'Klac'))
        lacSeq = msa(contains({msa.Header}, 'Klac')).Sequence;
        [~, idxLac] = ismember(lacSeq', AAorder);
        [~, idxCon] = ismember(conSeq', AAorder);
        score = B62(sub2ind(size(B62), idxLac,idxCon));
        gaps = strfind(lacSeq, '-');
        score(gaps) = [];
        lactisMSAscore.(extractBefore(ancestorMSA(a).name,'_')) = score';
    end
end

for a = 1:numel(ancestorMSA)
    msa = multialignread([ancestorMSA(a).folder,'/',ancestorMSA(a).name]);
    msa = msa(contains({msa.Header}, 'x'));
    [conSeq, ~] = seqconsensus(msa,'Gaps','all');
    if any(contains({msa.Header}, 'Zrou'))
        zrouSeq = msa(contains({msa.Header}, 'Zrou')).Sequence;
        [~, idxZrou] = ismember(zrouSeq', AAorder);
        [~, idxCon] = ismember(conSeq', AAorder);
        score = B62(sub2ind(size(B62), idxZrou,idxCon));
        gaps = strfind(zrouSeq, '-');
        score(gaps) = [];
        ZrouMSAscore.(extractBefore(ancestorMSA(a).name,'_')) = score';
    end
end

for a = 1:numel(ancestorMSA)
    msa = multialignread([ancestorMSA(a).folder,'/',ancestorMSA(a).name]);
    msa = msa(contains({msa.Header}, 'x'));
    [conSeq, ~] = seqconsensus(msa,'Gaps','all');
    if any(contains({msa.Header}, 'Egos'))
        EgosSeq = msa(contains({msa.Header}, 'Egos')).Sequence;
        [~, idxEgos] = ismember(EgosSeq', AAorder);
        [~, idxCon] = ismember(conSeq', AAorder);
        score = B62(sub2ind(size(B62), idxEgos,idxCon));
        gaps = strfind(EgosSeq, '-');
        score(gaps) = [];
        EgosMSAscore.(extractBefore(ancestorMSA(a).name,'_')) = score';
    end
end


%% Figure 6B - Ixr1 Abf2 tree and sequence alignment
DBDallParas = readtable('./allDBDpara.xls');
DBDallParas = DBDallParas(contains(DBDallParas.queryName, {'Scer','Klac','Zrou','Egos'}),:);
examplePairs = {'Ixr1','Abf2'};
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y','X', '-'};
B62 = blosum(62,'Order', cat(2,AAorder{1:end-1}));
B62(22,:) = -8;
B62(:,22) = -8;
colMapAlignment =brighten([1 1 1; flipud(bone)],0.8);
scatterCol = [51 51 51]/255%[148 175 211]/255;
markerCol = brighten(scatterCol, -0.5);
DBDboxCol = [135 48 146]/255;
intBases = intBasesVec(1:6701);
cMapCorr = brighten(flipud(bone),0.2);
colSpot = [[140 113 180]/255;[148 175 211]/255];

Walign = 0.07;
Wtree = 0.018;
Htree = 0.08;
xspacer = 0.005;
yspacer = 0.05;
xPos = [0.05];
yPos = [0.1, 0.4 0.7];
Wscatter = ((Wtree+Walign+5*xspacer)-xspacer)/3;
Hscatter = Wscatter*1.8;
Halign = (Htree-yspacer/5)/3;
colSpot = [[140 113 180]/255;[148 175 211]/255];

figure('Units','normalized','Position', [0 0 0.5 1])
    for i = 1:size(examplePairs,1)
    % Tree
    [order,b] = plotTrees(examplePairs(i,:),[xPos+Walign+6*xspacer yPos(i)+Hscatter+yspacer Wtree+0.08 Htree], gcf,'colSpot',colSpot);
    set(gcf, 'color','w')
    TF1 = examplePairs{i,order(1)};
    TF2 = examplePairs{i,order(2)};

    % sequence alignment
    if ismember(upper(TF1), extractfield(paraSeqs,'name'))
        [~, idx] = ismember(upper(TF1), extractfield(paraSeqs,'name'));
        pSeq = paraSeqs(idx).seq([22,12]) ;
    elseif ismember(upper(TF2), extractfield(paraSeqs,'name'))
        [~, idx] = ismember(upper(TF2), extractfield(paraSeqs,'name'));
        pSeq = paraSeqs(idx).seq([12,22]) ;
    end
    
    if ~isempty(paraSeqs(idx).seq{18})
        ancSeq = paraSeqs(idx).seq{18};
        MSAanc = 'lactisMSAscore';
        ancLabel = 'Klac';
        ancIdx = 18;
    elseif ~isempty(paraSeqs(idx).seq{20})
        ancSeq = paraSeqs(idx).seq{20};
        MSAanc = 'ZrouMSAscore';
        ancLabel = 'Zrou';
        ancIdx = 20;
    else
        ancSeq = paraSeqs(idx).seq{17};
        MSAanc = 'EgosMSAscore';
        ancLabel = 'Egos';
        ancIdx = 17;
    end
    maxLength = max(cellfun('prodofsize', [ancSeq;pSeq]));
    if ~isempty(ancSeq)
        clear posScore gaps
        for p =1:2
            [~,alignmentLac] = nwalign(pSeq{p}(1:end-1),ancSeq(1:end-1));
            [~, idxPara] = ismember(alignmentLac(1,:)', AAorder);
            [~, idxLac] = ismember(alignmentLac(3,:)', AAorder);
            posScore{p} = B62(sub2ind(size(B62), idxPara,idxLac));
            gaps{p} = strfind(alignmentLac(1,:), '-');
            posScoreNoGaps{p} = posScore{p}(setdiff(1:length(posScore{p}), gaps{p}));
        end
        imageMat = nan(2, maxLength);
        for p =1:2
            imageMat(p,1:length(posScoreNoGaps{p})) = posScoreNoGaps{p};
        end
        axes('Position', [xPos yPos(i)+Hscatter+1.2*yspacer+Halign Walign Halign*2])
        paraLength = sum(~isnan(imageMat),2);
        imageMat = movmean(imageMat,20,2, 'omitnan').*((imageMat+50)./(imageMat+50));
        imagesc(imageMat);
        colorLabel = sprintfc('\\color[rgb]{%.2f,%.2f,%.2f}',colSpot);
        set(gca, 'YTick',[1,2], 'YTickLabel', strcat(colorLabel','\fontsize{12}\bf', {TF1,TF2}), 'XTick', [], 'XColor', [1 1 1 0])
        
        % DBD annotations Scer
        hold on
        for p = 1:2
            currDBD = DBDallParas(strcmp(DBDallParas.queryName, sprintf('Scer_%s',upper(eval(sprintf('TF%d',p))))),:);
            for d = 1:size(currDBD,1)
                xVec = [currDBD.from_2(d),currDBD.to_2(d),currDBD.to_2(d), currDBD.from_2(d), currDBD.from_2(d)];
                yVec = [1.5, 1.5, 0.5, 0.5, 1.5]+p-1;
                plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
            end
            if p == 1
                text(mean([currDBD.from_2(1),currDBD.to_2(2)]), 0.3, 'DBD', 'color',DBDboxCol, 'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontWeight','bold')
            end
        end
        
        plot([0,max(paraLength)]+0.5, [1.5 1.5], 'k');
        plot(paraLength(1)*[1 1], [0.5, 1.5], 'k');
        plot(paraLength(2)*[1 1], [1.5, 2.5], 'k');
        plot([0,paraLength(1)]+0.5, [0.5 0.5], 'k');
        plot([0,paraLength(2)]+0.5, [2.5 2.5], 'k');
        colormap(gca, colMapAlignment);
        caxis([-8.5 7]);
        box off;
              
        % lactis - alignment
        axes('Position', [xPos yPos(i)+Hscatter+yspacer Walign Halign])
        if isfield(eval(MSAanc), upper(TF1))
            lacScore = eval(sprintf('%s.(upper(TF1))', MSAanc));
        elseif isfield(eval(MSAanc), upper(TF2))
            lacScore = eval(sprintf('%s.(upper(TF2))', MSAanc));
        end
        minYlim = min(movmean(lacScore,10,'omitnan'));
        maxYlim = max(movmean(lacScore,10,'omitnan'));
        imagesc(movmean(lacScore,20,'omitnan'));
        xticks([]);
        axis tight
        caxis([-8.5 7]);
        yticks(1)
        yticklabels(sprintf('\\it%s(%d)',ancLabel, maxLength))
        set(gca,'fontSize',12);
        colormap(gca, colMapAlignment);
        box on;
        
        % DBD annotations lac
        hold on
        currDBD = DBDallParas(contains(DBDallParas.queryName, {sprintf('%s_%s/%s',ancLabel, upper(TF1),upper(TF2)),sprintf('%s_%s/%s',ancLabel, upper(TF2),upper(TF1))}),:);
        for d = 1:size(currDBD,1)
            xVec = [currDBD.from_2(d),currDBD.to_2(d),currDBD.to_2(d), currDBD.from_2(d), currDBD.from_2(d)];
            yVec = [1.5, 1.5, 0.5, 0.5, 1.5];
            plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
        end
        
        cbr = colorbar()
        cbr.Position = [xPos+Walign+0.002 yPos(i)+Hscatter+yspacer 0.0048 Htree]
        set(cbr, 'Ticks', [min(caxis) max(caxis)], 'TickLabels', [])
        ylabel(cbr, 'Alignment\newline    score', 'fontSize',10)
    end
end
set(gcf, 'color','w')


%% Figures 6C and 7B - Examples: demonstrating neo/sub, sequence alignment, tree, scatters paralog deletion mutants vs K. lactis and corr matrix (sumProm,7mers)
DBDallParas = readtable('./allDBDpara.xls');
DBDallParas = DBDallParas(contains(DBDallParas.queryName, {'Scer','Klac','Zrou','Egos'}),:);
examplePairs = {'Rph1','Gis1'; 'Vhr1','Vhr2'; 'Tda9','Rsf2'};
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y','X', '-'};
B62 = blosum(62,'Order', cat(2,AAorder{1:end-1}));
B62(22,:) = -8;
B62(:,22) = -8;
colMapAlignment =brighten([1 1 1; flipud(bone)],0.8);
scatterCol = [51 51 51]/255%[148 175 211]/255;
markerCol = brighten(scatterCol, -0.5);
DBDboxCol = [135 48 146]/255;
intBases = createIntBasesForMotif();
cMapCorr = brighten(flipud(bone),0.2);
colSpot = [[140 113 180]/255;[148 175 211]/255];
%neoCol = [148 175 211]/255; %light blue
neoCol = [204 51 0]/255; % red

thCol = [.3 .3 .3];
colMapParalogsTargets = cbrewer2('BuPu');

Walign = 0.1;
Wtree =0.1;
Htree = 0.05;
xspacer = 0.005;
yspacer = 0.05;
xPos = [0.05];
yPos = [0.1, 0.4 0.7];
Wscatter = ((Wtree+Walign+5*xspacer)-xspacer)/3;
Hscatter = Wscatter*1.8;
Halign = (Htree-yspacer/5)/3;
zscoreTH = 3.5;
zscoreTHL = 2.5;
quantileTH = 1;

figure('Units','normalized','Position', [0 0 0.5 1], 'color','w')
for i = 1:size(examplePairs,1)
    % Tree
    [order,b] = plotTrees(examplePairs(i,:),[xPos+Walign+5*xspacer yPos(i)+Hscatter+yspacer Wtree+0.08 Htree], gcf,'colSpot',colSpot);
    set(gcf, 'color','w')
    TF1 = examplePairs{i,order(1)};
    TF2 = examplePairs{i,order(2)};

    % sequence alignment
    if ismember(upper(TF1), extractfield(paraSeqs,'name'))
        [~, idx] = ismember(upper(TF1), extractfield(paraSeqs,'name'));
        pSeq = paraSeqs(idx).seq([22,12]) ;
    elseif ismember(upper(TF2), extractfield(paraSeqs,'name'))
        [~, idx] = ismember(upper(TF2), extractfield(paraSeqs,'name'));
        pSeq = paraSeqs(idx).seq([12,22]) ;
    end
    if ~isempty(paraSeqs(idx).seq{18})
        ancSeq = paraSeqs(idx).seq{18};
        MSAanc = 'lactisMSAscore';
        ancLabel = 'Klac';
        ancIdx = 18;
    elseif ~isempty(paraSeqs(idx).seq{20})
        ancSeq = paraSeqs(idx).seq{20};
        MSAanc = 'ZrouMSAscore';
        ancLabel = 'Zrou';
        ancIdx = 20;
    else
        ancSeq = paraSeqs(idx).seq{17};
        MSAanc = 'EgosMSAscore';
        ancLabel = 'Egos';
        ancIdx = 17;
    end
    maxLength = max(cellfun('prodofsize', [ancSeq;pSeq]));
    if ~isempty(ancSeq)
        clear posScore gaps
        for p =1:2
            [~,alignmentLac] = nwalign(pSeq{p}(1:end-1),ancSeq(1:end-1));
            [~, idxPara] = ismember(alignmentLac(1,:)', AAorder);
            [~, idxLac] = ismember(alignmentLac(3,:)', AAorder);
            posScore{p} = B62(sub2ind(size(B62), idxPara,idxLac));
            gaps{p} = strfind(alignmentLac(1,:), '-');
            posScoreNoGaps{p} = posScore{p}(setdiff(1:length(posScore{p}), gaps{p}));
        end
        imageMat = nan(2, maxLength);
        for p =1:2
            imageMat(p,1:length(posScoreNoGaps{p})) = posScoreNoGaps{p};
        end
        axes('Position', [xPos yPos(i)+Hscatter+yspacer Walign Halign*2])
        paraLength = sum(~isnan(imageMat),2);
        imageMat = movmean(imageMat,20,2, 'omitnan').*((imageMat+50)./(imageMat+50));
        imagesc(imageMat);
        colorLabel = sprintfc('\\color[rgb]{%.2f,%.2f,%.2f}',colSpot);
        %set(gca, 'YTick',[1,2], 'YTickLabel', strcat(colorLabel','\fontsize{12}\bf', {TF1,TF2}), 'XTick', [], 'XColor', [1 1 1 0])
        set(gca, 'YTick',[1,2], 'YTickLabel', strcat(colorLabel','\fontsize{12}', {TF1,TF2}), 'XTick', [], 'XColor', [1 1 1 0])
        
        % DBD annotations Scer
        hold on
        for p = 1:2
            currDBD = DBDallParas(strcmp(DBDallParas.queryName, sprintf('Scer_%s',upper(eval(sprintf('TF%d',p))))),:);
            for d = 1:size(currDBD,1)
                xVec = [currDBD.from_2(d),currDBD.to_2(d),currDBD.to_2(d), currDBD.from_2(d), currDBD.from_2(d)];
                yVec = [1.5, 1.5, 0.5, 0.5, 1.5]+p-1;
                plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
            end
        end
        
        plot([0,max(paraLength)]+0.5, [1.5 1.5], 'k');
        plot(paraLength(1)*[1 1], [0.5, 1.5], 'k');
        plot(paraLength(2)*[1 1], [1.5, 2.5], 'k');
        plot([0,paraLength(1)]+0.5, [0.5 0.5], 'k');
        plot([0,paraLength(2)]+0.5, [2.5 2.5], 'k');
        colormap(gca, colMapAlignment);
        caxis([-8.5 7]);
        box off;
        
        % lactis - alignment
        axes('Position', [xPos yPos(i)+Hscatter+1.2*yspacer+2*Halign Walign Halign])
        if isfield(eval(MSAanc), upper(TF1))
            lacScore = eval(sprintf('%s.(upper(TF1))', MSAanc));
        elseif isfield(eval(MSAanc), upper(TF2))
            lacScore = eval(sprintf('%s.(upper(TF2))', MSAanc));
        end
        minYlim = min(movmean(lacScore,10,'omitnan'));
        maxYlim = max(movmean(lacScore,10,'omitnan'));
        imagesc(movmean(lacScore,20,'omitnan'));
        xticks([]);
        axis tight
        caxis([-8.5 7]);
        yticks(1)
        %yticklabels(sprintf('\\it%s(%d)',ancLabel, maxLength))
        yticklabels(sprintf('%s(%d)',ancLabel, maxLength))
        set(gca,'fontSize',12);
        colormap(gca, colMapAlignment);
        box on;
        
        % DBD annotations lac
        hold on
        currDBD = DBDallParas(contains(DBDallParas.queryName, {sprintf('%s_%s/%s',ancLabel, upper(TF1),upper(TF2)),sprintf('%s_%s/%s',ancLabel, upper(TF2),upper(TF1))}),:);
        for d = 1:size(currDBD,1)
            xVec = [currDBD.from_2(d),currDBD.to_2(d),currDBD.to_2(d), currDBD.from_2(d), currDBD.from_2(d)];
            yVec = [1.5, 1.5, 0.5, 0.5, 1.5];
            plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
        end
        
        cbr = colorbar()
        cbr.Position = [xPos+Wtree+0.002 yPos(i)+Hscatter+yspacer 0.0048 Htree]
        set(cbr, 'Ticks', [min(caxis) max(caxis)], 'TickLabels', [])
        ylabel(cbr, 'Alignment\newline     score', 'fontSize',10)
    end
       
    % scatter1
    axes('Position', [xPos yPos(i) Wscatter Hscatter])
    m1 = [TF1,'_d',upper(TF2)];
    m2 = [TF2,'_d',upper(TF1)];
    lac = [TF1, '_lactis'];
    zscoreMat = nan(6701,3);
    zscoreMat(:,1) =  nanZscore(checWTdelLactisSwap.sumProm.(m1))';
    zscoreMat(:,2) =  nanZscore(checWTdelLactisSwap.sumProm.(m2));
    zscoreMat(:,3) =  nanZscore(checWTdelLactisSwap.sumProm.(lac));
    neoLog = (zscoreMat(:,1) > zscoreTH) &  (zscoreMat(:,3) < zscoreTHL);
    nonLog = (zscoreMat(:,1) <= zscoreTH) &  (zscoreMat(:,3) <= zscoreTHL);
    oldLog=~neoLog&~nonLog;
    absThLac=mean(checWTdelLactisSwap.sumProm.(lac),'omitnan')+zscoreTHL*std(checWTdelLactisSwap.sumProm.(lac),[],'omitnan');
    absThM1=mean(checWTdelLactisSwap.sumProm.(m1),'omitnan')+zscoreTH*std(checWTdelLactisSwap.sumProm.(m1),[],'omitnan');
    absThM2=mean(checWTdelLactisSwap.sumProm.(m2),'omitnan')+zscoreTH*std(checWTdelLactisSwap.sumProm.(m2),[],'omitnan');
    maxX = max(checWTdelLactisSwap.sumProm.(m1));
    maxY = max(checWTdelLactisSwap.sumProm.(lac));   
%         scatter(checWTdelLactisSwap.sumProm.(m1)(~nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(~nonLog)/maxY,...
%         20, checWTdelLactisSwap.sumProm.(m2)(~nonLog)/max(checWTdelLactisSwap.sumProm.(m2)), 'filled', 'MarkerEdgeColor', 'k')
        scatter(checWTdelLactisSwap.sumProm.(m1)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,...
        20, checWTdelLactisSwap.sumProm.(m2)(oldLog)/max(checWTdelLactisSwap.sumProm.(m2)), 'filled', 'MarkerEdgeColor', 'k')
        hold on
        scatter(checWTdelLactisSwap.sumProm.(m1)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,...
        20, checWTdelLactisSwap.sumProm.(m2)(neoLog)/max(checWTdelLactisSwap.sumProm.(m2)), 'filled', 'MarkerEdgeColor', neoCol, 'LineWidth', 1)
        scatter(checWTdelLactisSwap.sumProm.(m1)(nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(nonLog)/maxY,...
        20, brighten(scatterCol,0.6), 'filled', 'MarkerFaceAlpha',0.5)
    %scatter(checWTdelLactisSwap.sumProm.(m1)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,20, scatterCol, 'filled','MarkerEdgeColor', markerCol)
    %scatter(checWTdelLactisSwap.sumProm.(m1)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,20, neoCol, 'filled','MarkerEdgeColor', markerCol)
    %xlabel(strrep(m1,'_d', ' \Delta'), 'fontSize',11)
    xlabel(strrep(m1,'_', ' '), 'fontSize',11)
    ylabel('K.lac','fontSize',11)
    crValue = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(lac),'rows','pairwise');
    Fneo = 100*sum(neoLog)/sum((zscoreMat(:,1) > min(quantile(zscoreMat(:,1), quantileTH), zscoreTH)));
        text(0.05, 1.1, ...
        sprintf('r=%.2f, \\color[rgb]{%.2f, %.2f, %.2f}new=%.0f%%', crValue,brighten(neoCol,-0.4), Fneo), 'FontSize',11,'VerticalAlignment','bottom', 'HorizontalAlignment', 'left')
    plot(xlim.*[1], absThLac.*[1 1]/maxY, '--','Color', thCol)
    plot(absThM1.*[1 1]/maxX, ylim.*[1], '--','Color', thCol)
    set(gca, 'XTick', [0:0.5:1], 'YTick', [0:0.5:1])
    %fill([absThM1.*[1 1]/maxX, 1,1], [0, absThLac.*[1 1]/maxY,0], brighten(neoCol,-0.4), 'FaceAlpha', 0.3, 'LineStyle', 'none')
    caxis([0 0.8])
    colormap(gca, colMapParalogsTargets)
     
%     scatter(checWTdelLactisSwap.sumProm.(m1)(nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(nonLog)/maxY,...
%         20, brighten(scatterCol,0.6), 'filled', 'MarkerFaceAlpha',0.5)
%     hold on
%     scatter(checWTdelLactisSwap.sumProm.(m1)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,20, scatterCol, 'filled','MarkerEdgeColor', markerCol)
%     scatter(checWTdelLactisSwap.sumProm.(m1)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,20, neoCol, 'filled','MarkerEdgeColor', markerCol)
%     xlabel(strrep(m1,'_d', ' \Delta'), 'fontSize',11)
%     ylabel('\itK.lac','fontSize',11)
%     ylim(1.2*ylim);
%     crValue = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(lac),'rows','pairwise');
%     Fneo = 100*sum(neoLog)/sum((zscoreMat(:,1) > min(quantile(zscoreMat(:,1), quantileTH), zscoreTH)));
%         text(0.05*max(xlim), 0.99*max(ylim), ...
%         sprintf('r=%.2f, \\color[rgb]{%.2f, %.2f, %.2f}new=%.0f%%', crValue,brighten(neoCol,-0.4), Fneo), 'FontSize',11,'VerticalAlignment','top')
%     
%     plot(xlim.*[0.9], absThLac.*[1 1]/maxY, '--','Color', thCol)
%     plot(absThM1.*[1 1]/maxX, ylim.*[0.9], '--','Color', thCol)
%     set(gca, 'XTick', [0:0.5:1], 'YTick', [0:0.5:1])

    % scatter2
    axes('Position', [xPos+Wscatter+xspacer yPos(i) Wscatter Hscatter])
    neoLog = (zscoreMat(:,2) > zscoreTH) &  (zscoreMat(:,3) < zscoreTHL);
    nonLog = (zscoreMat(:,2) <= zscoreTH) &  (zscoreMat(:,3) <= zscoreTHL);
    oldLog=~neoLog&~nonLog;
    maxX = max(checWTdelLactisSwap.sumProm.(m2));
    maxY = max(checWTdelLactisSwap.sumProm.(lac));
    
%     scatter(checWTdelLactisSwap.sumProm.(m2)(~nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(~nonLog)/maxY,...
%         20, checWTdelLactisSwap.sumProm.(m1)(~nonLog)/max(checWTdelLactisSwap.sumProm.(m1)), 'filled', 'MarkerEdgeColor', 'k')
    scatter(checWTdelLactisSwap.sumProm.(m2)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,...
        20, checWTdelLactisSwap.sumProm.(m1)(oldLog)/max(checWTdelLactisSwap.sumProm.(m1)), 'filled', 'MarkerEdgeColor', 'k')    
hold on
    scatter(checWTdelLactisSwap.sumProm.(m2)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,...
        20, checWTdelLactisSwap.sumProm.(m1)(neoLog)/max(checWTdelLactisSwap.sumProm.(m1)), 'filled', 'MarkerEdgeColor', neoCol, 'LineWidth',1)    
    scatter(checWTdelLactisSwap.sumProm.(m2)(nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(nonLog)/maxY,20, brighten(scatterCol,0.6), 'filled', 'MarkerFaceAlpha',0.5)
    %scatter(checWTdelLactisSwap.sumProm.(m2)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,20, scatterCol, 'filled','MarkerEdgeColor', markerCol)
    %scatter(checWTdelLactisSwap.sumProm.(m2)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,20, neoCol, 'filled','MarkerEdgeColor', markerCol)
    %xlabel(strrep(m2,'_d', ' \Delta'), 'fontSize',11)
    xlabel(strrep(m2,'_', ' '), 'fontSize',11)
    yticks([])
    Fneo = 100*sum(neoLog)/sum((zscoreMat(:,2) > min(quantile(zscoreMat(:,2), quantileTH), zscoreTH)));
    crValue = corr(checWTdelLactisSwap.sumProm.(m2), checWTdelLactisSwap.sumProm.(lac),'rows','pairwise');
    text(0.05, 1.1, ...
        sprintf('r=%.2f, \\color[rgb]{%.2f, %.2f, %.2f}new=%.0f%%', crValue,brighten(neoCol,-0.4), Fneo), 'FontSize',11,'VerticalAlignment','bottom', 'HorizontalAlignment', 'left')
    
    plot(xlim.*[1], absThLac.*[1 1]/maxY,'--','Color',thCol)
    plot(absThM2.*[1 1]/maxX,ylim.*[1],'--','Color',thCol)
    %fill([absThM2.*[1 1]/maxX, 1,1], [0, absThLac.*[1 1]/maxY,0], brighten(neoCol,-0.4), 'FaceAlpha', 0.3, 'LineStyle', 'none')
    caxis([0 0.8])
    colormap(gca, colMapParalogsTargets)
    
%     scatter(checWTdelLactisSwap.sumProm.(m2)(nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(nonLog)/maxY,20, brighten(scatterCol,0.6), 'filled', 'MarkerFaceAlpha',0.5)
%     hold on
%     scatter(checWTdelLactisSwap.sumProm.(m2)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,20, scatterCol, 'filled','MarkerEdgeColor', markerCol)
%     scatter(checWTdelLactisSwap.sumProm.(m2)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,20, neoCol, 'filled','MarkerEdgeColor', markerCol)
%     
%     xlabel(strrep(m2,'_d', ' \Delta'), 'fontSize',11)
%     yticks([])
%     ylim(1.2*ylim);
%     Fneo = 100*sum(neoLog)/sum((zscoreMat(:,2) > min(quantile(zscoreMat(:,2), quantileTH), zscoreTH)));
%        crValue = corr(checWTdelLactisSwap.sumProm.(m2), checWTdelLactisSwap.sumProm.(lac),'rows','pairwise')
%     text(0.05*max(xlim), 0.99*max(ylim), ...
%         sprintf('r=%.2f, \\color[rgb]{%.2f, %.2f, %.2f}new=%.0f%%', crValue,brighten(neoCol,-0.4), Fneo), 'FontSize',11,'VerticalAlignment','top')
%     
%     plot(xlim.*[0.9], absThLac.*[1 1]/maxY,'--','Color',thCol)
%     plot(absThM2.*[1 1]/maxX,ylim.*[0.9],'--','Color',thCol)  
    
    % corr matrix
    axes('Position', [xPos+2*(Wscatter+xspacer)+0.04 yPos(i) Wscatter Hscatter])
    intStrains = {lac, TF1,m1,TF2,m2};
    % WITH REPEATS
    %     [corrSumProm, idxSumProm] = getRepeatsCorr(intStrains, 'dataType','sumProm');
    %     [corr7mer, ~] = getRepeatsCorr(intStrains, 'dataType','7mer');
    %     combineMat = tril(corrSumProm) + triu(corr7mer,1);
    %     imagesc(combineMat)
    %     borders = [0; cumsum(accumarray(idxSumProm,1))]+0.5;
    %     tickPos = movmean(borders,2,'Endpoints','discard');
    %     hold on
    %     plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
    %     plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
    
    % WITHOUT REPEATS
    corrSumProm = plotSumPromCorr(intStrains, checWTdelLactisSwap,0);
    normProfile = chromosomes2fullProfile(checWTdelLactisSwap, intStrains);
    mers = mer_occupancy(normProfile, 7, 'method','else','intBases',intBases);
    corr7mer = corr(mers.score);
    combinedMat = tril(corrSumProm) + triu(corr7mer,1);
    imagesc(combinedMat)
    plotgrid(combinedMat)
    plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)
    colormap(gca, cMapCorr)
    caxis([0.3 1])
    %set(gca, 'YTick', [1:numel(intStrains)], 'YTickLabel',regexprep(intStrains, {'_d','.*_lactis'},{' \\Delta', '\\itK.lac'}), 'XTick',[], 'fontSize',9)
    set(gca, 'YTick', [1:numel(intStrains)], 'YTickLabel', strrep(intStrains,'_',' '), 'XTick',[], 'fontSize',9)
    hold on
    plot(xlim, ylim, 'color',[0.7 0.7 0.7])
    cbr = colorbar()
    cbr.Position = [xPos+3*(Wscatter+xspacer)+0.04 yPos(i) 0.0048 Hscatter]
    set(cbr, 'Ticks', [0.4:0.2:1], 'fontSize',11)
    ylabel(cbr,'correlation', 'fontSize',15)    
        if i ==1
            axes('Position', [xPos yPos(i)-0.06 2*Wscatter+xspacer 0.05])
            caxis([0 0.8])
            colormap(gca, colMapParalogsTargets);
            cbr = colorbar('Location', 'southoutside')
            axis off
            set(cbr, 'Ticks', [0:0.4:0.8], 'fontSize',11)
            ylabel(cbr, 'Paralog signal (in mutant background)', 'fontSize',12)            
        end
end
% saveas(gcf, 'Fig4Examples.svg')


%% Figure 7C - summary plot
% calculating parameters based on trees 
load('promCorrSort.mat')
groups = fieldnames(promCorrSort);
clear m1Closer t1 t2 zscoreMat m1Neo m2Neo
FN = fieldnames(checWTdelLactisSwap.sumProm);
lacSamples = FN(find(contains(FN,'_lactis') & ~contains(FN, promCorrSort.specialCases.TFsSort(:))&~contains(FN,'Sut1')));
zscoreTH = 4.5;
zscoreTHL = 3.5;

m1LacCorr = nan(length(lacSamples),1);
m2LacCorr = nan(length(lacSamples),1);
m1m2Corr = nan(length(lacSamples),1);
WTsCorr = nan(length(lacSamples),1);
meanNonDBDancesVec = nan(length(lacSamples),2);
WTflag = [];
t = {};
c=1;
for g = 1:length(groups)-1
    nParalogs = size(promCorrSort.(groups{g}).TFsSort,1);
    for i = 1:nParalogs
        currF1 = promCorrSort.(groups{g}).TFsSort{i,1};
        currF2 = promCorrSort.(groups{g}).TFsSort{i,2};
        m1 = [currF1, '_d', upper(currF2)];
        m2 = [currF2, '_d', upper(currF1)];
        lacIdx =  [find(strcmp(lacSamples, [currF1,'_lactis'])), find(strcmp(lacSamples, [currF2,'_lactis']))];
        curr_m1LacCorr = [];
        curr_m2LacCorr = [];
        curr_y1Err = [];
        curr_y2Err = [];
        %xErr = [];
        WTsCorr(c) =  corr(checWTdelLactisSwap.sumProm.(currF1), checWTdelLactisSwap.sumProm.(currF2), 'rows','pairwise');
        if ~isempty(lacIdx)
            currLac = lacSamples{lacIdx};
            zscoreMat(3,:) = nanZscore(checWTdelLactisSwap.sumProm.(currLac));
            if isfield(checWTdelLactisSwap.sumProm,m1)
                % curr_m1LacCorr = corr(checWTdelLactisSwap.sumProm.(currLac), checWTdelLactisSwap.sumProm.(m1), 'rows','pairwise');
                [curr_m1LacCorr, curr_y1Err] = corrBetweenRepeatsOf2Samples({currLac,m1},'dataType', 'sumProm');
                zscoreMat(1,:) = nanZscore(checWTdelLactisSwap.sumProm.(m1));
                WTflag(c,1)=0;
            else
                % curr_m1LacCorr = corr(checWTdelLactisSwap.sumProm.(currLac), checWTdelLactisSwap.sumProm.(currF1), 'rows','pairwise');
                [curr_m1LacCorr, curr_y1Err] = corrBetweenRepeatsOf2Samples({currLac,currF1},'dataType', 'sumProm');
                WTflag(c,1)=1;
                zscoreMat(1,:) = nanZscore(checWTdelLactisSwap.sumProm.(currF1));
                % m1m2Corr(c) = corr(checWTdelLactisSwap.sumProm.(currF1), checWTdelLactisSwap.sumProm.(m2), 'rows','pairwise');
                 [m1m2Corr(c), xErr(c)] = corrBetweenRepeatsOf2Samples({currF1,m2},'dataType', 'sumProm');
            end
            if isfield(checWTdelLactisSwap.sumProm,m2)
                % curr_m2LacCorr = corr(checWTdelLactisSwap.sumProm.(currLac), checWTdelLactisSwap.sumProm.(m2), 'rows','pairwise');
                [curr_m2LacCorr, curr_y2Err] = corrBetweenRepeatsOf2Samples({currLac,m2},'dataType', 'sumProm');
                zscoreMat(2,:) = nanZscore(checWTdelLactisSwap.sumProm.(m2));
            else
                % curr_m2LacCorr = corr(checWTdelLactisSwap.sumProm.(currLac), checWTdelLactisSwap.sumProm.(currF2), 'rows','pairwise');
               [curr_m2LacCorr, curr_y2Err] = corrBetweenRepeatsOf2Samples({currLac,currF2},'dataType', 'sumProm');
                zscoreMat(2,:) = nanZscore(checWTdelLactisSwap.sumProm.(currF2));
                WTflag(c,2)=1;
                % m1m2Corr(c) = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(currF2), 'rows','pairwise');
                 [m1m2Corr(c), xErr(c)] = corrBetweenRepeatsOf2Samples({m1,currF2},'dataType', 'sumProm');
            end
            if sum(WTflag(c,:)) == 0
                % m1m2Corr(c) = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(m2), 'rows','pairwise');
                 [m1m2Corr(c), xErr(c)] = corrBetweenRepeatsOf2Samples({m1,m2},'dataType', 'sumProm');
            end
            m1Neo(c) = sum(zscoreMat(1,:)>zscoreTH & zscoreMat(3,:)<zscoreTHL)/sum(zscoreMat(1,:)>zscoreTH);
            m2Neo(c) = sum(zscoreMat(2,:)>zscoreTH & zscoreMat(3,:)<zscoreTHL)/sum(zscoreMat(2,:)>zscoreTH);
            
            m1LacCorr(c) = curr_m1LacCorr;
            m2LacCorr(c) = curr_m2LacCorr;
            y1Err(c) = curr_y1Err
            y2Err(c) = curr_y2Err
            t1{c} = currF1;
            t2{c} = currF2;
            m1Closer(c) = ismember(upper(extractBefore(m1,'_')), extractAfter(treeFiles.Name1(treeFiles.order == 3),'_'));
            c = c+1;
        end
    end
end


% summary plot with Std
figure
subplot(1,2,1)
for i = 1: length(m1m2Corr)
    plot([m1m2Corr(i), m1m2Corr(i)], [m1LacCorr(i),m2LacCorr(i)],...
        'Color', [0.9 0.9 0.9 ], 'LineWidth', 1.5);
        hold on
end
cMap = brighten(flipud(bone),0.3);
capSize=3;
[maxCorr,maxIdx] = max([m1LacCorr,m2LacCorr],[],2) ;
[minCorr,minIdx] = min([m1LacCorr,m2LacCorr],[],2) ;
maxCloser = (maxIdx == 1 & m1Closer') |  (maxIdx == 2 & ~m1Closer');
sc = scatter(m1m2Corr,maxCorr, 40+120*maxCloser, WTsCorr, 'filled', 'MarkerEdgeColor','k');
% errorbar(m1m2Corr, maxCorr, y1Err,y1Err, xErr, xErr, 'LineStyle','none', 'color','k', 'lineWidth', 1.5)
errorbar(m1m2Corr, maxCorr, y1Err,y1Err,xErr,xErr,  'LineStyle','none', 'color','k', 'lineWidth', 1.1,'CapSize',capSize)
hold on
sc = scatter(m1m2Corr,minCorr, 40+120*(~maxCloser), WTsCorr, 'filled','MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0.1,'MarkerEdgeColor','k');
% errorbar(m1m2Corr, minCorr, y2Err,y2Err, xErr, xErr, 'LineStyle','none', 'color','k', 'lineWidth', 1.5)
errorbar(m1m2Corr, minCorr, y2Err,y2Err,xErr,xErr, 'LineStyle','none', 'color',[0.5 0.5 0.5], 'lineWidth', 1.1,'CapSize',capSize)
axis square
set(gcf, 'Color','w')
xlim([0 1])
ylim([0 1])
hold on
plot([1 0], [1 0], '--k')
set(gca, 'FontSize',15, 'XTick',[0:0.2:1], 'YTick',[0:0.2:1])
cbr = colorbar();
ylabel(cbr, sprintf('Paralog-Paralog\n (WT background)'), 'FontSize',20);
ylabel('Paralog-Ortholog Corr', 'FontSize',20);
xlabel('Paralog-Paralog Corr', 'FontSize',20);
set(cbr, 'YTick',[0:0.2:1])
colormap(gca, cMap)
box on
t12=[t1;t2];
text(m1m2Corr(m1m2Corr<0.8)+0.015,maxCorr(m1m2Corr<0.8), t12(sub2ind(size(t12),maxIdx(m1m2Corr<0.8),find(m1m2Corr<0.8))), 'FontSize',15)
text(m1m2Corr(m1m2Corr<0.8)+0.015,minCorr(m1m2Corr<0.8),  t12(sub2ind(size(t12),minIdx(m1m2Corr<0.8),find(m1m2Corr<0.8))), 'FontSize',15, 'Color', [0.7 0.7 0.7])
Rsf2Tda9Idx = find(contains(t1,{'Rsf2','Tda9'}));
text(m1m2Corr(Rsf2Tda9Idx)+0.015,maxCorr(Rsf2Tda9Idx), t1(Rsf2Tda9Idx), 'FontSize',15)
text(m1m2Corr(Rsf2Tda9Idx)+0.015,minCorr(Rsf2Tda9Idx), t2(Rsf2Tda9Idx), 'FontSize',15, 'Color', [0.7 0.7 0.7])
% saveas(gcf, 'Fig4SummaryPlotWithSTD.svg')


%% Figure 7D - sub vs neo plot
xvalue = abs(m1LacCorr - m2LacCorr);
% barColM =  cbrewer2('BuPu');
barColM =  flipud(gray);
colNeoSub = [0.2 0.2 0.2; 0.8 0.8 0.8];
colByNeoFraction = (m2LacCorr>m1LacCorr)'.*m1Neo + (m1LacCorr>m2LacCorr)'.*m2Neo;

figure('Units','normalized','Position', [0 0 0.5 1], 'Color','w')
ax = axes('Position', [0.05 0.6 0.2 0.05])
%scatter(xvalue, repmat(0, numel(xvalue),1), (rescale(1-m1m2Corr, 3,15)).^2, 1-WTsCorr, 'filled', 'MarkerEdgeColor', barColM(150,:))
scatter(xvalue, repmat(0, numel(xvalue),1), (rescale(1-m1m2Corr, 3,15)).^2, colByNeoFraction*100, 'filled', 'MarkerEdgeColor', barColM(150,:))
ylim([-0.5 0.5])
xlabel('\Delta correlations with \itK.lac', 'fontSize',15)
colormap(gca,barColM)
% cbr = colorbar
% ylabel(cbr, 'WTs divergence')
% set(cbr, 'Ticks', [0:0.4:1])
caxis([0 100])
set(gca, 'XTick', [0:0.1:max(xlim)], 'fontSize',12, 'YTick', [],'TickLength',[0 0])
box on
cbr = colorbar();
cbr.Position = [0.26 0.6 0.005 0.05]
ylabel(cbr, '% new', 'fontSize',15)

axes('Position', ax.Position+[0 ax.Position(4)+0.01 0 -ax.Position(4)/2])
fill([0 1 1 0], [0 0 1 0], colNeoSub(1,:), 'EdgeColor', [1 1 1])
hold on
fill([0 1 0 0], [1 1 0 1], colNeoSub(2,:),'EdgeColor', [1 1 1])
axis off
text(0, 1, ' sub', 'Color', [1 1 1], 'HorizontalAlignment','left', 'VerticalAlignment','top','fontSize',11,'FontWeight','bold')
text(1, 0, 'neo ', 'Color', [1 1 1], 'HorizontalAlignment','right', 'VerticalAlignment','bottom','fontSize',11,'FontWeight','bold')
