%% Figure 6—figure supplement 1 and Figure 7—figure supplement 1
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
load('SC_genome.mat')
load('Abf2Ixr1.mat')
load('./treeFiles.mat');
load('./paraSeqs.mat');
load('./paraSeqs.mat');

%% Figure 7—figure supplement 1A: histogram promoter binding signal correlations all samples vs repeats
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
lacNames = allSamples(contains(allSamples, '_lactis')& ~contains(allSamples,{'Sut1','Rsf2','Vhr1','Rph1','Ixr1','Nhp6B','Nfi1'}));
[corrMat, idxVec, fullMat, idxRep] = getRepeatsCorr(lacNames, 'dataType','sumProm');
[~, firstOcc] = unique(idxRep,'stable');
redFullMat = fullMat(:, firstOcc);
allSamplesCorr = 1-squareform(1-corr(redFullMat, 'rows','pairwise'));

clear repeatsCorr
for i = 1:numel(lacNames)
    repeatsCorr{i} = 1-squareform(1-corr(redFullMat(:,idxVec(firstOcc)==i),'rows','pairwise'));
end
repeatsCorr = cat(2,repeatsCorr{:});

% histograms - promoter signal
col = [[156 186 217]/255; [76 0 128]/255];
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
axes('Position', [0.1 0.5 0.25 0.3])
h1 = histogram(allSamplesCorr, 'FaceColor', col(1,:));
hold on
h2 = histogram(repeatsCorr, 'FaceColor', col(2,:));
h1.Normalization = 'probability';
h1.BinWidth = 0.02
h2.Normalization = 'probability';
h2.BinWidth = 0.02
title('Promoter signal','FontWeight','bold','FontSize',18)
legend({'All profiles', 'Repeats'}, 'Location','northwest');
xlabel('Correlation','FontSize',18);
ylabel('Probability', 'FontSize',18);
set(gca, 'FontSize',18, 'TickLength',[0 0])
xlim([-0.1 1])


%% Figure 6—figure supplement 1B: all protein sequence alignment with K.lactis and Z.rouxii - plot trees & alignments
% scoring lactis sequnce by multiple alignment sequence with ancestors consensus 
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

% plot trees + alignments
DBDallParas = readtable('./allDBDpara.xls');
DBDallParas = DBDallParas(contains(DBDallParas.queryName, {'Scer','Klac','Zrou','Egos'}),:);
examplePairs = [summaryTable.p1,summaryTable.p2];
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
Wtree = 0.02;
Htree = 0.06;
xspacer = 0.005;
yspacer = 0.05;
xPos = repmat([0.05, 0.3,0.55],10,1);
yPos = repmat([0.9:-0.08:0.18]',1,3);
Wscatter = ((Wtree+Walign+5*xspacer)-xspacer)/3;
Hscatter = Wscatter*1.8;
Halign = (Htree-yspacer/5)/3;
colSpot = [[140 113 180]/255;[148 175 211]/255];

figure('Units','normalized','Position', [0 0 0.5 1])
    for i = 1:size(examplePairs,1)
    % Tree
    [order,b] = plotTrees(examplePairs(i,:),[xPos(i)+Walign+3*xspacer yPos(i) Wtree+0.08 Htree], gcf,'colSpot',colSpot);
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
        axes('Position', [xPos(i) yPos(i) Walign Halign*2])
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
                %plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
            end
            if p == 1
                %text(mean([currDBD.from_2(1),currDBD.to_2(2)]), 0.3, 'DBD', 'color',DBDboxCol, 'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontWeight','bold')
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
        axes('Position', [xPos(i) yPos(i)+Halign*2+0.005 Walign Halign])
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
            %plot(xVec, yVec, 'color', DBDboxCol, 'LineWidth',2)
        end
        
%         cbr = colorbar()
%         cbr.Position = [xPos(i)+Walign+0.002 yPos(i) 0.0048 Htree]
%         set(cbr, 'Ticks', [min(caxis) max(caxis)], 'TickLabels', [])
%         ylabel(cbr, 'Alignment\newline    score', 'fontSize',10)
    end
end
set(gcf, 'color','w')
axes('Position', [0.04 0.14 0.2 0.05])
cbr = colorbar('Location', 'south')
colormap(gca,colMapAlignment )
set(cbr,'AxisLocation','out')
set(cbr, 'Ticks', [min(caxis) max(caxis)], 'TickLabels', {'min','max'}, 'fontSize',12)
ylabel(cbr, 'Alignment score', 'fontSize',15)
axis off


%% Figure 7—figure supplement 1B: scatters deletion mutants against K.lactis and corr matrix (sumProm,motifs)
DBDallParas = readtable('./allDBDpara.xls');
DBDallParas = DBDallParas(contains(DBDallParas.queryName, {'Scer','Klac','Zrou','Egos'}),:);
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
lacNames = allSamples(contains(allSamples, '_lactis')& ~contains(allSamples,{'Sut1','Ixr1','Nhp6B','Nfi1'}));
examplePairsIdx = find(ismember(summaryTable.p1, extractBefore(lacNames,'_lactis')) | ismember(summaryTable.p2, extractBefore(lacNames,'_lactis')));
examplePairs = [summaryTable.p1(examplePairsIdx),summaryTable.p2(examplePairsIdx)];

DBDallParas = readtable('./allDBDpara.xls');
DBDallParas = DBDallParas(contains(DBDallParas.queryName, {'Scer','Klac','Zrou','Egos'}),:);
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y','X', '-'};
B62 = blosum(62,'Order', cat(2,AAorder{1:end-1}));
B62(22,:) = -8;
B62(:,22) = -8;
colMapAlignment =brighten([1 1 1; flipud(bone)],0.8);
scatterCol = [51 51 51]/255%[148 175 211]/255;
markerCol = brighten(scatterCol, -0.5);
DBDboxCol = [135 48 146]/255;
intBases = intBasesVec(1:6701);
cMapCorr = brighten(flipud(bone(512)),0.2);
colSpot = [[140 113 180]/255;[148 175 211]/255];
neoCol = [148 175 211]/255;
goodPos = createIntBasesForMotif();

thCol = [.3 .3 .3];
zscoreTH = 4.5;
zscoreTHL = 3.5;
quantileTH = 1;

% Walign = 0.1;
% Wtree =0.1;
% Htree = 0.05;
% xspacer = 0.005;
% yspacer = 0.05;
xPos = repmat([0.05,0.32,0.59],6,1);
%Wscatter = ((Wtree+Walign+5*xspacer)-xspacer)/3;
Wscatter = 0.06;
Hscatter = Wscatter*1.8;
yPos = repmat([0.85:-(Hscatter+0.05):0.9-(Hscatter+0.05)*6]',1,3);
Halign = (Htree-yspacer/5)/3;

figure('Units','normalized','Position', [0 0 0.5 1],'color','w')
for i = 1:size(examplePairs,1)
    % Tree
    %[order,b] = plotTrees(examplePairs(i,:),[xPos+Walign+5*xspacer yPos(i)+Hscatter+yspacer Wtree+0.08 Htree], gcf,'colSpot',colSpot);
    TF1 = examplePairs{i,1};
    TF2 = examplePairs{i,2};
    lac = allSamples{contains(allSamples, {TF1,TF2})&contains(allSamples, '_lactis')};
    if corr(checWTdelLactisSwap.sumProm.(TF2),checWTdelLactisSwap.sumProm.(lac),'rows','pairwise')>corr(checWTdelLactisSwap.sumProm.(TF1),checWTdelLactisSwap.sumProm.(lac),'rows','pairwise')
        TF1 = examplePairs{i,2};
        TF2 = examplePairs{i,1};
    end
    
      % scatter1
      axes('Position', [xPos(i) yPos(i) Wscatter Hscatter])
      m1 = [TF1,'_d',upper(TF2)];
      if ~any(strcmp(allSamples, m1))
          m1 = TF1;
      end
      m2 = [TF2,'_d',upper(TF1)];
      if ~any(strcmp(allSamples, m2))
          m2 = TF2;
      end 
    lac = allSamples{ismember(allSamples, [TF1,'_lactis']) | ismember(allSamples, [TF2,'_lactis'])};
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
    
    scatter(checWTdelLactisSwap.sumProm.(m1)(nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(nonLog)/maxY,...
        20, brighten(scatterCol,0.6), 'filled', 'MarkerFaceAlpha',0.5)
    hold on
    scatter(checWTdelLactisSwap.sumProm.(m1)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,20, scatterCol, 'filled','MarkerEdgeColor', markerCol)
    scatter(checWTdelLactisSwap.sumProm.(m1)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,20, neoCol, 'filled','MarkerEdgeColor', markerCol)
    xlabel(strrep(m1,'_d', ' \Delta'), 'fontSize',10)
    ylabel('\itK.lac','fontSize',10)
    ylim(1.2*ylim);
    crValue = corr(checWTdelLactisSwap.sumProm.(m1), checWTdelLactisSwap.sumProm.(lac),'rows','pairwise');
    Fneo = 100*sum(neoLog)/sum((zscoreMat(:,1) > min(quantile(zscoreMat(:,1), quantileTH), zscoreTH)));
        text(0.05, 1.2, ...
        sprintf('r=%.2f, \\color[rgb]{%.2f, %.2f, %.2f}new=%.0f%%', crValue,brighten(neoCol,-0.4), Fneo), 'FontSize',10,'VerticalAlignment','top', 'HorizontalAlignment', 'left')
    
    plot([0 1], absThLac.*[1 1]/maxY, '--','Color', thCol)
    plot(absThM1.*[1 1]/maxX, [0 1], '--','Color', thCol)
    set(gca, 'XTick', [0:0.5:1], 'YTick', [0:0.5:1])
    
    % scatter2
    axes('Position', [xPos(i)+Wscatter+0.005 yPos(i) Wscatter Hscatter])
    neoLog = (zscoreMat(:,2) > zscoreTH) &  (zscoreMat(:,3) < zscoreTHL);
    nonLog = (zscoreMat(:,2) <= zscoreTH) &  (zscoreMat(:,3) <= zscoreTHL);
    oldLog=~neoLog&~nonLog;
    maxX = max(checWTdelLactisSwap.sumProm.(m2));
    maxY = max(checWTdelLactisSwap.sumProm.(lac));
    scatter(checWTdelLactisSwap.sumProm.(m2)(nonLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(nonLog)/maxY,20, brighten(scatterCol,0.6), 'filled', 'MarkerFaceAlpha',0.5)
    hold on
    scatter(checWTdelLactisSwap.sumProm.(m2)(oldLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(oldLog)/maxY,20, scatterCol, 'filled','MarkerEdgeColor', markerCol)
    scatter(checWTdelLactisSwap.sumProm.(m2)(neoLog)/maxX, checWTdelLactisSwap.sumProm.(lac)(neoLog)/maxY,20, neoCol, 'filled','MarkerEdgeColor', markerCol)
    xlabel(strrep(m2,'_d', ' \Delta'), 'fontSize',10)
    yticks([])
    ylim(1.2*ylim);
    Fneo = 100*sum(neoLog)/sum((zscoreMat(:,2) > min(quantile(zscoreMat(:,2), quantileTH), zscoreTH)));
       crValue = corr(checWTdelLactisSwap.sumProm.(m2), checWTdelLactisSwap.sumProm.(lac),'rows','pairwise');
    text(0.05, 1.2, ...
        sprintf('r=%.2f, \\color[rgb]{%.2f, %.2f, %.2f}new=%.0f%%', crValue,brighten(neoCol,-0.4), Fneo), 'FontSize',10,'VerticalAlignment','top', 'HorizontalAlignment', 'left')  
    plot([0 1], absThLac.*[1 1]/maxY,'--','Color',thCol)
    plot(absThM2.*[1 1]/maxX, [0 1],'--','Color',thCol)
      
    % corr matrix
    axes('Position', [xPos(i)+2*(Wscatter)+0.05 yPos(i) Wscatter Hscatter])
    m1 = [TF1,'_d',upper(TF2)];
    m2 = [TF2,'_d',upper(TF1)];
    intStrains = {lac, TF1,m1,TF2,m2};
    [~, idxVec, sumPromRep, repIdx] = getRepeatsCorr(intStrains,'dataType','sumProm');
    [~, ~, mer7Rep, ~] = getRepeatsCorr(intStrains,'dataType','7mer');
    samplesWOrepeats = unique(idxVec(all(isnan(sumPromRep),1)));
    
    for z = samplesWOrepeats'
        if isfield(checWTdelLactisSwap.sumProm, intStrains{z})
        sumPromRep(:,idxVec==z) = repmat(checWTdelLactisSwap.sumProm.(intStrains{z}),1,sum(idxVec==z));
        currNorm = chromosomes2fullProfile(checWTdelLactisSwap, intStrains(z));
        currMer = mer_occupancy(currNorm,7,'intBases', goodPos,'method','normal');
        mer7Rep(:,idxVec==z) = repmat(currMer.score,1,sum(idxVec==z));
        end
    end
    
    sumPromCorrMat = corr(sumPromRep,'rows','pairwise');
    mer7CorrMat = corr(mer7Rep,'rows','pairwise');
    combineMat = tril(sumPromCorrMat) + triu(mer7CorrMat,1);
    
    imagesc(combineMat)
    borders = [0; cumsum(accumarray(idxVec,1))]+0.5;
    tickPos = movmean(borders,2,'Endpoints','discard');
    nRepeats = accumarray(idxVec, repIdx, [], @(x)numel(unique(x(~isnan(x)))));
    
    hold on
    plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
    plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
    plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)
    %     corrSumProm = plotSumPromCorr(intStrains, checWTdelLactisSwap,0);
    %     normProfile = chromosomes2fullProfile(checWTdelLactisSwap, intStrains);
%     mers = mer_occupancy(normProfile, 7, 'method','else','intBases',intBases);
%     corr7mer = corr(mers.score);
%     combinedMat = tril(corrSumProm) + triu(corr7mer,1);
%     imagesc(combinedMat)
%     plotgrid(combinedMat)
    colormap(gca, [0.7 0.7 0.7; cMapCorr])
    caxis([-0.045 1])
    YLcorrMat = regexprep(intStrains, {'_d','.*_lactis'},{' \\Delta', '{\\itK.lac}'});
    if strcmp(TF1, 'Fkh2')
        YLcorrMat{1} = [YLcorrMat{1}, ' (2)'];
    else
        YLcorrMat{1} = [YLcorrMat{1}, ' (',num2str(nRepeats(1)), ')'];
    end
    set(gca, 'YTick', tickPos, 'YTickLabel',YLcorrMat, 'XTick',[])
    hold on
end
axes('Position', [xPos(18) 0.1 3*(Wscatter)+0.05 0.1])
caxis([-0.045 1])
colormap([0.7 0.7 0.7; cMapCorr])
cbr = colorbar('Location','southoutside')
axis off
set(cbr,'Ticks', [0:0.2:1], 'fontSize',12)
ylabel(cbr, 'Correlation', 'fontSize',15)

