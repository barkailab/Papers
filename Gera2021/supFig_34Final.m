%% Figure 3-figure supplement 1 and Figure 4-figure supplement 1
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
load('promCorrSort.mat')
if ~exist('checWTdelLactisSwap','var')
    load('checWTdelLactisSwap.mat')
end

%% Figure 3-figure supplement 1: DBDs sequence alignmnet between paralogs ("barcodes") and mean signal around In-Vitro motif
clear xPos yPos
pfamAlignFiles = dir('./ScerDBD/*.psi')
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y', '-'};
groupSize = [2,3,8,7,1];
baseColor = [255, 92, 51; 128, 191, 255; 255, 214, 51; 55, 153, 102]/255;
clear yPos xPos
conservationTH = 0.4;
AAtype = [repmat(1,1,2), repmat(2,1,3), repmat(3,1,8), repmat(4,1,7), 5];
typeCol =  [brighten(baseColor, 0.2); [194 194 184]/255; [255 255 204]/255; 1 1 1];
typeCol = brighten(typeCol , 0.3);

group{13} = 'ZCs_bZs';
group{14} = 'ZCs_bZs';
group{15} = 'ZFs';
group(1:12) = {'others'};

ySpace = 0.05;
xPos{15} = 0.08;
xPos{13} = 0.38;
xPos{14} = xPos{13};
xPos{10} = 0.68;
xPos{5} = xPos{10};
xPos{2} = xPos{10};
xPos{4} = xPos{10};
xPos{1} = xPos{10};
xPos{6} = xPos{10};
xPos{11} = xPos{10};
xPos{9} = xPos{10};
xPos{12} = xPos{10};

yPos{15} = [0.8788:-ySpace:0.1242];
yPos{13} = yPos{15}(1:6);
yPos{14} = yPos{15}(7:10);
yPos{10} = yPos{15}(3);
yPos{5} = yPos{15}([4,7]);
yPos{2} = yPos{15}(5);
yPos{4} = yPos{15}(6);
yPos{1} = yPos{15}(9);
yPos{6} = yPos{15}(2);
yPos{11} = yPos{15}(1);
yPos{9} = yPos{15}(8);
yPos{12} = yPos{15}(10);

W =  0.165;
H = 0.03;
familyIdx = find(~cellfun('isempty', yPos));
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
for i = familyIdx
    DBDalign = multialignread([pfamAlignFiles(i).folder,'/', pfamAlignFiles(i).name]);
    [~, sortPos] = ismember(extractBetween(extractfield(DBDalign, 'Header'),'_','_'), upper(promCorrSort.(group{i}).TFsSort'));
    if numel(sortPos) ==0
            [~, sortPos] = ismember(extractAfter(extractfield(DBDalign, 'Header'),'_'), upper(promCorrSort.(group{i}).TFsSort'));
    end
    DBDalign(sortPos == 0) = [];
    [~,sortPos] = sort(sortPos(sortPos>0));
    DBDalign = DBDalign(sortPos);
    
       clear alignMat
    for z = 1: size(DBDalign,2)
        [~,idx] = ismember(upper(DBDalign(z).Sequence)', AAorder);
        alignMat(z, :) = idx';
    end
    nPara = size(alignMat,1)/2;
    colorMap = alignMat;
    currTable = detectImportOptions('similarRegression.xlsx','sheet', extractBefore(pfamAlignFiles(i).name, '.psi'));
    currTable =  setvaropts(currTable,'weights', 'Type',  'double');
    currTable = readtable('similarRegression.xlsx', currTable);
    importanceTH = 1.5*mean(currTable{:,3},'omitnan');
    conservationTH = max(currTable{:,1})*0.5;

    colorMap = AAtype(colorMap);
%     colorMap(:, currTable.weights>=importanceTH) = colorMap(:, currTable.weights>=importanceTH)+5;
    colorMap(:,~(currTable.conservation>=conservationTH| currTable.weights>=importanceTH)) = 6;
    colorMap(:,currTable.conservation>=conservationTH & ~(currTable.weights>=importanceTH)) = 5;
    colorMap(alignMat == 21) = 7;

    for n = 1:nPara
        axes('Position', [xPos{i}, yPos{i}(n), W, H])
        imagesc(colorMap(n*2-1:n*2,:));
        %plotgrid(colorMap(n*2-1:n*2,:))
        hold on
        plot(xlim,[1.5,1.5],'k', 'LineWidth',0.5);
        colormap(gca, typeCol);
        set(gca, 'YTick', [1:2], 'YTickLabel', extractBetween(extractfield(DBDalign(n*2-1:n*2), 'Header'),'_','_'), 'TickLength', [0.001 0.001], 'FontSize',13)
        xticks([]);
        caxis([0.5 size(typeCol,1)+0.5]);
        set(gcf, 'Color', 'w');
        changedPos = find(diff(alignMat(n*2-1:n*2,:)));
        scatter([changedPos,changedPos], [ones(1, length(changedPos)), 2*ones(1, length(changedPos))], 20, 'k.');
        missingPos = find(all(alignMat(n*2-1:n*2,:) == 21));
        if n==1
            if i == 10
                titleStr = 'SANT/Myb';
            else
                titleStr = extractBefore(pfamAlignFiles(i).name, '.');
            end
            title(titleStr, 'fontSize',10)
        end
        for l = get(gca, 'Children')'
            if strcmp(class(l), 'matlab.graphics.chart.primitive.Line')
                if diff(l.XData) == 0 & ismember(l.XData(1) -0.5, missingPos) & ismember(l.XData(1) +0.5, missingPos)
                    delete(l)
                elseif diff(l.XData) == 0
                    l.LineWidth = 0.1;
                end
            end
        end
             out = meanSignalAroundInVitroMotif(checWTdelLactisSwap, string(extractBetween(DBDalign(n*2-1).Header,'_','_')) , 'nmer',5,'windowSize',150);
             if numel(out) > 0
                 for p =  find(~cellfun('isempty',out.motif))   
                        axes('Position', [xPos{i}+W+0.002*p+(p-1)*W/4, yPos{i}(n), W/4, H])
                        plot([-150:150], out.aroundMotif{p}(1,:),'Color', [0 0 0],'LineWidth',1,'DisplayName','Mig2')
                        hold on
                        plot([-150:150], out.aroundMotif{p}(2,:),'Color', [0.7 0.7 0.7],'LineWidth',1,'DisplayName','Mig3')
                        set(gca, 'XTick', 0, 'XTickLabel', strrep(out.motif{p},'.','N'), 'YTick',[]);
                        axis tight
                 end
             end
    end
end



%% Figure 4-figure supplement 1A: histogram sumProm corr all samples vs repeats
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
swapNames = allSamples(contains(allSamples, '_DBD')& ~contains(allSamples,'Rlm1'));
[corrMat, idxVec, fullMat, idxRep] = getRepeatsCorr(swapNames, 'dataType','sumProm');
[~, firstOcc] = unique(idxRep,'stable');
redFullMat = fullMat(:, firstOcc);
allSamplesCorr = 1-squareform(1-corr(redFullMat, 'rows','pairwise'));

for i = 1:numel(swapNames)
    repeatsCorr{i} = 1-squareform(1-corr(redFullMat(:,idxVec(firstOcc)==i),'rows','pairwise'));
end
repeatsCorr = cat(2,repeatsCorr{:});

% histograms - promoter signal
typeCol =  [brighten(baseColor, 0.2); [194 194 184]/255; [255 255 204]/255; 1 1 1];

col = typeCol([end-1, end-2],:);
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
axes('Position', [0.1 0.5 0.3 0.25])
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

%% Figure 4-figure supplement 1C: binding signal corr matrix for each pair, including repeats 
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
swapNames = allSamples(contains(allSamples, '_DBD')& ~contains(allSamples,{'Rlm1'}));
goodPos = createIntBasesForMotif();
swapPairsIdx = find(any(ismember([summaryTable.p1,summaryTable.p2], extractBefore(swapNames,'_')),2));
colMap = brighten(flipud(bone),0.2);

W = 0.055;
H = W*1.8;
xspacer = 0.065;
yspacer = 0.03;
xPos = repmat([0.1+(W+xspacer).*[0,1,2,3,4]],1,2);
yPos = repmat(0.6-(H+yspacer).*[0,1],5,1);
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
c=1;
for i = swapPairsIdx'
    axes('Position', [xPos(c), yPos(c), W, H])
    intStrains = {summaryTable.p1{i}, swapNames{startsWith(swapNames, summaryTable.p1{i})},...
        summaryTable.p2{i}, swapNames{startsWith(swapNames, summaryTable.p2{i})}} ;
    [currSumPromMat, idxVec, ~, repIdx] = getRepeatsCorr(intStrains, 'dataType','sumProm');
    [curr7merMat, ~, ~, ~] = getRepeatsCorr(intStrains, 'dataType','7mer');
    currCombinedMat = tril(currSumPromMat,-1) + triu(curr7merMat);
    imagesc(currCombinedMat)
    title([summaryTable.p1{i}, ', ', summaryTable.p2{i}])
    borders = [0; cumsum(accumarray(idxVec,1))]+0.5;
    tickPos = movmean(borders,2,'Endpoints','discard');
    nRepeats = accumarray(idxVec, repIdx, [], @(x)numel(unique(x(~isnan(x)))));
    yLabels = strcat(strrep(intStrains,'_',' ')', ' (',num2str(nRepeats), ')');
    yLabels = strrep(yLabels, ' DBD', '_{DBD}');
    set(gca, 'YTick', tickPos, 'YTicklabel', yLabels,'XTick',[], 'TickLength',[0 0])
    hold on
    plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
    plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
    plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)
    colormap(gca, colMap)
    %cbr = colorbar()
    caxis([0.4 1])
    c=c+1; 
end

axes('Position', [xPos(c-1)+W+0.01, yPos(c-1), 0.02, H])
axis off
caxis([0.4 1])
cbr = colorbar()
colormap(gca, colMap)
set(cbr,'Ticks', [0.4:0.2:1], 'FontSize',12)
ylabel(cbr, 'Correlation','FontSize',15)



%% Figure 4-figure supplement 1B: protein scheme with DBD annotations
load('DBDdefForSwapping.mat')
exampleStrains = {'Gis1','Rph1';'Dot6','Tod6'; 'Yrr1','Pdr8'}%'Fkh1', 'Fkh2'};

allSamples = fieldnames(checWTdelLactisSwap.sumProm);
swapSamples = allSamples(contains(allSamples,'DBD') & ~contains(allSamples,'Rlm1'));
swapIdx = find(contains(summaryTable.p1, [regexp(swapSamples, '^[A-Z][a-z]{2}\d+', 'match','once');regexp(swapSamples, '(?<=_)[A-Z][a-z]{2}\d+', 'match','once')]));
maxL = max(summaryTable.proteinLength(swapIdx,:),[],'all');

figure('Units','normalized','Position', [0 0 0.5 1])
ax(1) = axes('Position', [0.05 0.1 0.1 0.23], 'XLim',  [0 maxL+50]);
ax(2) = axes('Position', [0.185 0.1 0.1 0.23],'XLim',  [0 maxL+50]);
ax(3) = axes('Position', [0.32 0.1 0.1 0.23],'XLim',  [0 maxL+50]);
DBDcol = [[207 213 232]; [235 210 208]]/255;
DBDtextCol = brighten(DBDcol,-0.8);
c=1;
yShift = [-0.2 0.2];
for i = swapIdx'
    axes(ax(ceil(c/3)))
    if ismember(summaryTable.p1(i),exampleStrains(:,1))
        pOrder = [1,2];
    else
        pOrder = [2,1];
    end
    z = 1;
    for p = pOrder
        plot([1,summaryTable.proteinLength(i,p)], yShift(z)*[1 1]+c, 'color', brighten(DBDcol(z,:),-0.7), 'LineWidth',3)
        text(summaryTable.proteinLength(i,p)+4, yShift(z)+c, num2str(summaryTable.proteinLength(i,p)),...
            'HorizontalAlignment','left', 'fontSize',12)
        hold on
        currDBD = DBDdefForSwapping.(upper(summaryTable.(sprintf('p%d',p)){i})).AA;
        rectangle('Position',[currDBD.start{1}, c+yShift(z)-0.1, currDBD.end{1}-currDBD.start{1}, 0.2], 'Curvature',0.2, 'FaceColor',DBDcol(z,:), 'EdgeColor', brighten(DBDcol(z,:),-0.5))
        text(-4, c+yShift(z), summaryTable.(sprintf('p%d',p)){i},'fontSize',12, 'HorizontalAlignment','right')
            text(currDBD.start{1}, c+yShift(z)-0.1, num2str(currDBD.start{1}), 'HorizontalAlignment','right', 'VerticalAlignment', 'baseline','fontSize',10, 'color',DBDtextCol(z,:))
            text(currDBD.end{1}, c+yShift(z)-0.1, num2str(currDBD.end{1}), 'HorizontalAlignment','left', 'VerticalAlignment', 'baseline','fontSize',10, 'color',DBDtextCol(z,:))
        if c == 1
            %text(mean([currDBD.start{1},currDBD.end{1}]), c+yShift(z)*2, 'DBD', 'color',  brighten(DBDcol(z,:),-0.5),...
                %'fontSize',12, 'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold')
        end
        z = z+1;
    end  
    c=c+1;
    axis off
end
set(gcf, 'color','w')
for i = 1:3
    set(ax(i), 'YDir','reverse')
end




