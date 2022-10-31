%% Figures 3 and 4
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
load('SC_genome.mat')
SC_genomeCat = upper(cat(2,SC_genome.Sequence));
if ~exist('checWTdelLactisSwap','var')
    load('checWTdelLactisSwap.mat')
end
% load('goodAllPWM.mat')

%% Figure 3B - horizontal bar plots of AA changes 
summaryTable=summaryTable([1:10,21:30,17:20,11:16],:)
YL = summaryTable.label;
PairsWithNoSR = {'Spt23','Smp1','Vhr1'};
intBasesMotif = createIntBasesForMotif();
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
swapSamples = allSamples(contains(allSamples, 'DBD') & ~contains(allSamples, 'Rlm1'));
baseColor = [255, 92, 51; 128, 191, 255; 255, 214, 51; 55, 153, 102]/255;
typeCol =  [brighten(baseColor, 0.2); [194 194 184]/255; [255 255 204]/255; 1 1 1];
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y', '-'};
typeCol = brighten(typeCol , 0.3);
clearvars changeTable
% calcualte fractions
allPsiFiles=dir('./ScerDBD/*.psi')
for i=1:numel(allPsiFiles)
    try
     dbdInPsi{i}=extractBetween(extractfield(multialignread(['./ScerDBD/',allPsiFiles(i).name]),'Header'),'_','_')
    end
end
[~,psiIdx]=ismember(upper(summaryTable{:,[1,2]}),cat(2,dbdInPsi{:}))
psiIdx=reshape(sum(psiIdx(:)>cumsum(cellfun('prodofsize',dbdInPsi)),2)+1,30,2);
psiIdx(contains(summaryTable.p1,'Dot6'),:)=find(contains(extractfield(allPsiFiles,'name'),'TD6.psi'))

for i = 1: size(summaryTable,1)
    % barcode
    DBDalign = multialignread(['./ScerDBD/',allPsiFiles(psiIdx(i,1)).name]);
    
    [~, sortPos] = ismember(upper(summaryTable{i,[1:2]}),extractBetween(extractfield(DBDalign, 'Header'),'_','_'));
    DBDalign = DBDalign(sortPos);
    clear alignMat
    for z = 1: size(DBDalign,2)
        [~,idx] = ismember(upper(DBDalign(z).Sequence)', AAorder);
        alignMat(z, :) = idx';
    end
    gapRes = all(alignMat == 21,1);
    alignMat = alignMat(:,~gapRes);
    changeTable.allChg(i,1)=sum(diff(alignMat)~=0);
    changeTable.allFrac(i,1)=mean(diff(alignMat)~=0);
    try      
        currTable = detectImportOptions('similarRegression.xlsx','sheet', extractBefore(allPsiFiles(psiIdx(i,1)).name, '.psi'));
        currTable =  setvaropts(currTable,'weights', 'Type',  'double');
        currTable = readtable('similarRegression.xlsx', currTable);
        specificTH = 1.5*mean(currTable{:,3}, 'omitnan');
        conservationTH = max(currTable{:,1})*0.5;
        
        currTable = currTable(~gapRes,:);

        changeTable.impChg(i,1)=sum(diff(alignMat(:,currTable{:,3}>=specificTH))~=0);
        changeTable.impFrac(i,1)=mean(diff(alignMat(:,currTable{:,3}>=specificTH))~=0);
    catch
        noSR(i)=true;
        changeTable.impChg(i,1)=NaN;
        changeTable.impFrac(i,1)=NaN;
    end
    try
        [crMat, idxVec, ~, repIdx] = getRepeatsCorr(summaryTable{i,[1,2]},'dataType','sumProm');
        [repIdx,uRpt]=unique(repIdx,'stable');
        idxVec=idxVec(uRpt);
        crMat=crMat(uRpt,uRpt);
        meanCrProm(i,1)=mean(crMat(idxVec==1,idxVec==2),'all');
        stdProm(i,1)=std(crMat(idxVec==1,idxVec==2),[],'all');
        
        [crMat, idxVec, ~, repIdx] = getRepeatsCorr(summaryTable{i,[1,2]},'dataType','7mer');
        [repIdx,uRpt]=unique(repIdx,'stable');
        idxVec=idxVec(uRpt);;
        crMat=crMat(uRpt,uRpt);
        merCorr(i,1)=mean(crMat(idxVec==1,idxVec==2),'all');
        stdMer(i,1)=std(crMat(idxVec==1,idxVec==2),[],'all');
    catch
        meanCrProm(i,1)=summaryTable.sumPromCorr(i);
        stdProm(i,1)=NaN;%std(crMat(idxVec==1,idxVec==2),[],'all');
        fullProfile = chromosomes2fullProfile(checWTdelLactisSwap,summaryTable{i,1:2});
        motifs = mer_occupancy(fullProfile, 7,'window',31, 'intBases', intBasesMotif, 'method','else');
        merCorr(i,1)=diag(corr(motifs.score,'rows','pairwise'),1);

        stdMer(i,1)=NaN;%std(crMat(idxVec==1,idxVec==2),[],'all');
    end
end

changeTable=struct2table(changeTable)
familyNames = {'Zinc finger', 'others','bZIP', 'Zinc cluster'};
familyDiv=[10,20.5,25]+0.75;
figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
axes('Position', [0.2 0.2 0.1 0.7])
barCol = [cbrewer2('Pastel1',2) ;[1,1,1]];
barCol(1,:) = typeCol(end-2,:);
barPos = [1:10, 11.5:20.5, 22:25, 26.5:31.5];
% vertical bars 3Plots
close all
figure('Units','normalized','Position',[0.5698 0.2324 0.3344 0.5019], 'color','w')
subplot(2,1,1)
bar(barPos, changeTable{:,[4]}*100,'FaceColor',barCol(2,:))
set(gca,'Ylim',[0 60],'Xlim',quantile(barPos,[0 1])+[-.5 .5])
hold on
plot( [1;1].*familyDiv,repmat(ylim',1,3), 'color',[1 1 1]*0.3, 'LineStyle', '--','LineWidt',1.5)
text(barPos([1,11,21,25]),repmat(0.99.*max(ylim),4,1), familyNames, 'fontSize',12, 'HorizontalAlignment','left','VerticalAlignment','top')
ylabel('mut. fraction imp AA')
set(gca, 'XTick', barPos, 'XTicklabel',YL, 'fontSize',10,'XAxisLocation','top','XTickLabelRotation',90)
subplot(6,1,4)
bar(barPos, changeTable{:,[2]}*100,'FaceColor',barCol(1,:))
set(gca,'Ylim',[0 60],'Xlim',quantile(barPos,[0 1])+[-.5 .5])
hold on
plot( [1;1].*familyDiv,repmat(ylim',1,3), 'color',[1 1 1]*0.3, 'LineStyle', '--','LineWidt',1.5)
set(gca, 'XTick', [], 'fontSize',10)
ylabel('mut. fraction all AA')
subplot(9,1,[7]+1)
crType=sum(meanCrProm>=[0,.4,.81],2)
cMapBars = flipud([126,8,119; 138,139,192; 183,205,227]/255);
hold off
for i=1:max(crType)
    bar(barPos-0.2,meanCrProm.*(crType==i),'BarWidth',.4,'FaceColor',cMapBars(i,:),'EdgeColor',[0 0 0])
    hold on
    bar(barPos+0.2,merCorr.*(crType==i),'BarWidth',.4,'FaceColor',cMapBars(i,:),'EdgeColor',[1 1 1].*0.5)
end
set(gca,'Ylim',[0 1],'Xlim',quantile(barPos,[0 1])+[-.5 .5])
plot( [1;1].*familyDiv,repmat(ylim',1,3), 'color',[1 1 1]*0.3, 'LineStyle', '--','LineWidt',1.5)
set(gca, 'XTick', [], 'TickLength',[0 0], 'fontSize',10,'XAxisLocation','bottom')
ylabel('binding divergence')
errorbar(barPos-0.2,meanCrProm,stdProm,'vertical','LineStyle','none','Color','k')
errorbar(barPos+0.2,merCorr,stdMer,'vertical','LineStyle','none','Color','k')
%save_gf(gcf,sprintf('horBars_2A'),'type',{'svg'},'paper','tamar21','size',[])


%% Figure 3A - DBD barcode examples, WTs scatters, signal arround inVitro motif (HeatMap)
exampleStrains = {'Gis1','Rph1';'Dot6','Tod6'; 'Yrr1','Pdr8'};
intPsiFiles = {'zf-C2H2.psi';  'TD6.psi'; 'Znclus.psi'}%'Forkhead.psi'};
familyNameTitle = {'C2H2 Zinc finger', 'SANT/myb', 'Zinc cluster'};
HLproms={[2370,2209,4765,4785,759,4718,5543,28],[3317 4498 1319 2113 1617 1865],[2804,3864   2344 933 6106]}
load('promoterLengthsORF.mat')
GP=load('./group_imp.mat')
patterns = {'CCCCT|AGGGG','CATCG|CGATG', 'CGGA.AT|AT.TCCG', }%'TGTTT|AAACA'}
pfamAlignFiles = dir('./ScerDBD/*.psi')
AAorder = {'D', 'E', 'K', 'R','H', 'A', 'W', 'V', 'I', 'L', 'P', 'F', 'M', 'T', 'S', 'N', 'G', 'Q', 'C', 'Y', '-'};
groupSize = [2,3,8,7,1];
baseColor = [255, 92, 51; 128, 191, 255; 255, 214, 51; 55, 153, 102]/255;
AAtype = [repmat(1,1,2), repmat(2,1,3), repmat(3,1,8), repmat(4,1,7), 5];
%typeCol =  [brighten(baseColor, 0.2); [191 128 255]/255; [255 255 204]/255; 1 1 1];
typeCol =  [brighten(baseColor, 0.2); [194 194 184]/255; [255 255 204]/255; 1 1 1];
typeCol = brighten(typeCol , 0.3);
cMap = flipud(bone);
cMapMotif=brewermap(128,'purples')
scatterCol = [51 51 51]/255%[148 175 211]/255;
markerCol = brighten(scatterCol, -0.5);
motifScatterCol = cMapMotif(64,:);%cbrewer2('Pastel1',1);
intBasesMotif = createIntBasesForMotif();

% Hbarcode =  0.03;
% spacer = 0.03;
% xPos = [0.05];
% yPos = [0.95:-(Hbarcode+2*spacer+Hscatter+0.04):0.05] ;
% Wbarcode =  0.15;
% Hscatter = Wbarcode*1.25;
% Wscatter = Hscatter/2;

Hbarcode =  0.03;
spacer = 0.03;
xPos = [0.05];
Wbarcode =  0.15;
Hscatter = Wbarcode*1.25;
yPos = [0.95:-(Hbarcode+2*spacer+Hscatter+0.04):0.05] ;
Wscatter = Hscatter/2;

for i = [1,2,3]%: size(exampleStrains,1)
    figure('Units','normalized','Position', [0.0523 0.3602 0.3386 0.4037], 'color','w')
    % barcode
    DBDalign = multialignread(['./ScerDBD/',intPsiFiles{i}]);
    [~, sortPos] = ismember(extractBetween(extractfield(DBDalign, 'Header'),'_','_'), upper(exampleStrains(i,:)));
    DBDalign(sortPos == 0) = [];
    [~,sortPos] = sort(sortPos(sortPos>0));
    DBDalign = DBDalign(sortPos);
    
    clear alignMat
    for z = 1: size(DBDalign,2)
        [~,idx] = ismember(upper(DBDalign(z).Sequence)', AAorder);
        alignMat(z, :) = idx';
    end
    gapRes = all(alignMat == 21,1);
    alignMat = alignMat(:,~gapRes);
    repType = (diff(alignMat)~=0)+(diff(AAtype(alignMat))~=0); % 0: same AA, 1: same Type 2:different Type
    colorMap = alignMat;
    currTable = detectImportOptions('similarRegression.xlsx','sheet', extractBefore(intPsiFiles{i}, '.psi'));
    currTable =  setvaropts(currTable,'weights', 'Type',  'double');
    currTable = readtable('similarRegression.xlsx', currTable);
    specificTH = 1.5*mean(currTable{:,3}, 'omitnan');
    conservationTH = max(currTable{:,1})*0.5;
    
    currTable = currTable(~gapRes,:);
%     importanceTH = min(maxk(currTable{:,3} ,10));
%     conservationTH = min(maxk(currTable{:,1} ,10));


    colorMap = AAtype(colorMap);
    colorMap(:,~(currTable.conservation>=conservationTH| currTable.weights>=specificTH)) = 6;
    colorMap(:,currTable.conservation>=conservationTH & ~(currTable.weights>=specificTH)) = 5;
    colorMap(alignMat == 21) = 7;
    
    %axes('Position', [xPos yPos(i) Wbarcode Hbarcode])
    subplot(5,3,[2,3])
    imagesc(colorMap);
    hold on
    plot(0.5+[1:size(colorMap,2)-1].*[1;1],repmat(ylim()',1,size(colorMap,2)-1),'-','Color',0*[1 1 1]) 
    plot(xlim,[1.5,1.5],'k', 'LineWidth',0.5);
    colormap(gca, typeCol);
    Lalignment = strcat(exampleStrains(i,:)',{'\fontsize{10}(P1)', '\fontsize{10}(P2)'}');
    set(gca, 'YTick', [1:2], 'YTickLabel', Lalignment, 'TickLength', [0.001 0.001], 'FontSize',13)
    xticks([]);
    caxis([0.5 size(typeCol,1)+0.5]);
    changedPos = find(repType==1);
    scatter([changedPos,changedPos], [ones(1, length(changedPos)), 2*ones(1, length(changedPos))], 80,[.5 .5 .5],'.');
    changedPos = find(repType==2);
    scatter([changedPos,changedPos], [ones(1, length(changedPos)), 2*ones(1, length(changedPos))], 80,0.*[1 1 1], '.');
    
%     for l = get(gca, 'Children')'
%         if strcmp(class(l), 'matlab.graphics.chart.primitive.Line')
%             if diff(l.XData) == 0 & ismember(l.XData(1) -0.5, missingPos) & ismember(l.XData(1) +0.5, missingPos)
%                 delete(l)
%             elseif diff(l.XData) == 0
%                 l.LineWidth = 0.1;
%             end
%         end
%     end
    text(1, 0.5, familyNameTitle{i}, 'HorizontalAlignment','left', 'fontSize',12, 'VerticalAlignment','bottom')

    % sumProm Scatter
    occByGenome = regexp(SC_genomeCat, patterns{i});
    disOccOrf = occByGenome-[GP.gene_infoR64.position(:,2)+GP.chrIdx(min(GP.gene_infoR64.position(:,1),18))];
    disOccOrf = disOccOrf.*GP.gene_infoR64.dir;
    occInGenes = sum(disOccOrf<0 & disOccOrf> -promoterLengthsORF,2);
    
    %axes('Position', [xPos-0.01 yPos(i)-spacer-Hscatter Wscatter  Hscatter])
    subplot(5,3,[4:3:13])
    maxX = max(checWTdelLactisSwap.sumProm.(exampleStrains{i,1}));
    maxY = max(checWTdelLactisSwap.sumProm.(exampleStrains{i,2}));
    scatter(checWTdelLactisSwap.sumProm.(exampleStrains{i,1})/maxX, checWTdelLactisSwap.sumProm.(exampleStrains{i,2})/maxY,...
        50 , scatterCol, 'filled', 'MarkerEdgeColor',markerCol,'MarkerFaceAlpha',.5)
        setAxisExponent()
    ylim([0 1.1])
    xlim([0 1.1])
    set(gca, 'XTick',[0:0.5:1], 'YTick',[0:0.5:1])
    text(0.05,1.1,sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(exampleStrains{i,1}), checWTdelLactisSwap.sumProm.(exampleStrains{i,2}), 'rows','pairwise')), 'fontSize',12, 'HorizontalAlignment','left','VerticalAlignment','top')
    set(gca, 'fontSize',12)
    xlabel(exampleStrains{i,1}, 'fontSize',15)
    ylabel(exampleStrains{i,2}, 'fontSize',15)
    title('Promoter binding', 'fontSize',12, 'FontWeight','normal')
    text(checWTdelLactisSwap.sumProm.(exampleStrains{i,1})(HLproms{i})/maxX, checWTdelLactisSwap.sumProm.(exampleStrains{i,2})(HLproms{i})/maxY,GP.gene_infoR64.nameNew(HLproms{i}))
    
    % signal arround motif
    [out,signalMats] = meanSignalAroundPattern(checWTdelLactisSwap, exampleStrains(i,:), patterns(i),'windowSize',150,'normalize',false);
    %axes('Position', [xPos+Wscatter+spacer/2 yPos(i)-spacer-0.4*Hscatter Wbarcode-Wscatter-spacer/2  Hscatter*0.45])
    subplot(5,6,3+6)
    plot([-150:150], out.aroundMotif{1}(1,:),'Color', [0 0 0],'LineWidth',1.5,'DisplayName',exampleStrains{i,1})
    set(gca, 'XTick', [-100 0 100], 'XTickLabel', {'-100',sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s',brighten(motifScatterCol,-0.5), extractBefore(strrep(patterns{i},'.','N'),'|')),'100'},'XtickLabelRotation',0,'Xtick',[], 'YTick',floor(ylim),'fontSize',10,'ylim',[0 max(out.aroundMotif{1},[],'all')*1.1],'xlim',[-150 150]);
    %ylabel('mean signal', 'fontWeight','normal', 'fontSize',11)
    
    subplot(5,6,4+6)
    plot([-150:150], out.aroundMotif{1}(2,:),'Color', [0.7 0.7 0.7],'LineWidth',2,'DisplayName',exampleStrains{i,2})
    set(gca, 'XTick', [-100 0 100], 'XTickLabel', {'-100',sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s',brighten(motifScatterCol,-0.5), extractBefore(strrep(patterns{i},'.','N'),'|')),'100'},'XtickLabelRotation',0,'Xtick',[], 'YTick',[],'fontSize',10,'ylim',[0 max(out.aroundMotif{1},[],'all')*1.1],'xlim',[-150 150]);
    %ylabel('mean signal', 'fontWeight','normal', 'fontSize',11)
    
    % try nicer clustering
    [cluIdx,cluCen]=kmeans(zscore(log2([sum(signalMats{1}(:,101:201),2), sum(signalMats{2}(:,101:201),2)]+1)),3,'MaxIter',500,'Replicates',100,'Options',statset('UseParallel',true))
    cluIdx=changem(cluIdx,4-ranking(cluCen(:,1)),1:3)
    [~,sIdx]=sortrows(table(cluIdx,sum(signalMats{1}(:,101:201),2)),[1,2],'descend')
    
    %[~,sIdx]=sort(sum(signalMats{1}(:,101:201),2),'descend');
    subplot(5,6,[15:6:30])
    imagesc(movmean(signalMats{1}(sIdx,:),3,2),'XData',[-150 150],[0 20])
    xlabel('')
    ylabel(sprintf('n=%d',size(signalMats{1},1)))
    set(gca,'Ytick',[])    
    colormap(gca,cMapMotif)
    set(gca, 'XTick', [-100 0 100], 'XTickLabel', {'-100',sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s',brighten(motifScatterCol,-0.5), extractBefore(strrep(patterns{i},'.','N'),'|')),'100'},'XtickLabelRotation',0, 'YTick',[ ],'fontSize',10,'xlim',[-150 150]);
 %   title(exampleStrains{i,1})
    subplot(5,6,[16:6:30])
    imagesc(movmean(signalMats{2}(sIdx,:),3,2),'XData',[-150 150],[0 20])
    set(gca, 'XTick', [-100 0 100], 'XTickLabel', {'-100',sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s',brighten(motifScatterCol,-0.5), extractBefore(strrep(patterns{i},'.','N'),'|')),'100'},'XtickLabelRotation',0, 'YTick',[ ],'fontSize',10,'xlim',[-150 150]);
    colormap(gca,cMapMotif)
  %  title(exampleStrains{i,2})
    
%     cbr = colorbar('Location', 'east')
%     set(cbr, 'Position',  [0.1884-xPos(1)+xPos(i) yPos-spacer-Hscatter 0.0090 0.1199], 'AxisLocation', 'out')
%     ylabel(cbr, '# motifs', 'fontSize',12)
%     colormap(gca, cMap)
%     caxis([0 5])

    % 7mer scatter
    %axes('Position', [xPos+Wscatter+spacer/2 yPos(i)-2*spacer-0.4*Hscatter-Hscatter*0.45 Wbarcode-Wscatter-spacer/2  Hscatter*0.41])
    subplot(5,3,[6:3:15])
    fullProfile = chromosomes2fullProfile(checWTdelLactisSwap,exampleStrains(i,:));
    motifs = mer_occupancy(fullProfile, 7,'window',31, 'intBases', intBasesMotif, 'method','else');
    motifs.score=motifs.score./max(motifs.score);
    motifsMers = table2struct(motifs.mers);
    mIdx = find(~cellfun('isempty', regexp({motifsMers.seq}, patterns{i},'once')));
    scatter(motifs.score(:,1), motifs.score(:,2), 50 , scatterCol, 'filled', 'MarkerEdgeColor',markerCol,'MarkerFaceAlpha',.5)
    hold on
    scatter(motifs.score(mIdx,1), motifs.score(mIdx,2), 50 , motifScatterCol, 'filled', 'MarkerEdgeColor',markerCol)
    %set(gca, 'XTick', xlim, 'XTickLabel', 'YTick',ylim, 'fontSize',8)
    set(gca, 'fontSize',12)
    xlabel(exampleStrains{i,1}, 'fontSize',15)
    ylabel(exampleStrains{i,2}, 'fontSize',15)
    xlim([0 1.05])
    ylim([0 1.05])

    title('7mer binding', 'fontSize',12, 'FontWeight','normal')
    title(sprintf('%.2f',corr(motifs.score(:,1), motifs.score(:,2))))
end


%% Figure 4B - binding signal corrMatrix
exampleStrains = {'Gis1','Rph1';'Dot6','Tod6'; 'Yrr1','Pdr8'}%'Fkh1', 'Fkh2'};
intBases = createIntBasesForMotif();
colMap = flipud(bone);

Hscatter = 0.25/2;
Wscatter = Hscatter/2;
spacer = 0.006;
xPos = [0.1: spacer+Wscatter:1-Wscatter];
yPos = 0.7;

figure
for i = 1: size(exampleStrains,1)
    WT1 = exampleStrains{i,1};
    WT2 = exampleStrains{i,2};
    s1 = [WT1,'_',WT2,'_DBD'];
    s2 = [WT2,'_',WT1,'_DBD'];
    intStrains = {WT1,s1,WT2,s2};
        corrSumProm = plotSumPromCorr(intStrains, checWTdelLactisSwap, 0);
        normProfile = chromosomes2fullProfile(checWTdelLactisSwap, intStrains);
        mers = mer_occupancy(normProfile, 7, 'method','else','intBases',intBases);
        corr7mer = corr(mers.score);
%         [corrSumProm, idxSumProm] = getRepeatsCorr(intStrains, 'dataType','sumProm');
%         [corr7mer, ~] = getRepeatsCorr(intStrains, 'dataType','7mer');

    axes('Position', [xPos(i) yPos Wscatter Hscatter])
    combineMat = tril(corrSumProm) + triu(corr7mer,1);
    imagesc(combineMat)
    %set(gca, 'YTick', [1:numel(intStrains)], 'YTicklabel', strrep(intStrains,'_',' '),'XTick',[])
    %borders = [0; cumsum(accumarray(idxSumProm,1))]+0.5;
    %tickPos = movmean(borders,2,'Endpoints','discard');
    plotgrid(combineMat)
    set(gca, 'YTick', [], 'YTicklabel', strrep(intStrains,'_',' '),'XTick',[])
    hold on
%     plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
%     plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
    plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)
    colormap(gca, colMap)
    caxis([0.4 1])
%     if i==1
%         set(gca, 'YTick', tickPos, 'YTickLabel', {'P1','P1 P2_{DBD}','P2','P2 P1_{DBD}'},'fontSize',10)
%     end
    title([WT1,', ',WT2], 'fontSize',15)
    if i == size(exampleStrains,1)
        cbr = colorbar('Location', 'east')
        set(cbr, 'Position', [xPos(size(exampleStrains,1))+Wscatter+0.005 yPos 0.0075 Hscatter])
        set(cbr, 'AxisLocation', 'out', 'Ticks', [0.4:0.2:1],'fontSize',12)
        ylabel(cbr, 'Correlation', 'fontSize',15)
    end
end
set(gcf, 'color','w')

saveas(gcf, 'Fig2corrMatExamples.svg')



%%  Figure 4C- DBD-swapping summary scatter 
allSamples = fieldnames(checWTdelLactisSwap.sumProm);
swapSamples = allSamples(contains(allSamples, 'DBD') & ~contains(allSamples, 'Rlm1'));
exampleStrains = {'Gis1','Rph1';'Dot6','Tod6'; 'Yrr1','Pdr8'}%'Fkh1', 'Fkh2'};
exampleStrainsIdx = find(contains(swapSamples,exampleStrains));
capSize=3;
zcColor=[151,178,213]/255
for s = 1:numel(swapSamples)
    donor = regexp(swapSamples{s}, '(?<=_).*(?=_DBD)', 'match');
    acceptor = regexp(swapSamples{s}, '^.*?(?=_)', 'match');
    paralogIdx(s) = find(contains(summaryTable.p1 , [donor,acceptor]));     
    [crMat, idxVec, ~, repIdx] = getRepeatsCorr([swapSamples(s),acceptor,donor],'dataType','sumProm');
    [~,uRpt]=unique(repIdx,'stable');
    idxVec=idxVec(uRpt);
    crMat=crMat(uRpt,uRpt);
    
    sumPromCorrSwapWT(s) = mean(crMat(idxVec==1,idxVec==2),'all');%corr(checWTdelLactisSwap.sumProm.(acceptor{1}), checWTdelLactisSwap.sumProm.(swapSamples{s}), 'rows','pairwise');
    accStd(s)=std(crMat(idxVec==1,idxVec==2),[],'all');
    
    sumPromCorrSwapParalog(s) = mean(crMat(idxVec==1,idxVec==3),'all');%corr(checWTdelLactisSwap.sumProm.(donor{1}), checWTdelLactisSwap.sumProm.(swapSamples{s}), 'rows','pairwise');
    donStd(s)=std(crMat(idxVec==1,idxVec==3),[],'all');
end

figure
cMap = brighten(flipud(bone),0.4);
cMapZC =  cbrewer2('OrRd',100);
cMapJ = [cMap;cMapZC];

axes('Position', [0.02+0.38, 0.2, 0.7, 0.7])
% scatter(sumPromCorrSwapParalog(exampleStrainsIdx), sumPromCorrSwapWT(exampleStrainsIdx), 120,...
%     summaryTable.sumPromCorr(paralogIdx(exampleStrainsIdx)),'MarkerEdgeColor','k', 'LineWidth', 2)
hold on
s1 = scatter(sumPromCorrSwapParalog, sumPromCorrSwapWT, 120, summaryTable.sumPromCorr(paralogIdx),...
    'filled', 'MarkerEdgeColor','k', 'DisplayName','rest')
errorbar(sumPromCorrSwapParalog, sumPromCorrSwapWT,donStd,donStd,accStd,accStd,'k-','Linestyle','none','CapSize',capSize)
colormap(gca,cMap)

hold on
selSwap = summaryTable.familyId(paralogIdx) == 2;
% s2 = scatter(sumPromCorrSwapParalog(selSwap), sumPromCorrSwapWT(selSwap), 120,...
%     cMapZC(round(rescale(summaryTable.sumPromCorr(paralogIdx(selSwap)), 0.5, 100.4, 'InputMin',0.2,'InputMax',1)),:),...
%     'filled','MarkerEdgeColor','k', 'DisplayName', 'zinc cluster')
s2 = scatter(sumPromCorrSwapParalog(selSwap), sumPromCorrSwapWT(selSwap), 120,...
    zcColor,'LineWidth',2, 'DisplayName', 'zinc cluster')
selSwap = sumPromCorrSwapWT<0.85;
% text(sumPromCorrSwapParalog(selSwap)-range(xlim)/50, sumPromCorrSwapWT(selSwap),...
%     regexp(swapSamples(selSwap), '^...\d+','match','once'), 'fontSize',14,'HorizontalAlignment','right')
tStrains = unique([find(selSwap), exampleStrainsIdx']);
text(sumPromCorrSwapParalog(tStrains)-range(xlim)/50, sumPromCorrSwapWT(tStrains),...
    regexp(swapSamples(tStrains), '^...\d+','match','once'), 'fontSize',15,'HorizontalAlignment','right')

axis square
xlim([0.2 1]); ylim([0.2 1])
set(gca, 'XTick',[0:0.2:1], 'YTick',[0:0.2:1], 'fontSize',18)
caxis([0.2 1])
colormap(gca,cMap)
xlabel('Correlation with DBD-donor (d1,d2)',  'fontSize',20)
ylabel('Correlation with DBD-acceptor (a1,a2)' ,'fontSize',20)
hold on
plot(xlim,xlim,'--k')
%legend([s2,s1],'Location','southeast', 'fontSize',15);

cbr = colorbar('Location','eastoutside')
set(cbr, 'YTick',[0.2:0.2:1], 'AxisLocation','out')
ylabel(cbr, 'Correlation between the WTs','fontSize',20 )
box on
set(gcf, 'color','w')
title('Promoter binding')

% Zinc cluster colormap
% axes('position', [0.5935 0.833 0.313 0.01])
% imagesc(0:0.01:1)
% caxis([0.2 1])
% colormap(gca, cMapZC)
% set(gca, 'XTick',[],'YTick','')

%% protein scheme with DBD annotations
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
        %     text(currDBD.start{1}, c+yShift(z)-0.1, num2str(currDBD.start{1}), 'HorizontalAlignment','right', 'VerticalAlignment', 'baseline','fontSize',10, 'color',DBDtextCol(z,:))
        %     text(currDBD.end{1}, c+yShift(z)-0.1, num2str(currDBD.end{1}), 'HorizontalAlignment','left', 'VerticalAlignment', 'baseline','fontSize',10, 'color',DBDtextCol(z,:))
        if c == 1
            text(mean([currDBD.start{1},currDBD.end{1}]), c+yShift(z)*2, 'DBD', 'color',  brighten(DBDcol(z,:),-0.5),...
                'fontSize',12, 'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold')
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

saveas(gcf, 'Fig2DBDconstructs.svg')

