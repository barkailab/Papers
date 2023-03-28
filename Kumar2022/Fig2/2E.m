descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
 myCurrTable(ismember(myCurrTable.TF, {'PIP2';'ARO80';'CAT8';'CHA4';'PHD1';'MSN1';'GCR1'}),:)= [];
 myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};
 
 subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll);
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumProm.(fldNames{i}) (subtelomereGenes) = nan;
end
clearvars fldNames i


 
for TFNo = 1:length(myCurrTable.TF)
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    WTDBD.name(TFNo,1) = {TFReq};
    WTDBD.family(TFNo,1) =  myCurrTable.FAMILY(TFNo);
    x =corr([medianSumProm.(TFReq) medianSumProm.(DBDName)], 'rows', 'pairwise');
    WTDBD.corr(TFNo,1) = x(2); clearvars x;
    WTDBD.intCorrDBD(TFNo, 1)= intCorrSumPromNewWOSubTel.(DBDName);
     WTDBD.intCorrWT(TFNo, 1)= intCorrSumPromNewWOSubTel.(TFReq);
    WTDBD.DBDtopSize(TFNo, 1)= sum(zscore(medianSumProm.(DBDName)(~isnan(medianSumProm.(DBDName)))) >2.5);
     WTDBD.WTtopSize(TFNo, 1)= sum(zscore(medianSumProm.(TFReq)(~isnan(medianSumProm.(TFReq)))) >2.5);
   WTDBD.intersectSize(TFNo, 1)= length(intersect(find(zscore(medianSumProm.(DBDName)(~isnan(medianSumProm.(DBDName)))) >2.5), find(zscore(medianSumProm.(TFReq)(~isnan(medianSumProm.(TFReq)))) >2.5))) ;
     
    clearvars  TFReq nonDBDName
end
WTDBD = struct2table(WTDBD);
 WTDBD = sortrows(WTDBD,'corr','descend');
%%
DBDcolor = [0.2130    0.6299    0.3328];
WTcolor = [0.2062    0.5313    0.7562];

figure(  'Renderer', 'painters', 'Color', [1 1 1]);
axes('Position',  [0.1300 0.600 0.7750 0.3150])
myBar =bar(WTDBD.corr, 'FaceColor', 'flat', 'BarWidth', 1);
xlim([0.5 length(myCurrTable.TF)+0.5])
xticks(1:length(myCurrTable.TF)); xticklabels(WTDBD.name)
% hold on; plot((1:length(myCurrTable.TF)) , WTDBD.intCorrDBD)
hold on; scatter((1:length(myCurrTable.TF)) , WTDBD.intCorrDBD, 'filled', 'MarkerFaceColor', DBDcolor)
hold on; scatter((1:length(myCurrTable.TF)) , WTDBD.intCorrWT, 'filled', 'MarkerFaceColor', WTcolor)
ylabel('DBD-WT corr');

axes('Position',  [0.1300 0.10 0.7750 0.2])
bar(WTDBD.DBDtopSize, 'EdgeColor', 'none', 'FaceColor', DBDcolor, 'FaceAlpha', 1);
% set(gca, 'XColor', 'none', 'YColor', 'none'); xticks([]); yticks([])
set(gca, 'YDir', 'reverse')
hold on; bar((WTDBD.intersectSize/2),  'EdgeColor', 'none', 'FaceColor', [0.9987    0.8134    0.3874], 'FaceAlpha', 1)
ylim([0 200])
xlim([0.5 length(myCurrTable.TF)+0.5])

axes('Position',  [0.1300 0.30 0.7750 0.2]);
bar(WTDBD.WTtopSize, 'EdgeColor', 'none', 'FaceColor', WTcolor , 'FaceAlpha', 1);
% set(gca, 'XColor', 'none', 'YColor', 'none')
xticks([]);
hold on; bar((WTDBD.intersectSize/2),  'EdgeColor', 'none', 'FaceColor', [0.9987    0.8134    0.3874], 'FaceAlpha', 1)
ylim([0 200])
xlim([0.5 length(myCurrTable.TF)+0.5])
% max(WTDBD.DBDtopSize)
% max(WTDBD.WTtopSize)