%needs pughFeatures, hahnFeatures, contOPNDPN

subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll);
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumProm.(fldNames{i}) (subtelomereGenes) = nan;
end
clearvars fldNames i

descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');


myWTTable=descTF;
myWTTable( ismember(myWTTable.TF, {'HAP4';'GCR1'}),:)=[]; %wtpoor

myDBDTable=descTF;
myDBDTable(cellfun(@isempty,myDBDTable.DBDname),:)=[];
myDBDTable( ismember(myDBDTable.TF, {'CAT8';'PIP2';'ARO80';'MSN1'}),:)=[]; %dbdpoor
myDBDTable.DBDname(ismember(myDBDTable.TF, 'YAP1'))= {'YAP1_d2_500'};

mynonDBDTable=descTF;
mynonDBDTable(cellfun(@isempty,mynonDBDTable.nonDBDName),:)=[];
mynonDBDTable(ismember(mynonDBDTable.TF, {'PHO4';'PHD1';'HAP4';'TYE7';'RTG3';'HCM1';'PUT3';'MGA1';'MOT3';'HOT1';'RTG1';'ZAP1';'ARO80'; 'MBP1';'HCM1';'INO2';'PDR3'}),:)=[]; %nondbdpoor

%%

pughClass(cellfun(@isempty, PughFeatures.FeatureClassLevel1)) = {'none'};
pughClass(~cellfun(@isempty, PughFeatures.FeatureClassLevel1)) =PughFeatures.FeatureClassLevel1(~cellfun(@isempty, PughFeatures.FeatureClassLevel1)) ;
pughClass(contains(pughClass, 'RP')) = {'RP'};
pughClass(contains(pughClass, 'STM')) = {'STM'};
pughClass(contains(pughClass, {'UNB'; 'TFO'})) = {'noTF'};
pughClass(contains(pughClass, {'05'; '06'; '07'; '08';'13';'none'})) = {'others'};
pughTypes = unique(pughClass);

hahnClass(cellfun(@isempty, HahnFeatures.coactivator)) = {'none'};
hahnClass(~cellfun(@isempty, HahnFeatures.coactivator)) =HahnFeatures.coactivator(~cellfun(@isempty, HahnFeatures.coactivator)) ;

pughClassNo(contains(pughClass, 'RP')) = 1;
pughClassNo(contains(pughClass, 'STM')) = 2;
pughClassNo(contains(pughClass, 'noTF')) = 3;
pughClassNo(contains(pughClass, 'others')) = 4;
hahnClassNo(contains(hahnClass, 'TFIID')) = 2;
hahnClassNo(contains(hahnClass, 'CR')) = 1;
hahnClassNo(contains(hahnClass, 'none')) = 3;

load('contOPNDPN')
contOPNDPN = -contOPNDPN;
contOPNDPNSelect = contOPNDPN;
contOPNDPNSelect(find(contains(pughClass, 'others')))=nan;
contOPNDPNSelect(find(contains(hahnClass, 'none')))=nan;
[~, promOrder] = sort(contOPNDPNSelect, 'descend');
% promOrder(end-sum(isnan(contOPNDPNSelect)):end) = [];
promOrder(1:sum(isnan(contOPNDPNSelect))) = [];

% myCurrTable = descTF;
% myCurrTable(cellfun(@isempty,myCurrTable.nonDBDName),:)=[];
% myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
% myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};

for TFNo=1:length(myWTTable.TF)
     TFReq = myWTTable.TF{TFNo}; 
    a(isnan(medianSumProm.(TFReq))) =nan;
    a(~isnan(medianSumProm.(TFReq))) =zscore(medianSumProm.(TFReq)(~isnan(medianSumProm.(TFReq))));
    wtMat(:,TFNo) = a(promOrder); clearvars a
end

for TFNo=1:length(myDBDTable.TF)
    DBDName = myDBDTable.DBDname{TFNo}; 
    a(isnan(medianSumProm.(DBDName))) =nan;
    a(~isnan(medianSumProm.(DBDName))) =zscore(medianSumProm.(DBDName)(~isnan(medianSumProm.(DBDName))));
    dbdMat(:,TFNo) = a(promOrder); clearvars a
end

for TFNo=1:length(mynonDBDTable.TF)
    nonDBDName = mynonDBDTable.nonDBDName{TFNo}; 
    a(isnan(medianSumProm.(nonDBDName))) =nan;
    a(~isnan(medianSumProm.(nonDBDName))) =zscore(medianSumProm.(nonDBDName)(~isnan(medianSumProm.(nonDBDName))));
    nondbdMat(:,TFNo) = a(promOrder); clearvars a
end
   
    clearvars TFReq DBDName nonDBDName


%% make heatmap figure
figure('Position', [886 42 760 953], 'Color', [1 1 1], 'Renderer', 'painters')

promAxes =axes('Position', [0.13 0.02500 0.45 0.95]);
imagesc([dbdMat wtMat nondbdMat])
colormap(promAxes, brewermap(100, 'YlGnBu'));
cb =colorbar('Location', 'northoutside'); caxis([0 2.5]); xlabel(cb, 'promoter signal zscore')
hold on; plot([length(myDBDTable.TF)+0.5 length(myDBDTable.TF)+0.5], ylim, 'k')
hold on; plot([length(myDBDTable.TF)+length(myWTTable.TF)+0.5 length(myDBDTable.TF)+length(myWTTable.TF)+0.5], ylim, 'k')
yticks([]); xticks([]);

opnAxes =axes('Position', [0.09 0.02500 0.02 0.8884]);
imagesc(contOPNDPN(promOrder))
yticks([2 sum(~isnan(contOPNDPNSelect))-2]); yticklabels(round( contOPNDPNSelect(promOrder(yticks)), 2, 'significant'))  ; xticks([]);ylabel('OPN-DPN');
colormap(opnAxes, flipud(brewermap(100, 'PuOr'))); caxis([-3.7 3])

pughAxes =axes('Position', [0.59 0.02500 0.03 0.8884]);
imagesc(pughClassNo(promOrder)'); xticks([]); yticks([])
colormap(pughAxes, [[0.9608    0.9608    0.9608] ; [0.8471    0.7020    0.3961]; [0.3529    0.7059    0.6745]]  );

hahnAxes =axes('Position', [0.63 0.02500 0.03 0.8884]);
imagesc( hahnClassNo(promOrder)'); xticks([]); yticks([])
colormap(hahnAxes, [[0.6402         0    0.4830]; [0.8457    0.7615    0.8738]] );hold on;
scatter(nan, nan, 'MarkerFaceColor',  [0.9608    0.9608    0.9608], 'DisplayName', 'RP', 'MarkerEdgeColor', 'k'); hold on;
scatter(nan, nan, 'MarkerFaceColor',  [0.8471    0.7020    0.3961], 'DisplayName', 'STM', 'MarkerEdgeColor', 'k'); hold on;
scatter(nan, nan, 'MarkerFaceColor',  [0.3529    0.7059    0.6745], 'DisplayName', 'noTF', 'MarkerEdgeColor', 'k'); hold on;

% scatter(nan, nan, 'MarkerFaceColor',  [0.6510    0.7412    0.8588], 'DisplayName', 'CR', 'MarkerEdgeColor', 'k'); hold on;
% scatter(nan, nan, 'MarkerFaceColor',  [0.2118    0.5647    0.7529], 'DisplayName', 'TFIID', 'MarkerEdgeColor', 'k'); hold on;

scatter(nan, nan, 'MarkerFaceColor',  [0.6402         0    0.4830], 'DisplayName', 'CR', 'MarkerEdgeColor', 'k'); hold on;
scatter(nan, nan, 'MarkerFaceColor',  [ 0.8457    0.7615    0.8738], 'DisplayName', 'TFIID', 'MarkerEdgeColor', 'k'); hold on;
% scatter(nan, nan, 'MarkerFaceColor',  [0.5961    0.3059    0.6392], 'DisplayName', 'others', 'MarkerEdgeColor', 'none')
legend('Position', [0.66 0.8205 0.1079 0.0824], 'Box', 'off', 'FontSize', 10)


clearvars cb wtMat nondbdMat promOrder promOrderName TFNo promAxes pughAxes hahnAxes p1 p2 prevEndIdx promClass hahnClass pughClass hahnClassNo pughClassNo cumPromType promTypes pughTypes
clearvars rebAxes zReb1 dbdMat opnDpnVals opnAxes
