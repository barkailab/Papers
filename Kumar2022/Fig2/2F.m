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
    clearvars  TFReq nonDBDName
end
WTDBD = struct2table(WTDBD);
WTDBD.colorCode(ismember(WTDBD.family, 'Zinc cluster')) = {[1 0 0 ]};
WTDBD.colorCode(ismember(WTDBD.family, 'C2H2 ZF')) = {[0 1 0 ]};
WTDBD.colorCode(ismember(WTDBD.family, 'bHLH')) = {[0 0 1 ]};
WTDBD.colorCode(~ismember(WTDBD.family, {'Zinc cluster';'C2H2 ZF';'bHLH'})) = {[ 0.5450    0.5149    0.7399]};

WTDBD = sortrows(WTDBD,'corr','descend');
figure(  'Renderer', 'painters', 'Color', [1 1 1]);
% axes('Position',  [0.500 0.10 0.4350 0.87])
myBar =bar(WTDBD.corr, 'FaceColor', 'flat', 'BarWidth', 1);
% set(gca, 'YDir', 'reverse')
xlim([0.5 length(myCurrTable.TF)+0.5])
xticks(1:length(myCurrTable.TF)); xticklabels(WTDBD.name)
fillColor = (brewermap(500, 'GnBu'));
fillIdx = round(rescale(WTDBD.intCorrDBD, 1, 500, 'InputMin', 0.7, 'InputMax', 1));
myBar.CData = fillColor(fillIdx, :);

% myBar.CData = cell2mat(WTDBD.colorCode);
ylabel('DBD-WT corr'); 

cb = colorbar('Location', 'eastoutside');colormap(fillColor); xlabel(cb, 'nonDBD repeat corr'); caxis([0.7 1]);cb.Ticks=[0.7 1];
%%
toCombine = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/toCombine.xlsx');
fullWTMat=[]; fullDBDMat=[];
for TFNo = 1:length(myCurrTable.TF)
    
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    
    
    %internal correlation WT
    repeatList = strrep(toCombine.sampleNames(~cellfun(@isempty, regexp(toCombine.sampleNames, [ '^' TFReq '_\d']))), '.mat', '');
    for i=1:length(repeatList)
        load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/checProfiles/' repeatList{i}], 'sumPromNew')
        miniMat(:,i) =sumPromNew; clearvars sumPromNew
    end
    corrMat = corr(miniMat, 'rows', 'pairwise'); corrMat(find(diag(diag(corrMat)))) = [];
    intCorr.WT(TFNo, 1) = mean(corrMat);
    fullWTMat(:, end+1:end+length(repeatList)) =miniMat;
    clearvars repeatList miniMat corrMat i
    
    %internal correlation DBD
    repeatList = strrep(toCombine.sampleNames(~cellfun(@isempty, regexp(toCombine.sampleNames, [ '^' DBDName '_\d']))), '.mat', '');
    for i=1:length(repeatList)
        load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/checProfiles/' repeatList{i}],  'sumPromNew')
        miniMat(:,i) = sumPromNew; clearvars sumPromNew
    end
    corrMat = corr(miniMat, 'rows', 'pairwise'); corrMat(find(diag(diag(corrMat)))) = [];
   intCorr.DBD(TFNo, 1) = mean(corrMat);
     fullDBDMat(:, end+1:end+length(repeatList)) =miniMat;
    clearvars repeatList miniMat corrMat i posMat motifPos   
    
end

%%
corrWTall= corr(fullWTMat, 'rows', 'pairwise');
a =1-corrWTall;
wtUpperMat =1-squareform(a); clearvars a
corrDBDall= corr(fullDBDMat, 'rows', 'pairwise');
a =1-corrDBDall;
dbdUpperMat =1-squareform(a); clearvars a
a =histogram([wtUpperMat dbdUpperMat], 'BinWidth', 0.025);
binCentersAllWT = mean([a.BinEdges(1:end-1); a.BinEdges(2:end)]);
binValsAllWT = a.Values/length([wtUpperMat dbdUpperMat]);
b =histogram([intCorr.WT intCorr.DBD], 'BinWidth', 0.025);
binCentersRepWT = mean([b.BinEdges(1:end-1); b.BinEdges(2:end)]);
binValsRepWT = b.Values/(length(intCorr.WT)*2); clearvars a b
a =histogram( WTDBD.corr ,'BinWidth', 0.05);
binCenters = mean([a.BinEdges(1:end-1); a.BinEdges(2:end)]);
binVals = a.Values/length(WTDBD.corr); clearvars a

figure(  'Renderer', 'painters', 'Color', [1 1 1]);
plot(binCenters, cumsum(binVals, 'reverse'), 'DisplayName', 'WT-DBD', 'Color', [0.4923    0.0320    0.4684])
xlabel('corr.'); ylabel('fraction of TFs');
hold on; plot(binCentersAllWT, cumsum(binValsAllWT, 'reverse'),'Color', [0.9867  0.5463  0.3453], 'DisplayName', 'random')
hold on; plot(binCentersRepWT, cumsum(binValsRepWT, 'reverse'),'Color', [ 0.5716  0  0 ], 'DisplayName', 'repeats')
legend('Location', 'bestoutside')
clearvars a b binCentersAllWT binCentersRepWT binValsAllWT binValsRepWT corrWTall wtUpperMat 
clearvars a b binCentersAllDBD binCentersRepDBD binValsAllDBD binValsRepDBD corrDBDall dbdUpperMat binCenters binVals

