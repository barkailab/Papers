descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};
 
 toCombine = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/toCombine.xlsx');
fullWTMat=[]; fullDBDMat=[];

for TFNo = 1:length(myCurrTable.TF)
    
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    
    
    %internal correlation WT
    repeatList = strrep(toCombine.sampleNames(~cellfun(@isempty, regexp(toCombine.sampleNames, [ '^' TFReq '_\d']))), '.mat', '');
    for i=1:length(repeatList)
        load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/checProfiles/' repeatList{i}], 'sumPromNew');
        sumPromNew(subtelomereGenes)=NaN;
        miniMat(:,i) =sumPromNew; clearvars sumPromNew
    end
    corrMat = corr(miniMat, 'rows', 'pairwise'); corrMat(find(diag(diag(corrMat)))) = [];
    lociCorrWTDBD.intWT(TFNo, 1) = mean(corrMat);
    fullWTMat(:, end+1:end+length(repeatList)) =miniMat;
    clearvars repeatList miniMat corrMat i
    
    %internal correlation DBD
    repeatList = strrep(toCombine.sampleNames(~cellfun(@isempty, regexp(toCombine.sampleNames, [ '^' DBDName '_\d']))), '.mat', '');
    for i=1:length(repeatList)
        load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/checProfiles/' repeatList{i}],  'sumPromNew');
        sumPromNew(subtelomereGenes)=NaN;
        miniMat(:,i) =sumPromNew; clearvars sumPromNew
    end
    corrMat = corr(miniMat, 'rows', 'pairwise'); corrMat(find(diag(diag(corrMat)))) = [];
    lociCorrWTDBD.intDBD(TFNo, 1) = mean(corrMat);
     fullDBDMat(:, end+1:end+length(repeatList)) =miniMat;
    clearvars repeatList miniMat corrMat i posMat motifPos   
    
end
clearvars CISBPnmer nmerRed n topN

%%

corrWTall= corr(fullWTMat, 'rows', 'pairwise');
a =1-corrWTall;
wtUpperMat =1-squareform(a); clearvars a

a =histogram(wtUpperMat, 'BinWidth', 0.025);
binCentersAllWT = mean([a.BinEdges(1:end-1); a.BinEdges(2:end)]);
binValsAllWT = a.Values/length(wtUpperMat);

b =histogram(lociCorrWTDBD.intWT, 'BinWidth', 0.025);
binCentersRepWT = mean([b.BinEdges(1:end-1); b.BinEdges(2:end)]);
binValsRepWT = b.Values/length(lociCorrWTDBD.intWT);

clearvars a b 


corrDBDall= corr(fullDBDMat, 'rows', 'pairwise');
a =1-corrDBDall;
dbdUpperMat =1-squareform(a); clearvars a

a =histogram(dbdUpperMat,'BinWidth', 0.025);
binCentersAllDBD = mean([a.BinEdges(1:end-1); a.BinEdges(2:end)]);
binValsAllDBD = a.Values/length(dbdUpperMat);

b =histogram(lociCorrWTDBD.intDBD, 'BinWidth', 0.025);
binCentersRepDBD = mean([b.BinEdges(1:end-1); b.BinEdges(2:end)]);
binValsRepDBD = b.Values/length(lociCorrWTDBD.intWT);
%%
figure('Renderer', 'painters', 'Color', [1 1 1], 'Position', [771 552 500 600]);
axes('Position', [0.200 0.5250 0.7 0.3])
bar(binCentersAllWT, binValsAllWT, 'DisplayName', 'random WT pairs', 'FaceColor', [0.9908    0.7793    0.5619], 'EdgeColor', 'none')
hold on
bar(binCentersRepWT, binValsRepWT, 'DisplayName', 'WT repeats', 'FaceColor', [0.1119    0.4913    0.2489 ], 'EdgeColor', 'none'); hold off
legend('Location', 'northwest'); %title('WT TFs')
ylabel('Fraction of samples'); %xlabel('correlation of signal over all promoters');
xticks([])
clearvars a b binCentersAllWT binCentersRepWT binValsAllWT binValsRepWT corrWTall wtUpperMat 


% figure('Renderer', 'painters', 'Color', [1 1 1], 'Position', [771 552 466 224]);
axes('Position', [0.200 0.200 0.7 0.3])
bar(binCentersAllDBD, binValsAllDBD, 'DisplayName', 'random DBD pairs', 'FaceColor', [0.9908    0.7793    0.5619], 'EdgeColor', 'none')%[0.9867  0.5463  0.3453]
hold on
bar(binCentersRepDBD, binValsRepDBD, 'DisplayName', 'DBD repeats', 'FaceColor', [0.1119    0.4913    0.2489 ] , 'EdgeColor', 'none')%[ 0.5716  0  0 ]
legend('Location', 'northwest'); %title('DBDs')
ylabel('Fraction of samples'); xlabel('correlation of signal over all promoters'); linkaxes;
clearvars a b binCentersAllDBD binCentersRepDBD binValsAllDBD binValsRepDBD corrDBDall dbdUpperMat 
