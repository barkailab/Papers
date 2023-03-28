descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
myCurrTable(ismember(myCurrTable.FAMILY, 'Zinc cluster'),:)= [];
 myCurrTable(ismember(myCurrTable.TF, {'PIP2';'ARO80';'CAT8';'CHA4';'PHD1';'MSN1';'GCR1';'HOT1'}),:)= [];
 %no in vitro motif available for HOT1
 myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};

%% new version tbl maker

% needs whichProm, promVecIdx, GP, CISBP, medianMotif
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
promoterIDXvec(ismember(whichProm,subtelomereGenes))=0;
for TFNo=1:length(myCurrTable.TF)
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    
     [ ~, topIVMIdxall] = maxk(CISBP7merDist.(TFReq), 10);
    [topIVMValsWT, topIVMIdxWT] = max(medianMotifNew.(TFReq)(topIVMIdxall));  topIVMIdxWT = topIVMIdxall( topIVMIdxWT );
    [topIVMValsDBD, topIVMIdxDBD] = max(medianMotifNew.(DBDName)(topIVMIdxall));  topIVMIdxDBD = topIVMIdxall( topIVMIdxDBD );
    x= zscore(medianMotifNew.(TFReq));
    tbl.selBiasWT(TFNo,1) =x(topIVMIdxWT); clearvars x;
    tbl.name(TFNo,1) ={TFReq};
    x= zscore(medianMotifNew.(DBDName));
    tbl.selBiasDBD(TFNo,1) = x(topIVMIdxDBD);clearvars x;
    inVitroZscore = CISBP7merDist.(TFReq);
    tbl.maxCISBPScoreWT(TFNo,1) =inVitroZscore(topIVMIdxWT);
        tbl.maxCISBPScoreDBD(TFNo,1) =inVitroZscore(topIVMIdxDBD);
    clearvars n topIVMIdx topIVMVals TFReq DBDName topIVMIdxall topIVMValsWT topIVMValsDBD topIVMIdxWT  topIVMIdxDBD inVitroZscore
    
end
%%
scatter((tbl.selBiasWT),(tbl.selBiasDBD), 200,mean([tbl.maxCISBPScoreWT  tbl.maxCISBPScoreDBD], 2), 'filled')
% addMrkr(tbl.selBiasWT',tbl.selBiasDBD', 1.25 , tbl.maxCISBPScoreDBD,  tbl.maxCISBPScoreWT)
hold on; text(tbl.selBiasWT(highlightIdx)+0.7 , tbl.selBiasDBD(highlightIdx)+0.7, tbl.name(highlightIdx),  'FontWeight' , 'bold');
scatter(tbl.selBiasWT(highlightIdx) , tbl.selBiasDBD(highlightIdx), 200, 'MarkerEdgeColor', 'k')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1], 'Position', [680 558 560 420])%, 'Position'  , [2185 292 828 672]);
xlabel('WT'); ylabel('DBD');
colormap(brewermap(100, 'Greens'))
cb =colorbar;  ylabel(cb, '{\it in vitro} score'); %caxis([3 6])
xlim([-1 45]); ylim(xlim); xticks([0  10  20  30 40]); yticks(xticks)
caxis([3 15]); cb.Ticks =[3 15]; hold on
plot(xlim, ylim, '--k')
