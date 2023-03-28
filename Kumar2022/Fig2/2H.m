subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll);
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumProm.(fldNames{i}) (subtelomereGenes) = nan;
end
clearvars fldNames i

descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
myCurrTable(ismember(myCurrTable.FAMILY, 'Zinc cluster'),:)= [];
 myCurrTable(ismember(myCurrTable.TF, {'PIP2';'ARO80';'CAT8';'CHA4';'PHD1';'MSN1';'GCR1'}),:)= [];
 myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};
 

load('/home/labs/barkailab/divyakr/matlab/Disorder tendency/disTendency.mat');

for TFNo = 1:length(myCurrTable.TF)
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    tbl.name(TFNo,1) = {TFReq};
    tbl.family(TFNo,1) =  myCurrTable.FAMILY(TFNo);
    x =corr([medianSumProm.(TFReq) medianSumProm.(DBDName)], 'rows', 'pairwise');
    tbl.corr(TFNo,1) = x(2);
     tbl.corr2(TFNo,1) = myCurrTable.DBDvsWT(TFNo); clearvars x TFReq DBDName;
    
end

myCurrTable.DBDname(ismember(myCurrTable.TF, 'SUM1'))= {'SUM1_d524_1062'};
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d141_650'};
for TFNo = 1:length(myCurrTable.TF)
    TFReq = myCurrTable.TF{TFNo}; 
    t = find(ismember(GP.gene_infoR64.name, TFReq));
    currDis = disTendency(:,t);
    
    dNos=   str2double(regexp(myCurrTable.DBDname{TFNo},'(?<=idr)\d+|(?<=d)\d+|(?<=_)\d{2,4}','match'));
    if contains(myCurrTable.DBDname{TFNo}, '_d')
    currDis(1:dNos(1)) =0; currDis(dNos(2):myCurrTable.ProteinLength(TFNo)) =0;
    end
    if contains(myCurrTable.DBDname{TFNo}, '_idr')
    currDis(dNos(1):dNos(2)) =0;
    end
    tbl.idrLength(TFNo,1) = sum(currDis>0.5);
    tbl.idrAmt(TFNo,1) = sum(currDis, 'omitnan');
    tbl.nonDBDLength(TFNo,1) =myCurrTable.ProteinLength(TFNo) -myCurrTable.DBDTaken(TFNo);
    clearvars t TFReq currDis
end
myCurrTable.DBDname(ismember(myCurrTable.TF, 'SUM1'))= {'SUM1_CT_del'};
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};

lowLightIdx =[4 39];
scatter( tbl.idrLength, 1-tbl.corr2,60, tbl.nonDBDLength-tbl.idrLength, 'filled', 'DisplayName', 'Others', 'MarkerEdgeColor', [0 0 0])
% hold on; scatter( tbl.idrLength, 1-tbl.corr2,60,  'MarkerEdgeColor', [1 1 1], 'HandleVisibility', 'off')
hold on; scatter( tbl.idrLength(lowLightIdx) , 1-tbl.corr2(lowLightIdx), 60, tbl.nonDBDLength(lowLightIdx)-tbl.idrLength(lowLightIdx),'filled','MarkerFaceColor', [0.5 0.5 0.5],  'MarkerFaceAlpha', 1)
ylabel('nonDBD effect'); xlabel('IDR length deleted')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
colorbar; ylabel(colorbar, 'str nonDBD length')
fillColor =flipud(brewermap(100, 'BrBG'));
colormap(gca, fillColor(1:100,:)); caxis([0 500]);  ylim([0 1])
xi =tbl.idrLength; xi(lowLightIdx) = [];
yi =1-tbl.corr2; yi(lowLightIdx) = [];
p =linortfit2(xi, yi);
x=linspace(47, 573);
hold on; plot(x, (p(1)*x) +p(2), ':k'); clearvars x p xi yi
%%
scatter( tbl.idrLength, 1-tbl.corr2,'filled', 'MarkerFaceColor', [0.4 0.4 0.4], 'DisplayName', 'Others')
ZnClIdx = ismember(tbl.family, 'Zinc cluster');
hold on; scatter( tbl.idrLength(ZnClIdx), 1-tbl.corr2(ZnClIdx),'filled', 'MarkerFaceColor', [0 0 1], 'DisplayName', 'Zinc cluster'); 
bHLHIdx = ismember(tbl.family, 'bHLH');
hold on; scatter( tbl.idrLength(bHLHIdx), 1-tbl.corr2(bHLHIdx),'filled', 'MarkerFaceColor', [0 1 0], 'DisplayName', 'bHLH'); 
ZnFnIdx = ismember(tbl.family, 'C2H2 ZF');
hold on; scatter( tbl.idrLength(ZnFnIdx), 1-tbl.corr2(ZnFnIdx),'filled', 'MarkerFaceColor', [1 0 0], 'DisplayName', 'Zn C2H2'); 
clearvars ZnClIdx bHLHIdx ZnFnIdx
ylabel('1- WT:DBD corr'); xlabel('idr length deleted')
legend('Location', 'bestoutside')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1])

figure;
scatter(tbl.idrAmt, 1-tbl.corr2, 'filled', 'MarkerFaceColor', [0.4 0.4 0.4], 'DisplayName', 'Others')
ZnClIdx = ismember(tbl.family, 'Zinc cluster');
hold on; scatter(tbl.idrAmt(ZnClIdx), 1-tbl.corr2(ZnClIdx), 'filled', 'MarkerFaceColor', [0 0 1], 'DisplayName', 'Zinc cluster'); 
bHLHIdx = ismember(tbl.family, 'bHLH');
hold on; scatter( tbl.idrAmt(bHLHIdx), 1-tbl.corr2(bHLHIdx),'filled', 'MarkerFaceColor', [0 1 0], 'DisplayName', 'bHLH'); 
ZnFnIdx = ismember(tbl.family, 'C2H2 ZF');
hold on; scatter( tbl.idrAmt(ZnFnIdx), 1-tbl.corr2(ZnFnIdx),'filled', 'MarkerFaceColor', [1 0 0], 'DisplayName', 'Zn C2H2'); 
clearvars ZnClIdx bHLHIdx ZnFnIdx
ylabel('1- WT:DBD corr'); xlabel('idr amount deleted')
legend('Location', 'bestoutside')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1])

figure;
scatter(tbl.nonDBDLength, 1-tbl.corr2, 'filled', 'MarkerFaceColor', [0.4 0.4 0.4], 'DisplayName', 'Others')
ZnClIdx = ismember(tbl.family, 'Zinc cluster');
hold on; scatter(tbl.nonDBDLength(ZnClIdx), 1-tbl.corr2(ZnClIdx), 'filled', 'MarkerFaceColor', [0 0 1], 'DisplayName', 'Zinc cluster'); 
bHLHIdx = ismember(tbl.family, 'bHLH');
hold on; scatter( tbl.nonDBDLength(bHLHIdx), 1-tbl.corr2(bHLHIdx),'filled', 'MarkerFaceColor', [0 1 0], 'DisplayName', 'bHLH'); 
ZnFnIdx = ismember(tbl.family, 'C2H2 ZF');
hold on; scatter( tbl.nonDBDLength(ZnFnIdx), 1-tbl.corr2(ZnFnIdx),'filled', 'MarkerFaceColor', [1 0 0], 'DisplayName', 'Zn C2H2'); 
clearvars ZnClIdx bHLHIdx ZnFnIdx
ylabel('1- WT:DBD corr'); xlabel('total deleted')
legend('Location', 'bestoutside')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1])