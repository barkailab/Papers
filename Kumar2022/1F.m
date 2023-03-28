figure('Position', [685 42 300 953])
myCurrTable = sortrows(myCurrTable,'DBDintCorr','descend');
myCurrTable(myCurrTable.DBDintCorr>0.77, :) = sortrows(myCurrTable(myCurrTable.DBDintCorr>0.77, :),'ProteinLength','descend');
myCurrTable(myCurrTable.DBDintCorr<=0.77, :) = sortrows(myCurrTable(myCurrTable.DBDintCorr<=0.77, :),'ProteinLength','descend');

skipList1 =find(contains(myCurrTable.TF, myCurrTable.TF(myCurrTable.DBDintCorr<=0.77, :))); 
skipList2 =find(contains(myCurrTable.TF, myCurrTable.TF(ismember(myCurrTable.FAMILY, 'Zinc cluster'))));
skipList3 =find(contains(myCurrTable.TF, {'ADR1'; 'SWI4';'GLN3'}));
skipList4 =find(contains(myCurrTable.TF, {'AZF1';'CBF1';'FKH2';'RFX1';'SOK2';'SKO1'}));
myCurrTable.DBDname(ismember(myCurrTable.TF, 'SUM1'))= {'SUM1_d524_1062'};


axes('Position', [0.2500 0.1100 0.7 0.8150])
myBar = barh(myCurrTable.ProteinLength, 'BarWidth', 1, 'FaceColor', 'flat', 'EdgeColor', [0.7498    0.5067    0.1773]);
ylim([0.5 length(myCurrTable.TF)+.5]);
xlim([0 max(myCurrTable.ProteinLength)+100]);
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d141_650'};
for TFNo=1:length(myCurrTable.TF)
    hold on
    a = [TFNo-0.5 TFNo+0.5];
    dNos=   str2double(regexp(myCurrTable.DBDname{TFNo},'(?<=idr)\d+|(?<=d)\d+|(?<=_)\d{2,4}','match'));
    if ~contains(myCurrTable.DBDname{TFNo}, 'idr')
        if dNos(2)<myCurrTable.ProteinLength(TFNo)
            dNos(1) =dNos(2)+1;
            dNos(2) =myCurrTable.ProteinLength(TFNo);
        else
            dNos(2) =dNos(1)-1;
            dNos(1) =1;
        end
    end
    fill([dNos(1), dNos(1), dNos(2), dNos(2)], [a a(2:-1:1) ],[0.1765  0 0.2941], 'EdgeColor', 'none' )%, 'FaceAlpha', 0.5)
    fill([myCurrTable.FinalStart(TFNo), myCurrTable.FinalStart(TFNo),  myCurrTable.FinalEnd(TFNo), myCurrTable.FinalEnd(TFNo)], [a a(2:-1:1) ],[0.5020 0.4510 0.6745] , 'FaceAlpha', 1, 'EdgeColor', 'none')
clearvars a  dNos
end
    

set(gca, 'YDir', 'reverse'); set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
yticks([1:length(myCurrTable.TF)]); yticklabels(myCurrTable.TF);
xlabel('TF size (amino acids)')
unusedBar= nan(size(myCurrTable.ProteinLength));
  unusedBar(skipList1)=myCurrTable.ProteinLength(skipList1);
  hold on; barh(unusedBar, 'BarWidth', 1, 'FaceColor', [1 1 1], 'EdgeColor', [0.7498    0.5067    0.1773], 'FaceAlpha', 0.7);


axes('Position', [0.0400 0.1100 0.1 0.8150])
scatter(linspace(1,1, length(skipList3)),length(myCurrTable.TF)-skipList3+1, [], 'filled', 'MarkerFaceColor', [0.1019    0.5187    0.2472],  'DisplayName', 'small'); 
ylim([0.5 length(myCurrTable.TF)+.5]);axis off
axes('Position', [0.0400 0.1100 0.1 0.8150])
hold on; scatter(linspace(1,1, length(skipList2)),length(myCurrTable.TF)-skipList2+1, [], 'filled', 'MarkerFaceColor',[0.7498    0.5067    0.1773],'DisplayName', 'ZnC6 TFs'); 
ylim([0.5 length(myCurrTable.TF)+.5]);axis off
axes('Position', [0.0400 0.1100 0.1 0.8150])
hold on; scatter(linspace(1,1, length(skipList4)),length(myCurrTable.TF)-skipList4+1, [], 'filled', 'MarkerFaceColor',[0.5020 0.4510 0.6745],'DisplayName', 'double sided del'); 
ylim([0.5 length(myCurrTable.TF)+.5]);axis off



fillColor = flipud(brewermap(500, 'BrBG'));
fillColor = fillColor(250:350,:);
myBar.CData=[0.9505    0.8788    0.6936];

myCurrTable.DBDname(ismember(myCurrTable.TF, 'SUM1'))= {'SUM1_CT_del'};
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};
