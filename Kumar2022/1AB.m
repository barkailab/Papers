descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];


finalAll =readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/FinalLengthsAll.xlsx', 'Sheet', 'Sheet2');
ZnClusAll = finalAll(ismember(finalAll.TFFAMILY, 'Zinc cluster'),:);
ZnClusAll = sortrows(ZnClusAll,'TOTALLENGTH','descend');
figure; barh(ZnClusAll.TOTALLENGTH, 'FaceColor', [0.8 0.8 0.8],'EdgeColor', [0.7 0.7 0.7], 'BarWidth', 1, 'DisplayName', 'struc. nonDBD')
% figure; barh(ZnClusAll.TOTALLENGTH, 'FaceColor', [0.8 0.8 0.8],'EdgeColor', [0.65    0.65    0.65], 'BarWidth', 1, 'DisplayName', 'struc. nonDBD')
hold on; barh(ZnClusAll.DBDLENGTH+ZnClusAll.nonDBDIdr, 'FaceColor', [0.0069    0.3040    0.5961], 'EdgeColor', [0.7 0.7 0.7], 'BarWidth', 1, 'DisplayName', 'IDR nonDBD')
hold on; barh(ZnClusAll.DBDLENGTH, 'FaceColor', [0.3712    0.7333    0.4273], 'EdgeColor', [0.7 0.7 0.7], 'BarWidth', 1, 'DisplayName', 'DBD')
hold on;scatter(ZnClusAll.TOTALLENGTH(ismember(ZnClusAll.SGDGENENAME, myCurrTable.TF))+20, find(ismember(ZnClusAll.SGDGENENAME, myCurrTable.TF)), '_', 'HandleVisibility', 'off');
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
 ylim([0.5 45.5]);xlim([0 1500]); xticks(xlim); xlabel('Amino Acids'); legend

finalAll(ismember(finalAll.TFFAMILY, 'Zinc cluster'),:)=[];
finalAll = sortrows(finalAll,'TOTALLENGTH','descend');
figure; barh(finalAll.TOTALLENGTH,'FaceColor', [0.8 0.8 0.8],'EdgeColor', [0.7 0.7 0.7], 'BarWidth', 1, 'DisplayName', 'struc. nonDBD')
hold on; barh(finalAll.DBDLENGTH+finalAll.nonDBDIdr, 'FaceColor', [0.0069    0.3040    0.5961], 'EdgeColor', [0.7 0.7 0.7], 'BarWidth', 1, 'DisplayName', 'IDR nonDBD')
hold on; barh(finalAll.DBDLENGTH, 'FaceColor', [0.3712    0.7333    0.4273], 'EdgeColor',[0.7 0.7 0.7], 'BarWidth', 1, 'DisplayName', 'DBD')
hold on;scatter(finalAll.TOTALLENGTH(ismember(finalAll.SGDGENENAME, myCurrTable.TF))+20, find(ismember(finalAll.SGDGENENAME, myCurrTable.TF)), '_', 'HandleVisibility', 'off');
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
ylim([0.5 102.5]);xlim([0 1500]); xticks(xlim); xlabel('Amino Acids'); legend