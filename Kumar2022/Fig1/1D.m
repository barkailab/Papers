descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];

finalAll =readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/FinalLengthsAll.xlsx', 'Sheet', 'Sheet2');
familyCat ={'Zinc cluster'};
% colorList = {[0.7725    0.1059    0.4902]; [0.9935    0.6214    0.2613]; [0.3967    0.2803    0.6210]};
colorList = {[0.9935    0.6214    0.2613]};

for i=1:1
    currList = ismember(finalAll.TFFAMILY, familyCat{i});
    hold on; scatter(finalAll.nonDBDStr(currList), finalAll.nonDBDIdr(currList), 50, colorList{i}, 'filled', 'DisplayName', familyCat{i});
    clearvars currList
end
otherList = ~ismember(finalAll.TFFAMILY, familyCat);
hold on; scatter(finalAll.nonDBDStr(otherList), finalAll.nonDBDIdr(otherList),50, 'filled', 'DisplayName', 'Others', 'MarkerFaceColor', [0.6 0.6 0.6] );
clearvars otherList
legend; set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
xlabel('nonDBD ordered'); ylabel('nonDBD disordered');

selectList =ismember(finalAll.SGDGENENAME, myCurrTable.TF);
hold on; scatter(finalAll.nonDBDStr(selectList), finalAll.nonDBDIdr(selectList),50, 'o','k', 'HandleVisibility', 'off', 'LineWidth', 1.5);
clearvars selectList

ylim(xlim); hold on; plot(xlim, ylim, '--k', 'HandleVisibility', 'off')