TFReq = 'GLN3'; 
DBDName = char(myCurrTable.DBDname(ismember(myCurrTable.TF, TFReq)));
zWT = zscore(medianMotifNew.(TFReq)); zDBD = zscore(medianMotifNew.(DBDName));
scatter(zWT, zDBD, 100, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7])
xlabel('7mer zscore WT'); ylabel('7mer zscore DBD');set(gcf, 'Renderer', 'painters', 'Color', [1 1 1], 'Position', [785 650 405 292]);
[ ~, topTenIdx] = maxk(CISBP7merDist.(TFReq), 10); title(TFReq)
[~, x]=max(medianMotifNew.(TFReq)(topTenIdx)); h(1) = topTenIdx(x); clearvars x
[~, x]=max(medianMotifNew.(DBDName)(topTenIdx)); h(2) = topTenIdx(x); clearvars x
hold on; scatter(zWT(h), zDBD(h),100, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
% hold on; text(zWT(h), zDBD(h),nmerRed7.seq(h));
clearvars topTenIdx TFReq DBDName
xlim([-3 40]); ylim(xlim)