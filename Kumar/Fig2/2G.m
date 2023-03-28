subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll);
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumProm.(fldNames{i}) (subtelomereGenes) = nan;
end
clearvars fldNames i


% famSelect = 'bHLH';
% TFList ={'ADR1';'MSN2';'CRZ1';'YAP1';'RLM1';'GIS1'};
% TFList ={'SFL1';'MIG1';'IXR1';'DOT6';'TOD6';'SOK2';'SKO1'};
% TFList ={'MIG1';'IXR1';'DOT6';'TOD6';'SOK2';'SKO1'};
TFList ={'ASH1';'SWI5';'ACE2';'FKH1';'FKH2';'MBP1';'SWI4'};

% famTable = myCurrTable(ismember(myCurrTable.FAMILY, famSelect), :);
famTable = myCurrTable(ismember(myCurrTable.TF, TFList), :);

famTable = sortrows(famTable,'DBDvsWT','descend');
k1=1;
for TFNo = 1:length(famTable.TF)
    currTF = famTable.TF{TFNo}; currDBD = famTable.DBDname{TFNo};
    zWT(~isnan(medianSumProm.(currTF))) = zscore(medianSumProm.(currTF)(~isnan(medianSumProm.(currTF))));
    zDBD(~isnan(medianSumProm.(currDBD))) = zscore(medianSumProm.(currDBD)(~isnan(medianSumProm.(currDBD))));
    
    [~,b]=maxk(max([zWT; zDBD]), 60);
[~,d]=sort(zWT(b)./zDBD(b), 'descend');
promOrder(k1:k1+59)= b(d);
k1=60+k1;
%     p = unique([find(zWT>2.5) find(zDBD>2.5)]);
%     k2 = length(unique([find(zWT>2.5) find(zDBD>2.5)]));
%     [~, temp ]= sort(zWT(p), 'descend');
%     promOrder(k1:k2+k1-1)= p(temp);
%     k1=k2+k1;
    clearvars k2 p zWT zDBD currTF currDBD temp b d
end

for TFNo = 1:length(famTable.TF)
currTF = famTable.TF{TFNo}; currDBD = famTable.DBDname{TFNo};
    zWT(~isnan(medianSumProm.(currTF))) = zscore(medianSumProm.(currTF)(~isnan(medianSumProm.(currTF))));
    zDBD(~isnan(medianSumProm.(currDBD))) = zscore(medianSumProm.(currDBD)(~isnan(medianSumProm.(currDBD))));
    imageMat(:, (TFNo*2)-1)= zWT(promOrder);
    imageMat(:, (TFNo*2))= zDBD(promOrder);
    clearvars zWT zDBD  currTF currDBD
end

imagesc(imageMat')
yticks(1.5:2:(length(famTable.TF)*2 -0.5))
yticklabels(famTable.TF)
% yticks(2:2:(sum(ismember(myCurrTable.FAMILY, famSelect))*2))
% yticklabels(famTable.DBDname)
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1])
colorbar; colormap(brewermap(100, 'Reds')); caxis([0 10])


for TFNo = 1:length(famTable.TF)
hold on; plot([60*TFNo 60*TFNo], ylim, 'k');
hold on; plot([xlim], [0.5+TFNo*2 0.5+TFNo*2], 'k');
end
clearvars famTable famSelect k1 imageMat promOrder TFNo

%%
% TFList = {'ADR1';'MSN2';'CRZ1';'YAP1';'RLM1';'GIS1'};
% TFList ={'MIG1';'IXR1';'DOT6';'TOD6';'SOK2';'SKO1'};
TFList ={'ASH1';'SWI5';'ACE2';'FKH1';'FKH2';'MBP1';'SWI4'};

famTable = myCurrTable(ismember(myCurrTable.TF, TFList), :);
famTable = sortrows(famTable,'DBDvsWT','descend');
TFList = famTable.TF; clearvars famTable
sumPromMini=nan(6701, length(TFList));

for i=1:length(TFList)
    sumPromMini(:,(i*2)-1) = medianSumProm.(TFList{i});
    DBDName = myCurrTable.DBDname(ismember(myCurrTable.TF, TFList{i}));
    sumPromMini(:,(i*2)) = medianSumProm.(char(DBDName)); 
    ylabs((i*2)-1,1) = TFList(i); ylabs(i*2, 1) = DBDName;clearvars DBDName
end
figure('Color', [1 1 1]);
imagesc(corr(sumPromMini, 'rows','pairwise'));
colormap((brewermap(200, 'OrRd')))
caxis([0 1]); xticks([])
yticks([1:length(ylabs)]); yticklabels(strrep(ylabs, '_' , ' '));
set(gca,'FontSize',11)
clearvars TFList sumPromMini i ylabs
