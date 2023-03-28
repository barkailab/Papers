descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable=descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
myCurrTable( ismember(myCurrTable.TF, {'CAT8';'PIP2';'ARO80';'CHA4';'MSN1'}),:)=[]; %dbdpoor
myCurrTable(cellfun(@isempty,myCurrTable.nonDBDName),:)=[];
myCurrTable(ismember(myCurrTable.TF, {'PHO4';'PHD1';'HAP4';'TYE7';'RTG3';'PUT3';'MGA1';'MOT3';'HOT1';'RTG1';'ZAP1';'ARO80'}),:)=[]; %nondbdpoor
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};

subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll);
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumProm.(fldNames{i}) (subtelomereGenes) = nan;
end
clearvars fldNames i

%% opn dpn tata

proxNuc = cell(6701,1); proxNuc(:)={''};
for i=1:12
    proxNuc(cell2mat(GP.groups{1,10}{1,2}(i))) = GP.groups{1,10}{1,1}(i);
end
promTata = cell(6701,1); promTata(:)={''};
for i=1:2
    promTata(cell2mat(GP.groups{1,11}{1,2}(i))) = GP.groups{1,11}{1,1}(i);
end
%%  promoter features- find top WT and DBD promoters and see their properties
 load('contOPNDPN')
contOPNDPN = -contOPNDPN;
expData = cat(2, GP.groups{4}{2}{:});
load('/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/newMedians2/promoterLengthsTS.mat')
load('/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/newMedians2/intGenicDist.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','dirStdClean')
for TFNo = 1:length(myCurrTable.TF)
    
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo}; nonDBDName = myCurrTable.nonDBDName{TFNo};
    
    wt(~isnan(medianSumProm.(TFReq))) = zscore(medianSumProm.(TFReq)(~isnan(medianSumProm.(TFReq))));
    dbd(~isnan(medianSumProm.(DBDName))) = zscore(medianSumProm.(DBDName)(~isnan(medianSumProm.(DBDName))));
    nondbd(~isnan(medianSumProm.(nonDBDName))) = zscore(medianSumProm.(nonDBDName)(~isnan(medianSumProm.(nonDBDName))));
    
    promNoWT= find(wt>2.5);
    promNoDBD= find(dbd>2.5);
     promNononDBD= find(nondbd>2.5);
    
    promProperties.TATA.name(TFNo,1) = {TFReq};
    promProperties.TATA.WT(TFNo,1) = sum(strcmp(promTata(promNoWT),{'TATA'}))/ length(promNoWT);
    promProperties.TATA.DBD(TFNo,1) = sum(strcmp(promTata(promNoDBD),{'TATA'}))/ length(promNoDBD);
     promProperties.TATA.nonDBD(TFNo,1) = sum(strcmp(promTata(promNononDBD),{'TATA'}))/ length(promNononDBD);
    
    promProperties.OPNDPN.name(TFNo,1) = {TFReq};
    promProperties.OPNDPN.WT(TFNo,1) = mean(contOPNDPN(promNoWT(isfinite(contOPNDPN(promNoWT)))))    ;
    promProperties.OPNDPN.DBD(TFNo,1) = mean(contOPNDPN(promNoDBD(isfinite(contOPNDPN(promNoDBD))))) ;
    promProperties.OPNDPN.nonDBD(TFNo,1) = mean(contOPNDPN(promNononDBD(isfinite(contOPNDPN(promNononDBD))))) ;
    

    
    promProperties.intGenDist.name(TFNo,1) = {TFReq};
    promProperties.intGenDist.WT(TFNo,1) =  median(intGeneDist(promNoWT), 'omitnan');
    promProperties.intGenDist.DBD(TFNo,1) =  median(intGeneDist(promNoDBD), 'omitnan');
    promProperties.intGenDist.nonDBD(TFNo,1) =  median(intGeneDist(promNononDBD), 'omitnan');
    
    promProperties.expGrp.WT(TFNo,1) = sum(sum(ismember(expData,promNoWT)).*(1:6))/ length(promNoWT);
    promProperties.expGrp.DBD(TFNo,1) = sum(sum(ismember(expData,promNoDBD)).*(1:6))/ length(promNoDBD);
    promProperties.expGrp.nonDBD(TFNo,1) = sum(sum(ismember(expData,promNononDBD)).*(1:6))/ length(promNononDBD);
    
    promProperties.expFlex.WT(TFNo,1) = median(abs(dirStdClean(promNoWT)), 'omitnan');
    promProperties.expFlex.DBD(TFNo,1) = median(abs(dirStdClean(promNoDBD)), 'omitnan');
    promProperties.expFlex.nonDBD(TFNo,1) = median(abs(dirStdClean(promNononDBD)), 'omitnan');
    
    clearvars maxLociDBD maxLociWT promNoDBD promNoWT TFReq DBDName wt dbd  promNononDBD  nonDBDName nondbd
end

%% plot OPN DPN tata data

useTable = promProperties.expFlex;

figure('Position', [2103 499 466 334], 'Color', [1 1 1], 'Renderer', 'painters', 'Visible', 'on')
scatter(useTable.WT, useTable.DBD, 'filled', 'DisplayName', 'DBD', 'MarkerFaceColor', [ 0.5648    0.8186    0.5141])
hold on; scatter(useTable.WT, useTable.nonDBD, 'filled', 'DisplayName', 'nonDBD', 'MarkerFaceColor', [  0.0372    0.3911    0.6582])
xlabel('WT'); ylabel('(non)DBD'); title('intGenic Dist')
xlim([min([useTable.WT; useTable.DBD; useTable.nonDBD]) max([useTable.WT; useTable.DBD; useTable.nonDBD])])
ylim([min([useTable.WT; useTable.DBD; useTable.nonDBD]) max([useTable.WT; useTable.DBD; useTable.nonDBD])])
% legend


% xlim([200 1600]); ylim([200 1600]);

% xlim([0.1 0.6]); ylim([0.1 0.6]);

% xlim([-1.5 1]); ylim([-1.5 1]);

 xlim([0.1 0.24]); ylim([0.1 0.24]);

xticks(xlim); yticks(ylim)
set(gcf, 'Position', [1340 511 341 243])
hold on; plot(xlim, ylim, '--k', 'HandleVisibility', 'off')
clearvars useTable
