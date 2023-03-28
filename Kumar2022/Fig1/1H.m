%% main text repeats example

TFNo = 17;
    
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    
 repeatList = strrep(toCombine.sampleNames(~cellfun(@isempty, regexp(toCombine.sampleNames, [ '^' DBDName '_\d']))), '.mat', '');
    for i=1:2
        load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/checProfiles/' repeatList{i}], 'sumPromNew');
         r(i, 1:6701) = nan;
        sumPromNew(subtelomereGenes)=NaN; r(i,(~isnan(sumPromNew))) = zscore(sumPromNew(~isnan(sumPromNew))); clearvars sumPromNew
    end
% scatter(r(1,:), r(2,:),200,  'filled', 'MarkerFaceColor', [0.1111    0.5191    0.2490], 'MarkerEdgeColor', [0 0 0]) %   0.7400    0.8938    0.5854]
% xlabel([TFReq ' rep.1']); ylabel([TFReq ' rep.2']);

scatter(r(1,:), r(2,:),200,  'filled', 'MarkerFaceColor', [0.7400    0.8938    0.5854], 'MarkerEdgeColor', [0 0 0]) %   0.7400    0.8938    0.5854]
xlabel([TFReq ' DBD rep.1']); ylabel([TFReq ' DBD rep.2']);
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
gname(GP.gene_infoR64.name)