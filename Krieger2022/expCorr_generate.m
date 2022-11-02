%% expected correlation

Nshuff = 10;
sizes = cell2mat(struct2cell(structfun(@size, dataS, 'uniformoutput', false)));
Nexpected = floor(min(sizes(1:4, 2))./2);
expCorr = struct;
for i = 1:length(F)
    data = dataS.(F{i});
    disp(['gene-gene correlation matrix for: ', F{i}]);
    expCorr.(F{i}) = expected_corr_general_fixedN(data, Nshuff, ...
        geneThrs, expThrsList(i), corrThrs, NCorrThrs, zeroDiag, ifFlt, Nexpected);
end

howManyGeneExp = nan(length(F), 1);
for i = 1:length(F)
    howManyGeneExp(i) = length(expCorr.(F{i}).genes);
end

intGExp = 1:len_allGenes;
for i = 1:length(F)
    intGExp = intersect(intGExp, expCorr.(F{i}).genes);
end

for i = 1:length(F)
    idx = find(ismember(expCorr.(F{i}).genes, intGExp));
    expCorr.(F{i}).genes = expCorr.(F{i}).genes(idx);
    expCorr.(F{i}).corrDiag = expCorr.(F{i}).corrDiag(idx, :);
    expCorr.(F{i}).PCMhalves{1} = expCorr.(F{i}).PCMhalves{1}(idx, idx);
    expCorr.(F{i}).PCMhalves{2} = expCorr.(F{i}).PCMhalves{2}(idx, idx);
end

intG = PCMs.(F{1}).genes;
intint = intersect(intG, intGExp);

for i = 1:length(F)
    idx = find(ismember(expCorr.(F{i}).genes, intint));
    expCorr.(F{i}).genes = expCorr.(F{i}).genes(idx);
    expCorr.(F{i}).corrDiag = expCorr.(F{i}).corrDiag(idx, :);
    expCorr.(F{i}).PCMhalves{1} = expCorr.(F{i}).PCMhalves{1}(idx, idx);
    expCorr.(F{i}).PCMhalves{2} = expCorr.(F{i}).PCMhalves{2}(idx, idx);
    
    idx = find(ismember(PCMs.(F{i}).genes, intint));
    PCMs.(F{i}).data = PCMs.(F{i}).data(:, idx);
    PCMs.(F{i}).PCM = PCMs.(F{i}).PCM(idx, idx);
    PCMs.(F{i}).PCMflt = PCMs.(F{i}).PCMflt(idx, idx);
    PCMs.(F{i}).pval = PCMs.(F{i}).pval(idx, idx);
    PCMs.(F{i}).genes = PCMs.(F{i}).genes(idx);
end
