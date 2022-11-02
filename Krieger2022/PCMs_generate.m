%% generate PCMs

% common filter parameters
zeroDiag = 0;

% observed correlations
PCMs = struct;
for i = 1:length(F)
        disp(['gene-gene correlation matrix for: ', F{i}]);
        PCMs.(F{i}) = gene_gene_corr(dataS.(F{i}), geneThrs, expThrsList(i), ...
            corrThrs, NCorrThrs, zeroDiag);
end

for i = 1:length(F)
howManyGenes(i) = length(PCMs.(F{i}).genes);
end

F = fields(PCMs);
intG = 1:len_allGenes;
for i = 1:length(F)
    intG = intersect(intG, PCMs.(F{i}).genes);
end

F = fields(PCMs);
for i = 1:length(F)
    idx = find(ismember(PCMs.(F{i}).genes, intG));
    PCMs.(F{i}).data = PCMs.(F{i}).data(:, idx);
    PCMs.(F{i}).PCM = PCMs.(F{i}).PCM(idx, idx);
    PCMs.(F{i}).PCMflt = PCMs.(F{i}).PCMflt(idx, idx);
    PCMs.(F{i}).pval = PCMs.(F{i}).pval(idx, idx);
    PCMs.(F{i}).genes = PCMs.(F{i}).genes(idx);
end

Flabels = F;
for i = 1:length(Flabels)
    Flabels{i} = [F{i}, ' (', num2str(size(PCMs.(F{i}).data, 1)), ')'];
end

% histogram
figure; hold on;
for i = 1:length(F)
    d = reshape(PCMs.(F{i}).PCM, size(PCMs.(F{i}).PCM,1)^2, 1);
    x = -1:.1:1.1;
    y = histcounts(d, x);
    x = x(1:end-1);
    plot(x,y, 'linewidth', 4);
end
plot([0 0], ylim, 'k');
legend(Flabels);
box on; set(gcf, 'color', 'w');
xlabel('Correlation coefficient'); ylabel('Count');
title('Gene-gene correlation coefficients');