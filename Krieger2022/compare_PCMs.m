% compare pairwise correlation matrices (PCMs)
%% load data
exp_per_sp_table = 'expression_per_species_test.xlsx';
load_expression_data;
F = fields(dataS);
%% parameters
global PCMs expCorr F intint ifFlt corrThrs NCorrThrs len_allGenes
ifFlt = 0; % whether to filter the pairwise correlation matrix
corrThrs = 0.2; % if filtering, take only |R| > 'corrThrs'
geneThrs = 10; % filter out genes that are experssed in only 'geneThrs'% of the samples (or less)
NCorrThrs = 20; % min number of genes to compute correlation
expThrs = min(min(dataS.(F{1}))); % minimal expression level
expThrsList = repmat(expThrs, 1, length(F));
len_allGenes = size(dataS.(F{1}), 1);
%% compare PCM of two datasets
disp('generating Pairwise correlation matrices...');
run PCMs_generate.m;

%% expected correlations (dataset-control)
disp('generating expected correlation values...');
run expCorr_generate.m;

Ngenes = length(intint);
disp(['# genes: ', num2str(Ngenes)]);
%% obsereved table
comps = nchoosek(F, 2);
obs = nan(Ngenes, length(comps));
comps_names = cell(length(comps), 1);
for i = 1:length(comps)
    spA = comps{i, 1};
    spB = comps{i, 2};
    nA = find(strcmp(F, spA));
    nB = find(strcmp(F, spB));
    obs(:, i) = corr_of_corr(nA, nB, ifFlt);
    comps_names{i} = [spA, '_', spB];
end
obsT = array2table(obs, 'variableNames', comps_names, 'rowNames', orfNames(intint));
disp('observed (comparative) table: obsT, is ready!');
%% expected table (dataset control)
exps = nan(Ngenes, length(F));
Ngenes_with_high_corr = nan(length(intint), 4);
for i = 1:length(F)
    exps(:, i) = nanmedian(expCorr.(F{i}).corrDiag, 2);
    Ngenes_with_high_corr(:, i) = sum(PCMs.(F{i}).PCMflt ~= 0, 2);
end
labels = cellfun(@(x) [x, '_Ngenes'], F, 'uniformoutput', false);
expT = array2table([exps, Ngenes_with_high_corr], 'variableNames', vertcat(F, labels), 'rowNames', orfNames(intint));
disp('expected (dataset-control) table: expT, is ready!');
%% nearest neighbors

nn_generate;

labels = cellfun(@(x) [x, '_geneIdx'], F, 'uniformoutput', false);
nnT = array2table([NN.corrs, NN.geneIdx], 'variableNames', vertcat(F, labels), 'rowNames', orfNames(intint));
disp('nearest neighbors table: nnT, is ready!');

labels = cellfun(@(x) [x, '_geneIdx'], crossSpNN.labels, 'uniformoutput', false);
nnCrossSpT = array2table([crossSpNN.maxCorrs, crossSpNN.geneIdx], 'variableNames',...
    vertcat(crossSpNN.labels, labels), 'rowNames', orfNames(intint));
disp('cross-species nearest neighbors table: nnCrossSpT, is ready!');

%% save PCM
name = datestr(datetime);
name = strrep(name(1:11), '-', '_');
ggc = struct;
ggc.PCMs = PCMs;
ggc.expCorr = expCorr;
ggc.intint = intint;
ggc.corrThrs = corrThrs;
ggc.ifFlt = ifFlt;
ggc.obsT = obsT;
ggc.expT = expT;
ggc.nnT = nnT;
ggc.nnCrossSpT = nnCrossSpT;
save(['ggc_', name, '.mat'], '-struct', 'ggc', '-v7.3');