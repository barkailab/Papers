function [expCorr] = expected_corr_general_fixedN(data, Nshuff, geneThrs, ...
    expThrs, corrThrs, NCorrThrs, zeroDiag, ifFlt, Nexpected)
global NgenesToOL corrThList len_allGenes
%% expected correlation
% Is the gene-gene correlation diverged or just noisy?
% calculate expected correlation between two halves of the data
% split the data into two and shuffle to get an expected distribution

corrDiag = cell(Nshuff, 1);
genes = cell(Nshuff, 1);

PCMhalf1 = cell(Nshuff, 1);
PCMhalf2 = cell(Nshuff, 1);
OL = cell(Nshuff, 1);
corrByTh = cell(Nshuff, 1);
%%
% parfor i = 1:Nshuff
for i = 1:Nshuff
    Nsamples = size(data, 2);
    sampleIdx = 1:Nsamples;
    if ceil(Nsamples/2) < Nexpected
        randHalf1 = randsample(sampleIdx, floor(Nsamples/2));
        randHalf2 = find(~ismember(sampleIdx, randHalf1));
    else
        randHalf1 = randsample(sampleIdx, Nexpected);
        randHalf2 = find(~ismember(sampleIdx, randHalf1));
        randHalf2 = randsample(randHalf2, Nexpected);
    end
    
    dataHalf1 = data(:, randHalf1);
    dataHalf2 = data(:, randHalf2);
    
    PCMhalf = struct;
    
    PCMhalf.rep1 = gene_gene_corr(dataHalf1, geneThrs, expThrs, corrThrs, NCorrThrs, zeroDiag);
    PCMhalf.rep2 = gene_gene_corr(dataHalf2, geneThrs, expThrs, corrThrs, NCorrThrs, zeroDiag);

    currintG = intersect(PCMhalf.rep1.genes, PCMhalf.rep2.genes);
    for j = 1:2
        currHalf = PCMhalf.(['rep', num2str(j)]);
        idx = find(ismember(currHalf.genes, currintG));
        currHalf.data = currHalf.data(:, idx);
        currHalf.PCM = currHalf.PCM(idx, idx);
        currHalf.PCMflt = currHalf.PCMflt(idx, idx);
        currHalf.pval = currHalf.pval(idx, idx);
        currHalf.genes = currHalf.genes(idx);
        PCMhalf.(['rep', num2str(j)]) = currHalf;
    end
    if ifFlt
        a = PCMhalf.rep1.PCMflt; b = PCMhalf.rep2.PCMflt;
    else
        a = PCMhalf.rep1.PCM; b = PCMhalf.rep2.PCM;
    end
    CCM = corr(a, b, 'rows', 'pairwise');
    cordiagRep = diag(CCM);
    
    corrDiag{i} = cordiagRep;
    genes{i} = currintG;
    PCMhalf1{i} = PCMhalf.rep1.PCM;
    PCMhalf2{i} = PCMhalf.rep2.PCM;
end

intAll = 1:len_allGenes;
for i = 1:Nshuff
    intAll = intersect(intAll, genes{i});
end

for i = 1:Nshuff
    idx = find(ismember(genes{i}, intAll));
    genes{i} = genes{i}(idx);
    corrDiag{i} = corrDiag{i}(idx);
    PCMhalf1{i} = PCMhalf1{i}(idx, idx);
    PCMhalf2{i} = PCMhalf2{i}(idx, idx);
end
expCorr = struct;
expCorr.genes = genes{1};
expCorr.corrDiag = cell2mat(corrDiag');
expCorr.PCMhalves = {PCMhalf1{1}, PCMhalf2{1}};
end