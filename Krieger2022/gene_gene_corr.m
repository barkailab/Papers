function [PCMs] = gene_gene_corr(data, geneThrs, expThrs, corrThrs, NCorrThrs, zeroDiag)
% INPUTS:
% data = expression data, samples are columns and genes are rows
% expThrs = take expression vaules > expThrs
% geneThrs = take only genes that have expression > expThrs in X % of the
% samples (X = geneThrs)
% corrThrs = filter the pairwise correlation matrix
% NCorrThrs = filter by number of genes with significant correlation.
% default is no filter
% OUTPUT:
% data = expression data, where samples are rows and genes are columns
% PCM = pairwise colrrelation matrix
% PCMflt = filtered PCM
% genes = gene list after filtering. corresponds to the genes that comprise
% the PCM
% pval = p-value for correlation coefficients

if nargin < 6
    geneThrs = 10;
    expThrs = min(min(data));
    corrThrs = 0.3;
    NCorrThrs = 0;
    zeroDiag = 0;
end
%%
highExp = sum(abs(data) > expThrs, 2);
notnanPerGene = (highExp./size(data, 2))*100;
idxG = find(notnanPerGene >= geneThrs);
data = data(idxG, :);
data(isnan(data)) = expThrs;
data = data';

% generate PCM
[PCM, pval] = corr(data, 'rows', 'pairwise');

if zeroDiag
    PCM(1:(size(PCM, 1)+1):numel(PCM)) = 0;
end

% filter correlation coefficient 
PCMflt = PCM;
PCMflt(abs(PCMflt) < corrThrs) = 0;
% filter out genes with low number of correlating genes
NgenesCorr = sum(PCMflt > 0, 2);
badGenes = find(NgenesCorr < NCorrThrs);
if ~isempty(badGenes)
    idxG(badGenes) = [];
    PCMflt(badGenes, :) = [];
    PCMflt(:, badGenes) = [];
    PCM(badGenes, :) = [];
    PCM(:, badGenes) = [];
    data(:, badGenes) = [];
    pval(badGenes, :) = [];
    pval(:, badGenes) = [];
end
%%
% save in PCMs struct
PCMs = struct;
PCMs.data = data; 
PCMs.PCM = PCM; 
PCMs.PCMflt = PCMflt; 
PCMs.genes = idxG; 
PCMs.pval = pval;

end
