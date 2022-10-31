function [intBasesLog, intBasesGeneIdx] = intBasesVec(genesIdx)
load('promoterLengthsORF.mat');
GP = load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
intBasesLog = false(12071326,1);
intBasesGeneIdx= zeros(12071326,1);

for g = 1:length(genesIdx)
    if ~isnan(promoterLengthsORF(genesIdx(g)))
        pos = GP.gene_infoR64.position(genesIdx(g),:);
        intBasesLog(GP.chrIdx(pos(1))+ pos(2)+[-promoterLengthsORF(genesIdx(g)):-1].*GP.gene_infoR64.dir(genesIdx(g))) = true;
        intBasesGeneIdx(GP.chrIdx(pos(1))+ pos(2)+[-promoterLengthsORF(genesIdx(g)):-1].*GP.gene_infoR64.dir(genesIdx(g))) = genesIdx(g);
    end
end