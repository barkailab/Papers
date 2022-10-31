function goodPos = createIntBasesForMotif()
GP = load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
load('promoterIDXvec.mat')

promBases = intBasesVec(find(~isnan(GP.gene_infoR64.felixTss(:,1))));

pos = GP.gene_infoR64.tamarTss;
pos(isnan(pos)) = GP.gene_infoR64.position(isnan(pos));
good = ~contains(GP.gene_infoR64.status, 'Dubious');
posAbs = nan(6701,2);
posAbs(~isnan(pos(:,1)),:) = GP.chrIdx(pos(~isnan(pos(:,1)),1))+pos(~isnan(pos(:,1)),[2,3]);
posAbsSort = sort(posAbs,2);
good = ~contains(GP.gene_infoR64.status, 'Dubious') & all(~isnan(posAbsSort),2);
genomeVec = zeros(sum(GP.chr_len),1);
for i = find(good)'
    genomeVec(posAbsSort(i,1):posAbsSort(i,2)) = 1;
end
% now is a vec where an orf is not count as part of promoters
goodPos = promBases > 0 & genomeVec(1:GP.chrIdx(17)) == 0 ;
end