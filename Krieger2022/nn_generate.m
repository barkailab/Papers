% generate 2nd order pairwise corerlation matrix
PCM2 = nan(Ngenes, Ngenes, length(F));
NN = struct;
NN.corrs = nan(Ngenes, length(F));
NN.geneIdx = nan(Ngenes, length(F));
% number of nearest neighbors: for cross-species comparison
nnum = 5;
NN.corrs_nnum = nan(Ngenes, nnum, length(F));
NN.geneIdx_nnum = nan(Ngenes, nnum, length(F));
for i = 1:length(F)
    currPCM = PCMs.(F{i}).PCM;
    disp(['generating second PCM for ', F{i}, '...']);
    currPCM2 = corr(currPCM, 'rows', 'pairwise');
    PCM2(:, :, i) = currPCM2;
    currPCM2(1:(Ngenes+1):end) = nan;
    
    [currMax2, maxIdx2] = max(currPCM2, [], 2);
    NN.corrs(:, i) = currMax2;
    NN.geneIdx(:, i) = maxIdx2;
    
    for j = 1:Ngenes
        currVector = currPCM2(:, j);
        [sortdVector, sortdidx] = sort(currVector, 'descend');
        NN.corrs_nnum(j, :, i) = sortdVector(2:nnum+1);
        NN.geneIdx_nnum(j, :, i) = sortdidx(2:nnum+1);   
    end 
end

%% cross-species NN
% get all combinations where order matters
allComps = {};
allComps_labels = {};
n=1;
for i = 1:length(F)
    for j = 1:length(F)
        if i==j; continue; end
        allComps{n,1} = F{i};
        allComps{n,2} = F{j};
        allComps_labels{n} = [F{i}, ' by ', F{j}];
        n=n+1;
    end
end
%%
crossSpNN = struct;
crossSpNN.maxCorrs = nan(Ngenes, size(allComps, 1));
crossSpNN.geneIdx = nan(Ngenes, size(allComps, 1));
for i = 1:size(allComps, 1)
    spA = allComps{i, 1};
    spB = allComps{i, 2};
    nA = find(strcmp(F, spA));
    nB = find(strcmp(F, spB));
    for j = 1:Ngenes
        geneIdxB = NN.geneIdx_nnum(j, :, nB);
        corrInA = PCM2(j, geneIdxB, nA);
        [maxCorr, maxIdx] = max(corrInA);
        crossSpNN.maxCorrs(j, i) = maxCorr;
        crossSpNN.geneIdx(j, i) = geneIdxB(maxIdx);
    end
end
crossSpNN.labels = allComps_labels;
        
