function [corrMat, idxVec,repeatDataMat,fieldIdx] = getRepeatsCorr(samples,varargin)
ip = inputParser;
ip.addParameter('dataType', 'sumProm')
ip.parse(varargin{:})

fileNames = dir('./indRepeats/TF ohnologs new Table*summary.xlsx');
keepCol = 'goodSamples';
for i = 1:numel(fileNames)
    temp = readtable([fileNames(i).folder,'/',fileNames(i).name]);
    goodRep{i} = temp.(keepCol);
end
goodRep = cat(1,goodRep{:});
goodRep = goodRep(~cellfun('isempty', goodRep));
goodRep = cellfun(@(x) strsplit(x,{',',' '}), goodRep,'UniformOutput',false);
goodRep = cat(2,goodRep{:});

if contains(ip.Results.dataType,'sumProm')
    sumPromRepeats = load('./indRepeats/sumPromRepeats.mat');
    FN = fieldnames(sumPromRepeats);
else
    mer7repeats = load('./indRepeats/mer7Repeats.mat');
    FN = fieldnames(mer7repeats);
end

for i = 1: numel(samples)
    if ~contains(samples{i}, {'_d','_lactis', '_DBD'})
        fieldIdx{i} = find(startsWith(upper(FN), upper(samples{i})) & ~contains(upper(FN), {'_D','_LACTIS','_DBD'}));
    else
        fieldIdx{i} = find(startsWith(upper(FN), upper(samples{i})));
    end
    fieldIdx{i} = fieldIdx{i}(ismember(FN(fieldIdx{i}), goodRep));
end
maxN = max(cellfun('prodofsize',fieldIdx));

for i = 1: numel(samples)
    if numel(fieldIdx{i}) == 0
        fieldIdx{i} = nan(maxN,1);
    else
        fieldIdx{i} = repmat(fieldIdx{i},ceil(maxN/numel(fieldIdx{i})),1);
        fieldIdx{i} = sort(fieldIdx{i}(1:maxN));
    end
    idxVec{i} = repmat(i, numel(fieldIdx{i}),1);
end
fieldIdx =  cat(1, fieldIdx{:});
idxVec = cat(1, idxVec{:});
if contains(ip.Results.dataType,'sumProm')
    for i =1:numel(fieldIdx)
        if ~isnan(fieldIdx(i))
            repeatDataMat(:,i) = sumPromRepeats.(FN{fieldIdx(i)});
        else
            repeatDataMat(:,i) = nan(6701,1);
        end
    end
else
    for i =1:numel(fieldIdx)
        if ~isnan(fieldIdx(i))
            repeatDataMat(:,i) = mer7repeats.(FN{fieldIdx(i)});
        else
            repeatDataMat(:,i) = nan(4^7/2,1);
        end
    end
end
corrMat = corr(repeatDataMat, 'rows','pairwise');


end