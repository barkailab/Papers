function [mean_of_corr, std_of_corr] = corrBetweenRepeatsOf2Samples(samples,varargin)
ip = inputParser;
ip.addParameter('dataType', 'sumProm')
ip.parse(varargin{:})

fileNames = dir('/home/labs/barkailab/tamarge/Master/TF ohnologs new Table*summary.xlsx');
keepCol = 'goodSamples';
for i = 1:numel(fileNames)
    temp = readtable([fileNames(i).folder,'/',fileNames(i).name]);
    goodRep{i} = temp.(keepCol);
end
goodRep = cat(1,goodRep{:});
goodRep = goodRep(~cellfun('isempty', goodRep));
goodRep = cellfun(@(x) strsplit(x,{',',' '}), goodRep,'UniformOutput',false);
gododRep = cat(2,goodRep{:});

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
    %fieldIdx{i} = fieldIdx{i}(ismember(FN(fieldIdx{i}), goodRep));
end
maxN = max(cellfun('prodofsize',fieldIdx));

if ~(isempty(fieldIdx{1}) | isempty(fieldIdx{2})) % both are not empty
    sprintf('%s \n%s', strjoin(FN(fieldIdx{1}), ', '), strjoin(FN(fieldIdx{2}), ', '))
    c=1;
    if contains(ip.Results.dataType,'sumProm')
        for r1 = fieldIdx{1}'
            for r2 =  fieldIdx{2}'
                curr_corr(c) = corr(sumPromRepeats.(FN{r1}), sumPromRepeats.(FN{r2}), 'rows','pairwise');
                c=c+1;
            end
        end
    else
        for r1 = fieldIdx{1}'
            for r2 =  fieldIdx{2}'
                curr_corr(c) = corr(mer7repeats.(FN{r1}), mer7repeats.(FN{r2}), 'rows','pairwise');
                c=c+1;
            end
        end
    end
    
elseif isempty(fieldIdx{1}) & isempty(fieldIdx{2}) % both are empty
    c=1;
    if ~contains(samples{i}, {'_d','_lactis', '_DBD'})
        load('summaryTable.mat')
        row_idx = find((contains(summaryTable.p1, samples{1}) |  contains(summaryTable.p1, samples{2})) & (contains(summaryTable.p2, samples{1}) |  contains(summaryTable.p2, samples{2})))
        if contains(ip.Results.dataType,'sumProm')
            curr_corr = summaryTable.sumPromCorr(row_idx);
        else
            curr_corr = summaryTable.WTs7merCorr(row_idx);
        end
    end
    
else % only one is empty
    load('./indRepeats/sumPromNormSamplesWithoutRepeats.mat')
    sampleWithoutRepeats = samples{cellfun(@length, fieldIdx) == 0};
    if contains(ip.Results.dataType,'sumProm')
        c=1;
        for r2 = fieldIdx{find(cellfun(@length, fieldIdx) ~= 0)}'
            curr_corr(c) = corr(sumPromNormSamplesWithoutRepeats.sumProm.(sampleWithoutRepeats),...
                sumPromRepeats.(FN{r2}), 'rows','pairwise');
            c=c+1;
        end
    else
        c=1;
        goodPos = createIntBasesForMotif();
        currNorm = chromosomes2fullProfile(sumPromNormSamplesWithoutRepeats, {sampleWithoutRepeats});
        currMer = mer_occupancy(currNorm,7,'intBases', goodPos,'method','normal');
        for r2 = fieldIdx{find(cellfun(@length, fieldIdx) ~= 0)}'
            curr_corr(c) = corr(currMer.score(:,1), mer7repeats.(FN{r1}), 'rows','pairwise');
            c=c+1;
        end
    end
end

curr_corr
mean_of_corr = mean(curr_corr);
std_of_corr = std(curr_corr);

end