

[status, sheetNames] = xlsfinfo(exp_per_sp_table);
numSheets = length(sheetNames);
dataS = struct;
dataLabels = struct;
for i = 1:numSheets
    t = readtable(exp_per_sp_table, 'sheet', sheetNames{i});
    orfNames = t.Row;
    sampleNames = t.Properties.VariableNames(2:end);
    dataS.(sheetNames{i}) = table2array(t(:, 2:end));
    dataLabels.(sheetNames{i}) = sampleNames;
end