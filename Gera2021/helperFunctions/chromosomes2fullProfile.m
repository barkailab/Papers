function [allSamplesProfile] = chromosomes2fullProfile(checStruct, sampleName)

allSamplesProfile = [];
for j = 1: numel(sampleName)
    fullProfile = [];
    for i = 1: size(checStruct.norm.(sampleName{j}),2)
        if size(checStruct.norm.(sampleName{j}){1,i},2) == 1
            fullProfile = [fullProfile, checStruct.norm.(sampleName{j}){1,i}'];
        else
            fullProfile = [fullProfile, checStruct.norm.(sampleName{j}){1,i}];
        end
    end
    allSamplesProfile = [allSamplesProfile, fullProfile'];
end
end