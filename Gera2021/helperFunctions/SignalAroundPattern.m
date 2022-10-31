function [signalMat,motifs,seqsMatch] =SignalAroundPattern(checStruct, TF, pattern, varargin)
ip = inputParser;
ip.addParameter('windowSize',150);
load('promoterIDXvec.mat', 'promoterIDXvec');
ip.addParameter('intBases',promoterIDXvec);
ip.addParameter('returnGenes',false);
ip.addParameter('intBasesGeneIdx',{});
ip.addParameter('normalize',false);

ip.parse(varargin{:});

intBasesGeneIdx =  ip.Results.intBasesGeneIdx;
intBases = ip.Results.intBases;

load('SC_genome.mat', 'SC_genome');
SC_genome = upper(cat(2,SC_genome(1:16).Sequence));
intBases = logical(intBases);
windowSize = ip.Results.windowSize;

[motifOccS,motifOccE, seqsMatch]  = regexp(SC_genome, pattern{1}, 'start','end','match');
motifOcc = round((motifOccS+motifOccE)/2);
goodMotifsIdx =  find(intBases(motifOcc));
seqsMatch = seqsMatch(goodMotifsIdx);
seqs = {};
for i = 1:length(goodMotifsIdx)
    seqs{i} = SC_genome(motifOccS(goodMotifsIdx(i)):motifOccE(goodMotifsIdx(i)));
    if ip.Results.returnGenes == true
        motifs.promoterIdx(i) = unique(intBasesGeneIdx(motifOccS(goodMotifsIdx(i)):motifOccE(goodMotifsIdx(i))));
    end  
end
if ~isempty(seqs)
    [a,b,c] =  unique(seqs);
    motifs.seq = a;
    motifs.motifIdx = c';
else
    motifs = {};
end
motifOcc = motifOcc(intBases(motifOcc));
fullProfile = chromosomes2fullProfile(checStruct, TF)';
if ip.Results.normalize
    fullProfile=fullProfile./max(movmean(fullProfile(promoterIDXvec==1),10));
end
if ~isempty(motifOcc)
    signalMat = fullProfile(motifOcc'+[-windowSize:windowSize]);
else
    signalMat = [];
end
end






