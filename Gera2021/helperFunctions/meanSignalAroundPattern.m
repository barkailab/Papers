function [out,signalMatColl] = meanSignalAroundPattern(checStruct, TFs, pattern, varargin)
ip = inputParser;
ip.addParameter('windowSize',150);
load('promoterIDXvec.mat', 'promoterIDXvec');
ip.addParameter('intBases',promoterIDXvec);
ip.addParameter('normalize',false);
ip.parse(varargin{:});

intBases = ip.Results.intBases;
windowSize = ip.Results.windowSize;

for p = 1:numel(pattern)
    for tf = 1:numel(TFs)
        signalMat = SignalAroundPattern(checStruct, TFs(tf), pattern(p), 'intBases',intBases,'windowSize',windowSize,'normalize',ip.Results.normalize);
        if size(signalMat,1)>1
            out.aroundMotif{p}(tf,:) =  movmean(mean(signalMat),10);
        elseif size(signalMat,1)==1 & ~isempty(signalMat)
            out.aroundMotif{p}(tf,:) =  movmean(signalMat,10);
        else
            out.aroundMotif{p}(tf,:) =  nan(1, size(signalMat,2));
        end
        signalMatColl{p,tf}=signalMat;
    end
end


% 
% 
% load('SC_genome.mat', 'SC_genome')
% SC_genome = upper(cat(2,SC_genome(1:16).Sequence));
% intBases = logical(intBases);
% windowSize = ip.Results.windowSize;
% 
% for p = 1:numel(pattern)
%     [motifOccS,motifOccE]  = regexp(SC_genome,   pattern{p}, 'start','end');
%     motifOcc = round((motifOccS+motifOccE)/2);
%     motifOcc = motifOcc(intBases(motifOcc));
%     for tf = 1:numel(TFs)
%         fullProfile = chromosomes2fullProfile(checStruct, TFs(tf))';
%         signalMat = fullProfile(motifOcc'+[-windowSize:windowSize]);
%         if size(signalMat,1)>1
%             out.aroundMotif{p}(tf,:) =  movmean(mean(signalMat),10);
%         elseif size(signalMat,1)<1
%             out.aroundMotif{p}(tf,:) =  movmean(signalMat,10);
%         else
%             out.aroundMotif{p}(tf,:) =  nan(1, size(signalMat,2));
%         end
%     end
% end
end


